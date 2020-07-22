#!/home/humberto/humberto/miniconda3/envs/roble3/bin/python

import sys
import screed
import os
import khmer
import textwrap
from khmer import khmer_args, Countgraph
from contextlib import contextmanager
from khmer.khmer_args import (build_counting_args, add_loadgraph_args,
                              report_on_config, calculate_graphsize,
                              sanitize_help, check_argument_range)
from khmer.khmer_args import FileType as khFileType
import argparse
from khmer.kfile import (check_space, check_space_for_graph,
                         check_valid_file_exists, add_output_compression_type,
                         get_file_writer, describe_file_handle)
from khmer.utils import (write_record, broken_paired_reader, ReadBundle,
                         clean_input_reads)
from khmer.khmer_logger import (configure_logging, log_info, log_error)
from scipy.stats import ttest_ind_from_stats

def get_parser():
    epilog = """\
    Discard sequences based on whether or not their median k-mer abundance lies
    above a specified cutoff. Kept sequences will be placed in <fileN>.keep.

    By default, paired end reads will be considered together; if
    either read should be kept, both will be kept. (This keeps both
    reads from a fragment, and helps with retention of repeats.)
    Unpaired reads are treated individually.

    If :option:`-p`/:option:`--paired` is set, then proper pairing is required
    and the script will exit on unpaired reads, although
    :option:`--unpaired-reads` can be used to supply a file of orphan
    reads to be read after the paired reads.

    :option:`--force_single` will ignore all pairing information and treat
    reads individually.

    With :option:`-s`/:option:`--savegraph`, the k-mer countgraph
    will be saved to the specified file after all sequences have been
    processed. :option:`-l`/:option:`--loadgraph` will load the
    specified k-mer countgraph before processing the specified
    files.  Note that these graphs are are in the same format as those
    produced by :program:`load-into-counting.py` and consumed by
    :program:`abundance-dist.py`.

    To append reads to an output file (rather than overwriting it), send output
    to STDOUT with `--output -` and use UNIX file redirection syntax (`>>`) to
    append to the file.

    Example::

        normalize-by-median.py -k 17 tests/test-data/test-abund-read-2.fa

    Example::

        normalize-by-median.py -p -k 17 \\
        tests/test-data/test-abund-read-paired.fa

    Example::

        normalize-by-median.py -p -k 17 -o - tests/test-data/paired.fq \\
        >> appended-output.fq

    Example::

        normalize-by-median.py -k 17 -f tests/test-data/test-error-reads.fq \\
        tests/test-data/test-fastq-reads.fq

    Example::

        normalize-by-median.py -k 17 -s test.ct \\
        tests/test-data/test-abund-read-2.fa \\
        tests/test-data/test-fastq-reads.fq"""
    parser = build_counting_args(
        descr="Do digital normalization (remove mostly redundant sequences)",
        epilog=textwrap.dedent(epilog),
        citations=['diginorm'])
    parser.add_argument('-q', '--quiet', dest='quiet', default=False,
                        action='store_true')
    parser.add_argument('-p', '--paired', action='store_true',
                        help='require that all sequences be properly paired')
    parser.add_argument('--force_single', dest='force_single',
                        action='store_true',
                        help='treat all sequences as single-ended/unpaired')
    parser.add_argument('-u', '--unpaired-reads',
                        metavar="unpaired_reads_filename",
                        help='include a file of unpaired reads to which '
                        '-p/--paired does not apply.')
    parser.add_argument('-s', '--savegraph', metavar="filename", default=None,
                        help='save the k-mer countgraph to disk after all '
                        'reads are loaded.')
    parser.add_argument('-R', '--report',
                        help='write progress report to report_filename',
                        metavar='report_filename', type=argparse.FileType('w'))
    parser.add_argument('--report-frequency',
                        metavar='report_frequency', type=int, default=100000,
                        help='report progress every report_frequency reads')
    parser.add_argument('-f', '--force', dest='force',
                        help='continue past file reading errors',
                        action='store_true')
    parser.add_argument('-o', '--output', metavar="filename",
                        type=khFileType('wb'),
                        default=None, dest='single_output_file',
                        help='only output a single file with '
                        'the specified filename; use a single dash "-" to '
                        'specify that output should go to STDOUT (the '
                        'terminal)')
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        help='Input FAST[AQ] sequence filename.', nargs='+')
    add_loadgraph_args(parser)
    parser.add_argument('-z', '--loadgraph2', metavar="filename", default=None, help='load a second k-mer graph')
    add_output_compression_type(parser)
    return parser


def main():  # pylint: disable=too-many-branches,too-many-statements
    parser = sanitize_help(get_parser())
    args = parser.parse_args()

    configure_logging(args.quiet)
    report_on_config(args)

    report_fp = args.report
    force_single = args.force_single

    # check for similar filenames
    # if we're using a single output file only check for identical filenames
    # otherwise, check for identical BASE names as well.
    filenames = []
    basenames = []
    for pathfilename in args.input_filenames:
        filenames.append(pathfilename)
        if args.single_output_file:
            continue  # nothing more to worry about

        basename = os.path.basename(pathfilename)
        if basename in basenames:
            log_error('ERROR: Duplicate filename--Cannot handle this!')
            log_error('** Exiting!')
            sys.exit(1)

        basenames.append(basename)

    # check that files exist and there is sufficient output disk space.
    check_valid_file_exists(args.input_filenames)
    check_space(args.input_filenames, args.force)
    if args.savegraph is not None:
        graphsize = calculate_graphsize(args, 'countgraph')
        check_space_for_graph(args.savegraph, graphsize, args.force)

    # load or create counting table.
    if args.loadgraph:
        log_info('loading k-mer countgraph from {graph}',
                 graph=args.loadgraph)
        countgraph1 = Countgraph.load(args.loadgraph)

    # load second counting table.
    if args.loadgraph2:
        log_info('loading k-mer countgraph from {graph}',
                 graph=args.loadgraph2)
        countgraph2 = Countgraph.load(args.loadgraph2)

    # make a list of all filenames and if they're paired or not;
    # if we don't know if they're paired, default to allowing but not
    # forcing pairing.
    files = []
    for element in filenames:
        files.append([element, args.paired])
    if args.unpaired_reads:
        files.append([args.unpaired_reads, False])

    #
    # main loop: iterate over all files given, do diginorm.
    #

    for filename, require_paired in files:
        if not args.single_output_file:
            output_name = os.path.basename(filename) + '.keep'
            outfp = open(output_name, 'wb')
            outfp = get_file_writer(outfp, args.gzip, args.bzip)

        screed_iter = clean_input_reads(screed.open(filename))
        reader = broken_paired_reader(screed_iter, min_length=args.ksize,
                                      force_single=force_single,
                                      require_paired=require_paired)

        # actually do diginorm
        for _, is_paired, read0, read1 in reader:
            for record in snarf(is_paired, read0, read1, countgraph1, countgraph2):
                if record is not None:
                    write_record(record, outfp)

def snarf(is_paired, read0, read1, countgraph1, countgraph2):
        batch = ReadBundle(read0, read1)

        # if any in batch have differential coverage, consume &yield
        #coverage1 = batch.coverages(countgraph1)
        #coverage2 = batch.coverages(countgraph2)

        #aread = list(batch.reads)[0]
        #print(countgraph1.get_median_count(aread.cleaned_seq))
        for record in batch.reads:
            median1, mean1, sd1 = countgraph1.get_median_count(record.cleaned_seq)
            median2, mean2, sd2 = countgraph2.get_median_count(record.cleaned_seq)
            nobs1 = nobs2 = len(record.cleaned_seq) - 31
            res = (ttest_ind_from_stats(mean1, sd1, nobs1, 
                                        mean2, sd2, nobs2))
            if res.pvalue < 0.05:
            #self.countgraph.consume(record.cleaned_seq)
                yield record


# Main
if __name__ == '__main__':
    main()
