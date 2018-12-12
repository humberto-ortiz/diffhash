# diffhash - find regions of differential expression by looking at kmers

Humberto Ortiz-Zuazaga
<humberto.ortiz@upr.edu>

# Introduction

Finding differentially expressed genes in RNASeq from non-model organisms
is a complex task. Usually a transcriptome is assembled from the reads, then
the reads are aligned to the assembled transcriptome, then the counts are
tested for differential expression.

Much of the reads are not differentially expressed, and assembly of non-model
organisms is not trivial.

Can we pull regions of the assembly graph that are differentially
expressed *before* assembling? Will it simplify the assembly?

This project seeks to answer these questions.

# Methods

`diffhash.Rmd` R markdown code, uses the `polyester` package to generate
simulated reads for a few genes, as demonstrated in the `polyester` vignette. Then calls the julia code below to find the occurrence of each kmer. With the count data, uses `edgeR` to find differentially expressed kmers and write them out to `diffkmers.txt`.

`diffhash.jl` Julia code to count the occurrence of each kmer in a collection
of sequence files. Persists the counts in a `JLD` file. Uses `BioSequences` to
read FASTA files and iterate over kmers.

`showhash.jl` Julia code to print the stored kmer counts.

# Results

Before removing kmers with fewer than 1 occurrence on average.

```
$ wc -l hashcounts.tsv
  964512 hashcounts.tsv
```

After filtering kmers

```
$ wc -l hashcounts.tsv
   36745 hashcounts.tsv
```

The result of running the .Rmd file is shown in [diffhash.md](diffhash.md).

## Filtering reads

We can include a read if any of it's kmers are differentially expressed or if
all of it's kmers are differentially expressed:

any:
```
$ du -sh simulated_reads any
106M	simulated_reads
 48M	any
```

all:
```
$ du -sh simulated_reads filtered
106M	simulated_reads
 15M	filtered
```

## Differentially expressed kmers

Using "any" results in more false positives at this stage. All 20 of the
simulated genes have at least one read with a differentially expressed kmer.


Real differentially expressed genes:

```
$ grep "^>" ../*.fa | cut -d\| -f 2 | head -n 4
424037187
424037186
209977002
52546690
```

Filtered set of genes (out of 20):

```
$ grep "read" *.fa | cut -d\| -f 2 | sort -u
209977002
213512645
29029549
29029551
424037186
424037187
52546690
```

That's 4 true positives and 3 false positives (out of 20 genes).

## "Assembling"

Use `bcalm` to build unitigs.

```
$ ls -1 filtered/sample_* > list_reads
$ bcalm -in list_reads -kmer-size 21 -abundance-min 2
```

produces a bunch of  `list_reads*` files.

```
$ ~/src/spacegraphcats/bcalm/scripts/convertToGFA.py list_reads.unitigs.fa test.gfa 21
```
Bandage can read GFA and blast fragments against `chr22_small.fa`.

The 4 differentially expressed genes are the three structures on the left of the
figure below, and the two small fragments on the right are small pieces of false
positive genes.

![Bandage graph of fragments](test.png)

# References

1. Frazee AC, Jaffe AE, Kirchner R, Leek JT (2018). polyester: Simulate RNA-seq
reads. R package version 1.18.0.

1.   Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor
  package for differential expression analysis of digital gene expression
  data. Bioinformatics 26, 139-140

1.  McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression analysis
  of multifactor RNA-Seq experiments with respect to biological variation.
  Nucleic Acids Research 40, 4288-4297

1. Chikhi, Rayan, Antoine Limasset, and Paul Medvedev. "Compacting de Bruijn
   graphs from sequencing data quickly and in low memory." Bioinformatics 32.12
   (2016): i201-i208.

1. Wick R.R., Schultz M.B., Zobel J. & Holt K.E. (2015). Bandage: interactive
   visualisation of de novo genome assemblies. Bioinformatics, 31(20), 3350-3352.
