Setup packages
--------------

Check to see if the `BiocManager` and `polyester` packages are loadable, or install them if not.

``` r
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
```

    ## Loading required namespace: BiocManager

``` r
if (!requireNamespace("polyester"))
  BiocManager::install("polyester")
```

    ## Loading required namespace: polyester

``` r
if (!requireNamespace("edgeR"))
  BiocManager::install("edgeR")
```

    ## Loading required namespace: edgeR

Setup a small experiment
------------------------

### `simulate_experiment` example

A FASTA file called `chr22.fa` is provided with `polyester`. This file contains sequences for 918 transcripts on chromosome 22, as annotated in hg19. For this very small example, we will only simulate from the first 20 of these transcripts.

We will set the first 2 transcripts to be overexpressed in group A and the next 2 transcripts to be overexpressed in group B, each at a fold change of 3. The way to do this in Polyester is to provide a "fold change matrix": for each transcript and each group, specify a fold change. Polyester will generate baseline read numbers (assuming no differential expression), and will then multiply those mean numbers by the fold change you specify for the replicates in that group. The fold change matrix for this simple 2-group experiment looks like this:

``` r
fold_changes = matrix(c(4,4,rep(1,18),1,1,4,4,rep(1,16)), nrow=20)
head(fold_changes)
```

    ##      [,1] [,2]
    ## [1,]    4    1
    ## [2,]    4    1
    ## [3,]    1    4
    ## [4,]    1    4
    ## [5,]    1    1
    ## [6,]    1    1

The matrix has two columns, since there will be two groups (cases and controls) in this experiment.

The rest of the experiment can be simulated with code like the chunk below.

``` r
library(polyester)
library(Biostrings)
# FASTA annotation
fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet(fasta_file)
# subset the FASTA file to first 20 transcripts
small_fasta = fasta[1:20]
writeXStringSet(small_fasta, 'chr22_small.fa')
# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(20 * width(small_fasta) / 100)
# simulation call:
simulate_experiment('chr22_small.fa', reads_per_transcript=readspertx, 
    num_reps=c(10,10), fold_changes=fold_changes, outdir='simulated_reads') 
```

Process the reads
-----------------

After simulating the reads, the julia code builds a dictionary of counts per file for each kmer, and kmers with more than 1 read per file (on average) are stored in a table.

``` bash
if [ ! -f hashcounts.tsv ]; then
  julia diffhash.jl
  julia showhash.jl > hashcounts.tsv
fi
```

Find kmers that are differentially expressed
--------------------------------------------

``` r
sim_rep_info <- read.delim("simulated_reads/sim_rep_info.txt")
hashcounts <- read.delim("hashcounts.tsv", header=FALSE, row.names=1)
```

I should build the design from `sim_rep_info` but I'm lazy and impatient.

``` r
design <- cbind(rep(1, 20), c(rep(0,10), rep(1,10)))
colnames(design) <- c("C", "CvsT")
```

See the limma user's guide for examples of analyzing rnaseq data.

``` r
library(edgeR)
```

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

``` r
dge <- DGEList(counts=hashcounts)
logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))
```

    ##                   logFC  AveExpr         t      P.Value    adj.P.Val
    ## CATGAACTAAAAG -2.144087 6.850431 -39.53762 5.727023e-51 6.970428e-47
    ## TACCATGAACTAA -2.169818 6.836911 -39.33501 8.190064e-51 6.970428e-47
    ## CCTAGTCTTGTCA -2.163425 6.824691 -39.23303 9.811922e-51 6.970428e-47
    ## CCATGAACTAAAA -2.147536 6.848546 -39.20808 1.025618e-50 6.970428e-47
    ## ACCATGAACTAAA -2.155501 6.842216 -39.03735 1.389480e-50 6.970428e-47
    ## CTACCATGAACTA -2.168371 6.842266 -38.92439 1.699780e-50 6.970428e-47
    ## CATGAAGACAGTA -2.209651 6.793377 -38.92409 1.700677e-50 6.970428e-47
    ## GACAAGACTAGGA -2.162630 6.829611 -38.92298 1.704072e-50 6.970428e-47
    ## CTAGTCTTGTCAG -2.150109 6.826680 -38.92192 1.707276e-50 6.970428e-47
    ## ATACTGTCTTCAT -2.193787 6.797756 -38.86242 1.898957e-50 6.977718e-47
    ##                      B
    ## CATGAACTAAAAG 106.0530
    ## TACCATGAACTAA 105.6987
    ## CCTAGTCTTGTCA 105.5197
    ## CCATGAACTAAAA 105.4758
    ## ACCATGAACTAAA 105.1749
    ## CTACCATGAACTA 104.9752
    ## CATGAAGACAGTA 104.9747
    ## GACAAGACTAGGA 104.9727
    ## CTAGTCTTGTCAG 104.9708
    ## ATACTGTCTTCAT 104.8654

We can find which k-mers are significantly differencially expressed, using `fdr` to correct for multiple testing.

``` r
testresults <- decideTests(fit[,2]$p.value, adjust.method = "fdr")
sum(testresults != 0)
```

    ## [1] 5131

Writing out the kmers to a file will allow us to filter the reads to find reads that contain these kmers.

``` r
diffkmers <- rownames(fit)[testresults != 0]
write(diffkmers, "diffkmers.txt")
```
