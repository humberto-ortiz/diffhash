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

What if we can pull regions of the assembly graph that are differentially
expressed *before* assembling? Will it simplify the assembly?

This project seeks to answer these questions.

# Methods

`diffhash.Rmd` R markdown code, uses the `polyester` package to generate
simulated reads for a few genes.

`diffhash.jl` Julia code to count the occurrence of each kmer in a collection
of sequence files. Persists the counts in a `JLD` file. Uses `BioSequences` to
read FASTA files and iterate over kmers.

`showhash.jl` Julia code to print the stored kmer counts.
