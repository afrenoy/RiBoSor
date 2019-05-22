The RiBoSor is a tool to create overlapping reading frames within an existing gene.

# Principle

Creation of a new reading frame is achieved using codon redundancy in the genetic code, so the function of the original gene is preserved.

Such overlapping reading frame permits to initiate the translation of two proteins in the same sequence.

# Usage

## Basic

`./ribosor.py input.fasta` will treat separately all gene sequences present in the input fasta file. The program will create a directory `input`, and one sub-directory per input sequence.

## Advanced

Full documentation is available with `./ribosor.py -h`

