The RiBoSor is a tool to create alternative reading frames within an existing gene.

# Principle

Creation of a new reading frame is achieved using codon redundancy in the genetic code, so the function of the original gene is preserved.

Such alternative reading frame permits to initiate the translation of two proteins in the same sequence.

# Usage

`./ribosor.py input.fasta` will treat separately all gene sequences present in the input fasta file. The program will create a directory `input`, and one sub-directory per input sequence.
