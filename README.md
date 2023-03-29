# SeqWell assembly workflow

This is a port of the SeqWell assembly workflow for plasmid assembly to LatchBio.


It performs assembly of raw short-reads and identify circularized sequences (plasmids)

----

## Inputs

- Paired-ends Fastq reads

## Outputs

- Summary CSV with assembly statistics
- Sequences in Fasta format
- Assembly graphs (GFA) 
- Coverage plots (PNG)
- Intermediate files for read aligns

## Workflow

1. Down sample reads if required
2. Preprocess raw reads with BBmap (adapter trimming, merge and clean)
3. Aseembly with Uniclycler
4. Circularize sequences
5. Align reads to assembled sequences
6. Generate statistics
