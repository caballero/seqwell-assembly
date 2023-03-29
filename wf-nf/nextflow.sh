#!/bin/bash

run="20210917_Example"
plates="SeqWell-07SEP21_SO10828_FASTQ"
downsample=10000

nextflow run \
 wf-nf/fastq_assembly.nf \
 -work-dir work \
 --plates $plates \
 --downsample $downsample \
 --run $run \
 --readsdir "data/$run/$plates" \
 -with-docker seqwell:latest
