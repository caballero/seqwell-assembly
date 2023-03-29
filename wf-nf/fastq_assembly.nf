#!/usr/bin/env nextflow

params.dev=false
params.readsdir="data/"
params.number_of_inputs = 48
params.kmers="27,47,63,77,89,99,107,115,121,127"
params.plates=""
params.run=""
params.analysis="assembly"
params.downsample=0
work_dir = file(workflow.workDir).toString()
bin_dir = "/opt/conda/envs/seqWell-nf/bin"

// get Fastq to process
fq_ch = Channel
     .fromFilePairs("${params.readsdir}/*_R{1,2}*.fastq.gz")
     .map{ pair -> tuple(pair[0], pair[0].substring(0,pair[0].length()-4), pair[1][0], pair[1][1])}
     .take( params.dev ? params.number_of_inputs : -1 )


// downsample fastq
process downsample {

     input:
     tuple val(pair_id), val(plate_id), path(read1), path(read2) from fq_ch

     output:
     tuple val(pair_id), val(plate_id), path('*_R1*.fastq.gz'), path('*_R2*.fastq.gz') into sample_fq

     """
     if [ $params.downsample -gt 0 ]; then
          seqtk sample -s 14 ${read1} $params.downsample | gzip > ${pair_id}.${params.downsample}_R1.fastq.gz
          seqtk sample -s 14 ${read2} $params.downsample | gzip > ${pair_id}.${params.downsample}_R2.fastq.gz
     else
          ln -s ${read1} ${pair_id}_full_R1.fastq.gz
          ln -s ${read2} ${pair_id}_full_R2.fastq.gz
     fi
     """
}


process bbmerge {

    input:
    tuple val(pair_id), val(plate_id), path(read1), path(read2) from sample_fq

    output:
    path("*")
    tuple val(pair_id), val(plate_id), path("*.norm.merged.fq") into LR_ch
    path("*.norm75.fq") into NORM_ch

    publishDir path: "$work_dir/${params.run}/${params.analysis}/${plate_id}_LR/", pattern: "*.norm.merged.fq", mode: 'copy'
    publishDir path: "$work_dir/${params.run}/${params.analysis}/${plate_id}_NORM/", pattern: "*.norm75.fq", mode: 'copy'


    """
    $bin_dir/bbduk.sh -Xmx2G \
    in=$read1 \
    in2=$read2 \
    out=${pair_id}.clean.fastq.gz \
    ref=/opt/conda/envs/seqWell-nf/opt/bbmap-38.90-3/resources/nextera.fa.gz \
    k=21 mink=15 ktrim=n hdist=2 tpe tbo \
    overwrite=true

    $bin_dir/bbmerge-auto.sh -Xmx2G \
    in=${pair_id}.clean.fastq.gz \
    out=${pair_id}.merged.fq \
    extend2=20 iterations=10 \
    k=60 ecc=true \
    ecctadpole=true \
    reassemble=true \
    rem=true \
    merge=true \
    strict=true

    $bin_dir/bbnorm.sh -Xmx2G \
    in=${pair_id}.merged.fq \
    target=100 \
    maxdepth=150 \
    fixspikes=t \
    overwrite=true \
    out=${pair_id}.norm.merged.fq

    $bin_dir/tadpole.sh -Xmx2G \
    in=${pair_id}.clean.fastq.gz \
    out=${pair_id}.clean.ext.fastq.gz \
    mode=extend \
    extendleft=100 \
    extendright=100 \
    ecc=t \
    overwrite=true

    $bin_dir/bbnorm.sh -Xmx2G \
    in=${pair_id}.clean.ext.fastq.gz \
    fixspikes=t \
    target=75 \
    maxdepth=100 \
    passes=2 \
    fixspikes=true \
    tossbadreads=true \
    prefilter=true \
    ecc=true \
    overwrite=true \
    out=${pair_id}.norm75.fq

    """
 
}

(fq_align, fq_cycle) = LR_ch.into(2)

process unicycler {
     
     publishDir path: "$work_dir/${params.run}/${params.analysis}/${plate_id}_Unicycler", pattern: '*.gfa', mode: 'copy'
 
     input:
     tuple val(pair_id), val(plate_id), path(norm_fq) from fq_cycle

     output:
     tuple val(pair_id), val(plate_id), path('*.gfa') into unicycler_gfa

     """
     if [ -s $norm_fq ]; then
     
          $bin_dir/unicycler  \
               -s ${norm_fq} \
               -o ${pair_id} \
               --kmers ${params.kmers} \
          || touch ${pair_id}.gfa

     fi

     if [ -f ${pair_id}/assembly.gfa ]; then
        mv ${pair_id}/assembly.gfa ${pair_id}.gfa
     else
        touch ${pair_id}.gfa
     fi
     """
}


process circularize {

     publishDir path: "$work_dir/${params.run}/${params.analysis}/${plate_id}_GFA", pattern: "*.gfa", mode: 'copy'
     publishDir path: "$work_dir/${params.run}/${params.analysis}/${plate_id}_FASTA", pattern: '*.fasta', mode: 'copy'

     input:
     tuple val(pair_id), val(plate_id), path(gfa) from unicycler_gfa

     output:
     tuple val(pair_id), path("*_???.final.gfa"), path("*_???.final.fasta") into circle_out
     tuple val(pair_id),  path('*_???.info.csv') into circle_csv

     """
     if [ -s $gfa ]; then
          /root/wf-nf/bin/Graph.py 
     else
          touch ${pair_id}.final.gfa
          touch ${pair_id}.final.fasta
          touch ${pair_id}.info.csv
     fi

     """
}

align_in = fq_align.join(circle_out)

process bwa {

     publishDir path: "$work_dir/${params.run}/${params.analysis}/${plate_id}_BWA", pattern: "*.csv", mode: 'copy'
    	publishDir path: "$work_dir/${params.run}/${params.analysis}/${plate_id}_BWA", pattern: "*.bam", mode: 'copy'

     input:
     tuple val(pair_id), val(plate_id), path(fq), path(gfa), path(fa) from align_in

     output:
     path "*bam" into bam
     tuple val(pair_id), val(plate_id), path ("*.csv") into metrics

     """
     if [ -s $fa ]; then
          $bin_dir/bwa index $fa

          $bin_dir/bwa mem $fa $fq \
          | $bin_dir/samtools view -bh -F2048 - \
          | $bin_dir/samtools sort > ${pair_id}.bam

          $bin_dir/samtools depth -a ${pair_id}.bam > ${pair_id}.depth.csv
          $bin_dir/samtools view -c ${pair_id}.bam >${pair_id}.count.csv
          
     else
          touch ${pair_id}.depth.csv
          touch ${pair_id}.count.csv
          touch ${pair_id}.bam
     fi
     """
}

summary_ch = metrics.join(circle_csv)
                    .groupTuple(by: 1)
                    .map{ plates -> tuple(plates[1], plates[2].flatten() + plates[3]) }

process summarize {

     publishDir path: "$work_dir/${params.run}/${params.analysis}/", pattern: "*.csv", mode: 'copy'
    	publishDir path: "$work_dir/${params.run}/${params.analysis}/${plate_id}_FIGS", pattern: "*.png", mode: 'copy'

	input: 
	tuple val(plate_id), path(metrics) from summary_ch

        output:
        path("*") into summary_output
//        path("*png") into coverage_figures

	"""
	/root/wf-nf/bin/SummarizeAssembly.py $plate_id
	"""

}
