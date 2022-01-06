################################################################################
# Redo mapping to an assembly with minimap2
# Create merged fastq file for the mapping.
################################################################################


# For mapping
import "minimap2.wdl" as minimap2

workflow redo_mapping {

  meta {
    version: '0.0.1'
    author: 'Robert Riley'
  }

  Array[File] input_fastq
  File assembly_fasta
  String minimap2_container
  String samtools_container
  String bbtools_container
  String combined_fastq_filename = "reads.mm2.fastq.gz"



  # Concatenate input_fastq files to one file for minimap2 and Racon
  call combine_fastq {
    input: input_fastq = input_fastq,
           combined_fastq_filename = combined_fastq_filename,
           bbtools_container = bbtools_container
  }

  # Mapping
  call minimap2.run_minimap2 as map {
    input: input_fastq=combine_fastq.combined_fastq,
           assembly_fasta=assembly_fasta,
           minimap2_container=minimap2_container
  }

  call minimap2.run_samtools as mm_sort {
    input: samtools_container=samtools_container,
           assembly_fasta=assembly_fasta,
           sam=map.output_sam
  }

  call minimap2.run_pileup as pileup {
    input: bbtools_container=bbtools_container,
           ref=assembly_fasta,
           sam=map.output_sam
  }
  # TODO: produce readme?
}


# Tasks
task combine_fastq {
  Array[File] input_fastq
  String combined_fastq_filename
  String bbtools_container
  command {
    zcat ${sep = " " input_fastq } | shifter --image=${bbtools_container} reformat.sh in=stdin.fq out=${combined_fastq_filename} int=f qin=33
  }

  runtime {

  }

  output {
    File combined_fastq=combined_fastq_filename
  }
}

