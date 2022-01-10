################################################################################
# Run mapping to an assembly with minimap2
# Stand-alone workflow
# Create merged fastq file for the mapping.
################################################################################


# For mapping
import "minimap2.wdl" as minimap2

workflow map {

  meta {
    version: '0.0.1'
    author: 'Robert Riley'
  }

  Array[File] input_fastq
  File assembly_fasta
  String minimap2_container
  String minimap2_parameters
  String samtools_container
  String bbtools_container

  # Mapping
  call minimap2.run_minimap2 as mm2 {
    input: input_fastq=input_fastq,
           assembly_fasta=assembly_fasta,
           minimap2_container=minimap2_container,
           minimap2_parameters=minimap2_parameters
  }

  call minimap2.run_samtools as mm_sort {
    input: samtools_container=samtools_container,
           assembly_fasta=assembly_fasta,
           sam=mm2.output_sam
  }

  call minimap2.run_pileup as pileup {
    input: bbtools_container=bbtools_container,
           ref=assembly_fasta,
           sam=mm2.output_sam
  }

  call minimap2.run_samtools_stats as stats {
    input: samtools_container=samtools_container,
           sam=mm2.output_sam
  }
}
