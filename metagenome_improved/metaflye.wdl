################################################################################
# Assemble PacBio metagenome CCS reads, produce release-formatted fasta file, 
# mapping files, etc.
#
# WIP - getting rid of combined_fastq references
#     - can pbmm2, racon use multiple fastq?
################################################################################

# For assembling reads
import "flye.wdl" as flye

# For polishing with Racon
import "index.wdl" as index
import "align_pbmm2.wdl" as align
import "sort.wdl" as sort
import "racon_polish.wdl" as racon

# For mapping
import "minimap2.wdl" as minimap2

workflow metaflye {

  meta {
    version: '0.0.1'
    author: 'Robert Riley'
  }

  Array[File] input_fastq
  String flye_container
  String flye_parameters
  String smrtlink_container
  String racon_container
  String minimap2_container
  String minimap2_parameters
  String samtools_container
  String bbtools_container
  String combined_fastq_filename = "reads.fastq.gz"


  # Assemble the metagenome reads
  call flye.run_flye as assy {
    input: input_fastq = input_fastq,
           flye_container = flye_container,
           flye_parameters = flye_parameters
  }

  # Concatenate input_fastq files to one file for ppmm2 and Racon
  # pbmm2 needs either a combined file or fofn
  # Racon needs a combined file, so we have to make one
  # For flye and minimap2, combining fastqs is not necessary
  # so put fastqs individually on command line, more traceable if we have issues
  # with combining fastqs e.g. GAA-14582
  call combine_fastq {
    input: input_fastq = input_fastq,
           combined_fastq_filename = combined_fastq_filename,
           bbtools_container = bbtools_container
  }

# Polish the assembly with Racon
# Racon requires three inputs: contigs, reads, and read-to-contig mapping SAM file
  call index.index as index_round1 {
   input: round = "rd1",
    ref = assy.assembly_fasta,
    container = smrtlink_container
  }

  call align.pbmm2 {
    input: index = index_round1.outidx,
    container = smrtlink_container,
    fastq = combine_fastq.combined_fastq,
    round = "rd1" #bam="rd1.bam"
  }

  call sort.sort {
    input: bam=pbmm2.outbam,
    container=smrtlink_container,
    round = "rd1" #filename_sam="round1.filtered.sam"
  }

  call racon.racon {
    input: ref=assy.assembly_fasta,
    container=racon_container,
    input_sam=sort.outsam,
    input_fastq=combine_fastq.combined_fastq,
    round = "rd1",
    filename_polished="polished_assembly.rd1.fasta"
  }

  call index.index as index_round2 {
    #input: round= "rd2" + ".mmi", ref=racon.outfasta
    input: round= "rd2",
    ref=racon.outfasta,
    container = smrtlink_container
  }

  call align.pbmm2  as align_round2 {
    input: index = index_round2.outidx,
    container = smrtlink_container,
    fastq = combine_fastq.combined_fastq,
    round = "rd2" #bam = "rd2.bam"
  }

  call sort.sort as sort_round2 {
    input: bam=align_round2.outbam,
    container=smrtlink_container,
    round= "rd2" #filename_sam="round2.filtered.sam"
  }

  call racon.racon as racon_round2 {
    input: ref=racon.outfasta,
    container=racon_container,
    input_sam=sort_round2.outsam,
    input_fastq=combine_fastq.combined_fastq,
    round = "rd2",
    filename_polished="polished_assembly.fasta"
  }

  # Format polished assembly for release using fungalrelease.sh
  call format_assembly {
    input: bbtools_container=bbtools_container,
           polished_assembly=racon_round2.outfasta  #polish.polished_assembly
  }

  # Mapping
  call minimap2.run_minimap2 as map {
    input: input_fastq=input_fastq,
           assembly_fasta=format_assembly.asm_contigs,
           minimap2_container=minimap2_container,
           minimap2_parameters=minimap2_parameters
  }

  call minimap2.run_samtools as mm_sort {
    input: samtools_container=samtools_container,
           assembly_fasta=format_assembly.asm_contigs,
           sam=map.output_sam
  }

  call minimap2.run_pileup as pileup {
    input: bbtools_container=bbtools_container,
           ref=format_assembly.asm_contigs,
           sam=map.output_sam
  }

  call minimap2.run_samtools_stats as stats {
    input: samtools_container=samtools_container,
           sam=map.output_sam
  }

  call check_read_counts as check {
    input: input_fastq = input_fastq,
           combined_fastq=combine_fastq.combined_fastq
  }
  # TODO: produce readme?
}


# Tasks
task format_assembly {
  String bbtools_container
  File polished_assembly
  String asm_scaffolds_filename="assembly.scaffolds.fa"
  String asm_contigs_filename="assembly.contigs.fa"
  String asm_agp_filename="assembly.agp"
  String asm_legend_filename="assembly.legend"
  command {
    shifter --image=${bbtools_container} fungalrelease.sh -Xmx100g \
    in=${polished_assembly} out=${asm_scaffolds_filename} \
    outc=${asm_contigs_filename} agp=${asm_agp_filename} legend=${asm_legend_filename} \
    mincontig=500 minscaf=500 sortscaffolds=t sortcontigs=t overwrite=t
  }

  runtime {

  }

  output {
    File asm_scaffolds=asm_scaffolds_filename
    File asm_contigs=asm_contigs_filename
    File asm_agp=asm_agp_filename
    File asm_legend=asm_legend_filename
  }
}


task combine_fastq {
  Array[File] input_fastq
  String combined_fastq_filename
  String bbtools_container
  command {
    zcat ${sep = " " input_fastq} | shifter --image=${bbtools_container} reformat.sh in=stdin.fq out=${combined_fastq_filename} int=f qin=33
  }

  runtime {

  }

  output {
    File combined_fastq=combined_fastq_filename
  }
}


task check_read_counts {
  Array[File] input_fastq
  File combined_fastq
  String output_filename = "read_count_report.txt"
  command <<<
    echo "Input reads:" >> ${output_filename}
    total=0
    for fq in ${sep = " " input_fastq}
    do
      count=`echo $(zcat $fq | wc -l) /4 | bc`
      echo $fq $count >> ${output_filename}
      total=`expr $total + $count`
    done
    echo "Total = $total" >> ${output_filename}
    echo "Combined reads:" >> ${output_filename}
    fq=${combined_fastq}
    merged_total=`echo $(zcat $fq | wc -l) /4 | bc`
    echo "Merged total = $merged_total" >> ${output_filename}
    if [ $total == $merged_total ]
    then
      echo "Read counts match" >> ${output_filename}
    else
      echo "WARNING: read counts don't match" >> ${output_filename}
    fi
  >>>

  runtime {

  }

  output {
    File output_file=output_filename
  }
}
