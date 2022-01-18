workflow minimap2 {
  meta {
    version: '0.0.1'
    author: 'Robert Riley'
  }

  #File input_fastq
  Array[File] input_fastq
  File assembly_fasta
  String minimap2_container
  String minimap2_parameters
  String samtools_container
  String bbtools_container

  call run_minimap2 {
    input: input_fastq=input_fastq,
           assembly_fasta=assembly_fasta,
           minimap2_container=minimap2_container,
           minimap2_parameters=minimap2_parameters
  }

  call run_samtools {
    input: samtools_container=samtools_container,
           assembly_fasta=assembly_fasta,
           sam=run_minimap2.output_sam
  }

  call run_pileup {
    input: bbtools_container=bbtools_container,
           ref=assembly_fasta,
           sam=run_minimap2.output_sam
  }

  call run_samtools_stats {
    input: samtools_container=samtools_container,
           sam=run_minimap2.output_sam
  }
}



task run_minimap2 {
  Array[File] input_fastq
  File assembly_fasta
  String minimap2_container
  String minimap2_parameters
  String output_sam_filename = "contigs.sam"
  command {
    shifter -i=${minimap2_container} minimap2 \
    ${minimap2_parameters} \
    ${assembly_fasta} ${sep = " " input_fastq} > ${output_sam_filename}
  }
    #-t 16 -a -x map-pb \
    # Was formerly:
    #${assembly_fasta} ${input_fastq} > ${output_sam_filename}

  runtime {

  }

  output {
    File output_sam = output_sam_filename
  }
}


task run_samtools {
  File assembly_fasta
  File sam
  String samtools_container
  String output_bam_filename = "contigs.sorted.bam"
  command {
    shifter -i=${samtools_container} samtools view \
    -bT \
    ${assembly_fasta} ${sam} | \
    shifter -i=${samtools_container} samtools sort \
    -@ 32 \
    -o ${output_bam_filename}

  }

  runtime {

  }

  output {
    File output_bam = output_bam_filename
  }
}


task run_pileup {
  File ref
  File sam
  String bbtools_container
  String output_hist_filename = "contigs.sorted.bam.pileup.hist" 
  String output_basecov_filename = "contigs.sorted.bam.pileup.basecov"
  String output_bincov_filename = "contigs.sorted.bam.pileup.bincov"
  String output_pileup_filename = "contigs.sorted.bam.pileup.out"
  String output_stderr_filename = "stderr"
  command {
    shifter -i=${bbtools_container} pileup.sh \
    ref=${ref} \
    in=${sam} \
    hist=${output_hist_filename} \
    basecov=${output_basecov_filename} \
    bincov=${output_bincov_filename} \
    out=${output_pileup_filename}

  }

  runtime {

  }

  output {
    File output_hist = output_hist_filename
    File output_basecov = output_basecov_filename
    File output_bincov = output_bincov_filename
    File output_pileup = output_pileup_filename
    File output_stderr = output_stderr_filename
  }

}


task run_samtools_stats {
  String samtools_container
  File sam
  String output_stats_filename = "contigs.sam.stats"
  String summary_stats_filename = "summary.stats"
  command <<<
    shifter -i=${samtools_container} samtools stats ${sam} > ${output_stats_filename}

    # Pull some stats out of file for readme
    mapped=`grep "reads mapped:" ${output_stats_filename} | awk '{print $NF}'`
    total=`grep "raw total sequences:" ${output_stats_filename} | awk '{print $NF}'`

    echo -e "Input mapping read count:"'\t'$total >> ${summary_stats_filename}
    echo -e "Mapped read count:"'\t\t'$mapped "("`echo "a=100*$mapped/$total; scale=2; (a+0.005)/1" | bc -l`"%)" >> ${summary_stats_filename}
  >>>

  runtime {

  }

  output {
    File output_stats = output_stats_filename
    File summary_stats = summary_stats_filename
  }
}
