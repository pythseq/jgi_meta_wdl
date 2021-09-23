workflow samtools_sort {
  File input_bam
  String smrtlink_container #="bryce911/smrtlink:10.0.0.108728"

  call sort {
    input: bam=input_bam,
    container=smrtlink_container,
    filename_sam="round1.filtered.sam"
  }
}

task sort {
  File bam
  String container
  String round
  String filename_errlog = "samtools_sort_" + round + "_stderr.log"
  String filename_sam = round + ".filtered.sam"

  command {
    shifter --image=${container} samtools view -F 1796 -q 20  ${bam} 1> ${filename_sam} 2> ${filename_errlog}
  }

  runtime {
        
  }

  output {
    File outsam=filename_sam
    File errlog=filename_errlog
  }
}
