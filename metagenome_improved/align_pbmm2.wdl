workflow align_pbmm2 {
  File input_index
  File input_fastq
  String container
  String round
  call pbmm2 {
    input: index = input_index,
    container = container,
    fastq = input_fastq,
    round = round  #bam="rd1.bam"
  }
}

task pbmm2 {
  String container
  File index
  File fastq
  String round
  String filename_outlog = "align_pbmm2_" + round + "_stdout.log"
  String filename_errlog = "align_pbmm2_" + round + "_stderr.log"
  String bam = round + ".bam"
  command {
    shifter --image=${container} pbmm2 align --preset="CCS" --sort ${index} ${fastq} ${bam} 1> ${filename_outlog} 2> ${filename_errlog}
  }
  runtime {
        
  }
  output {
	  File outbam = bam
    File outlog = filename_outlog
    File errlog = filename_errlog
  }
}
