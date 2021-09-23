workflow racon_polish {
  File input_reference
  File input_sam
  File input_fastq
  String racon_container = "quay.io/biocontainers/racon:1.4.13--he513fc3_0"
  call racon {
    input: ref = input_reference,
           container = racon_container,
           input_sam = input_sam,
           input_fastq = input_fastq,
           round = "rd1",
           filename_polished = "polished_assembly.fasta"
  }
}

task racon {
  String container
  File ref
  File input_sam
  File input_fastq
  String round
  String filename_polished
  String filename_errlog = "racon_" + round + "_stderr.log"

  command {
    shifter --image=${container} racon -u -t 36  ${input_fastq} ${input_sam} ${ref} 1> ${filename_polished} 2> ${filename_errlog}
  }

  runtime {
        
  }

  output {
    File  outfasta = filename_polished
    File errlog = filename_errlog
  }
}
