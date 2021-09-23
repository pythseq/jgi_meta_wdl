workflow flye {
  meta {
    version: '0.0.1'
    author: 'Robert Riley'
  }
  #File input_fastq
  Array[File] input_fastq
  String flye_container
  String flye_parameters
  
  call run_flye {
    input: input_fastq=input_fastq,
           flye_container=flye_container,
           flye_parameters=flye_parameters
  }
}


task run_flye {
  Array[File] input_fastq
  String flye_container
  String flye_parameters
  String single = if (length(input_fastq) == 1 ) then "1" else "0"
  String combined_fastq_filename = "reads.flye.fastq.gz"
  String assembly_fasta_filename = "flye/assembly.fasta"
  String flye_log_filename = "flye/flye.log"
  String assembly_info_filename = "flye/assembly_info.txt"

  command {
    # Combine seq units, if multiple, into one, or just rename to reads.input.fastq.gz
    if [ ${single} == 0 ]
    then
      cat ${sep = " " input_fastq } > ${combined_fastq_filename}
    else
      ln -s ${input_fastq[0]} ./${combined_fastq_filename}
    fi

    # Run flye
    shifter --image=${flye_container} flye ${flye_parameters} ${combined_fastq_filename}
  }

  runtime {

  }

  output {
    File combined_fastq = combined_fastq_filename
    File assembly_fasta = assembly_fasta_filename
    File flye_log = flye_log_filename
    File assembly_info = assembly_info_filename
  }
}

