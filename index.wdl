workflow run_index {
  String round
  File ref
  String container
  call index {
    input: round = round,
           ref = ref,
           container = container
  }
}

task index {
  String container  #="bryce911/smrtlink:10.0.0.108728"
  File ref
  String round
  String outidx_filename = round + ".mmi"

  command {
    shifter --image=${container} pbmm2 index --preset="CCS" ${ref} ${outidx_filename} 1> index_${round}_stdout.log 2> index_${round}_stderr.log
  }

  runtime {
        
  }

  output {
    File outidx=outidx_filename
  }
}
