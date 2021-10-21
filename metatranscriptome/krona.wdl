workflow krona {
    Array[File] input_files

    String krona_container="foster505/krona:v1.1.7"
    call run_krona {
    	 input: reads=input_files, container=krona_container
    }
}

task run_krona {
     Array[File] reads
     String container

     String single = if(length(reads)==1) then "1" else "0"
     String input_reads="input.fastq.gz"
     runtime { backend: "Local"}

    command{
	cat ${sep=" " reads} > ${input_reads}
	metagenome_krona.py --fastq ${input_reads}
    }
    output{
	Array[File] outrqcfilter = glob("*.html")
    }
}
