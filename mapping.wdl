workflow mapping {
    Array[File] input_files
    File input_reference
    Boolean dofinalize_bams=true
    Boolean dotar_bams=false
    String bbtools_container="bryce911/bbtools:38.86"

    if (length(input_files) == 1 ){
       call mappingtask as single_run {
           input: reads=input_files[0], reference=input_reference, container=bbtools_container
       }
    }
    if (length(input_files) > 1 ){
       scatter (input_file in input_files) {
           call mappingtask as multi_run {
               input: reads=input_file, reference=input_reference, container=bbtools_container
           }
       }
    }
    if(dofinalize_bams){
    	call finalize_bams{
        	input: insing=single_run.outbamfile, inmult=multi_run.outbamfile, container=bbtools_container
    	}
    }
    if (dotar_bams){
    	call tar_bams{
    		input: insing=single_run.outbamfile, inmult=multi_run.outbamfile, container=bbtools_container
    	}
    }
}

task tar_bams {
    File? insing
    Array[File]? inmult
    String container

    String single = if(defined(insing)) then "1" else "0"
    String filename_tarbam="bamfiles.tar"
    String tarlist="bamfiles.txt"
    String dollar="$"

runtime {
docker: container
backend: "r5-120D-ceq"
memory: "120 GiB"
cpu:  16
}

    command{
        if [ ${single} == "1" ]
        then
                echo ${insing} > ${tarlist}
        else
                echo ${sep=" " inmult} | tr " " "\n" > ${tarlist}
        fi
	tar -cvf ${filename_tarbam} -T ${tarlist}
    }
    output{
        File outtarbam = filename_tarbam
    }
}   

task finalize_bams {
    File? insing
    Array[File]? inmult
    String container

    String single = if(defined(insing)) then "1" else "0"
    String java="-Xmx20g"
    String filename_outsam="pairedMapped.sam.gz"
    String filename_sorted="pairedMapped_sorted.bam"
    String filename_sorted_idx="pairedMapped_sorted.bam.bai"
    String filename_cov="pairedMapped_sorted.bam.cov"
    String filename_flagstat="pairedMapped_sorted.bam.flagstat"    
    String dollar="$"
    #runtime { docker: container}
runtime {
docker: container
backend: "r5-120D-ceq"
memory: "120 GiB"
cpu:  16
}

    command{
        if [ ${single} == "1" ]
        then
                cp ${insing} ${filename_sorted}      
        else
                samtools merge ${filename_sorted} ${sep=" " inmult}
        fi
        samtools index ${filename_sorted}
        reformat.sh threads=${dollar}(nproc) ${java} in=${filename_sorted} out=${filename_outsam} overwrite=true
        pileup.sh in=${filename_sorted} out=${filename_cov}
        samtools flagstat ${filename_sorted} 1>| ${filename_flagstat} 2>| ${filename_flagstat}.e          
    }
    output{
        File outsam = filename_outsam
        File outbam = filename_sorted
        File outbamidx = filename_sorted_idx
        File outcov = filename_cov
        File outflagstat = filename_flagstat    
    }
}

task mappingtask {
    File reads
    File reference
    String container

    String filename_unsorted="pairedMapped.bam"
    String filename_sorted="pairedMapped_sorted.bam"
    String dollar="$"

runtime {
docker: container
backend: "r5-120D-ceq"
memory: "120 GiB"
cpu:  16
}

    command{
        bbmap.sh threads=${dollar}(nproc)  nodisk=true \
        interleaved=true ambiguous=random rgid=filename \
        in=${reads} ref=${reference} out=${filename_unsorted}
        samtools sort -m100M -@ ${dollar}(nproc) ${filename_unsorted} -o ${filename_sorted};
  }
  output{
      File outbamfile = filename_sorted
   }
}

