workflow rqcfilter2 {
    Array[File] input_files
    Boolean nersc = true

    String bbtools_container="bryce911/bbtools:38.86"
    call rqcfilter {
    	 input: reads_files=input_files, is_nersc=nersc, container=bbtools_container
    }
    output {
        Array[File] final_filtered = rqcfilter.outrqcfilter
    }
}

task rqcfilter{
    Array[File] reads_files

    String container
    Boolean is_nersc
    String run_prefix = if(is_nersc) then "shifter --image=" + container + " -- " else ""    
    String single = if (length(reads_files) == 1 ) then "1" else "0"

    String rqcfilter_input = "rqcfilter.input.fastq.gz"
    String rqcfilter_output = "rqcfilter.output.fastq.gz"    

     String filename_readlen="readlen.txt"
     String filename_outlog="stdout.log"
     String filename_errlog="stderr.log"

     String java="-Xmx20g"
     String dollar="$"
     runtime {
     	     docker: container
     } 

    command {
        if [ ${single} == 0 ]
	then
	    cat ${sep = " " reads_files } > ${rqcfilter_input}
	else
	    cp ${reads_files[0]} ./${rqcfilter_input}
	fi
	touch ${filename_readlen}
	${run_prefix} readlength.sh -Xmx1g in=${rqcfilter_input} out=${filename_readlen} overwrite ;
        ${run_prefix} rqcfilter2.sh -Xmx45g jni=t in=${rqcfilter_input} rqcfilterdata=/data/RQCFilterData/ \
        path=filter rna=f trimfragadapter=t qtrim=r trimq=0 maxns=3 maq=3 minlen=51 \
        mlf=0.33 phix=t removehuman=t removedog=t removecat=t removemouse=t khist=t \
        removemicrobes=t sketch kapa=t clumpify=t tmpdir= barcodefilter=f trimpolyg=5 usejni=f \
        out=${rqcfilter_output} ;
     }
     output {
     	    Array[File] outrqcfilter = glob("filter/${rqcfilter_output}")
            File outreadlen = filename_readlen

     }

}
#            File stdout = filename_outlog
#            File stderr = filename_errlog
