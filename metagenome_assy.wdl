workflow metagenome_assy {
    Array[File] input_files

    String bbtools_container="bryce911/bbtools:38.90"
    String spades_container="bryce911/spades:3.15.2"

    call bbcms {
    	 input: reads_files=input_files, container=bbtools_container
    }
    call assy {
    	 input: infile1=bbcms.out1, infile2=bbcms.out2, container=spades_container    
    }
    call create_agp {
         input: scaffolds_in=assy.out, container=bbtools_container
    }
    output {
        File final_contigs = create_agp.outcontigs
        File final_scaffolds = create_agp.outscaffolds
        File final_spades_log = assy.outlog
	    File final_readlen = bbcms.outreadlen
	    File final_counts = bbcms.outcounts
    }
}

task bbcms{
    Array[File] reads_files

    String container
    String single = if (length(reads_files) == 1 ) then "1" else "0"

    String bbcms_input = "bbcms.input.fastq.gz"
    String filename_counts="counts.metadata.json"
  
    String filename_outfile1="input.corr.left.fastq.gz"
    String filename_outfile2="input.corr.right.fastq.gz"
    
     String filename_readlen="readlen.txt"
     String filename_outlog="stdout.log"
     String filename_errlog="stderr.log"
     String filename_kmerfile="unique31mer.txt"

     String java="-Xmx100g"
     String dollar="$"
     runtime { docker: container} 

     command {
        if [ ${single} == 0 ]
	then
	    cat ${sep = " " reads_files } > ${bbcms_input}
	else
	    ln -s ${reads_files[0]} ./${bbcms_input}
	fi
	touch ${filename_readlen}
	readlength.sh -Xmx1g in=${bbcms_input} out=${filename_readlen} overwrite 
        bbcms.sh ${java} metadatafile=${filename_counts} mincount=2 highcountfraction=0.6 \
	    in=${bbcms_input} out1=${filename_outfile1} out2=${filename_outfile2} \
	    1> ${filename_outlog} 2> ${filename_errlog}

     }
     output {
            File out1 = filename_outfile1
            File out2 = filename_outfile2
            File outreadlen = filename_readlen
            File stdout = filename_outlog
            File stderr = filename_errlog
            File outcounts = filename_counts

     }

}

task assy{
    File infile1
    File infile2    

    String container

    String outprefix="spades3"
    String filename_outfile="${outprefix}/scaffolds.fasta"
    String filename_spadeslog ="${outprefix}/spades.log"
    String dollar="$"
    runtime { docker: container}     

    command{
       spades.py -m 2000 --tmp-dir ${dollar}PWD -o ${outprefix} --only-assembler -k 33,55,77,99,127 --meta -t ${dollar}(nproc) -1 ${infile1} -2 ${infile2}
    }
    output {
           File out = filename_outfile
           File outlog = filename_spadeslog
    }
}

task create_agp {
    File scaffolds_in
    String container
    String prefix="assembly"
    
    String filename_contigs="${prefix}.contigs.fasta"
    String filename_scaffolds="${prefix}.scaffolds.fasta"
    String filename_agp="${prefix}.agp"
    String filename_legend="${prefix}.scaffolds.legend"
    String java="-Xmx40g"

    runtime { docker: container}         


    command{
        fungalrelease.sh ${java} in=${scaffolds_in} out=${filename_scaffolds} \
        outc=${filename_contigs} agp=${filename_agp} legend=${filename_legend} \
        mincontig=200 minscaf=200 sortscaffolds=t sortcontigs=t overwrite=t
    }
    output{
        File outcontigs = filename_contigs
        File outscaffolds = filename_scaffolds
        File outagp = filename_agp
        File outlegend = filename_legend
    }
}
