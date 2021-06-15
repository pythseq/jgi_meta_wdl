workflow metatranscriptome_assy {
    Array[File] input_files

    String bbtools_container="bryce911/bbtools:38.86"
    String megahit_container_prod="vout/megahit:release-v1.2.9"

    call readstats_raw {
    	 input: reads_files=input_files, container=bbtools_container
    }
    call assy {
         input: reads_files=input_files, container=megahit_container_prod
    }
    call create_agp {
         input: contigs_in=assy.out, container=bbtools_container
    }
    output {
        File final_contigs = create_agp.outcontigs
        File final_scaffolds = create_agp.outscaffolds
        File final_log = assy.log
	    File final_readlen = readstats_raw.outreadlen
    }
}


task readstats_raw {
     Array[File] reads_files
     String container

     String single = if (length(reads_files) == 1 ) then "1" else "0"

     String reads_input="reads.input.fastq.gz"
     String outfile="readlen.txt"
     #runtime {backend: "Local"}
     command {
        if [ ${single} == 0 ]
	then
	    cat ${sep = " " reads_files } > ${reads_input}
	else
	    ln -s ${reads_files[0]} ./${reads_input}
	fi

        readlength.sh in=${reads_input} 1>| ${outfile}
     }
     output {
         File outreadlen = outfile
     }
}

task assy {
     Array[File] reads_files
     String container

     String single = if (length(reads_files) == 1 ) then "1" else "0"

     String reads_input="reads.input.fastq.gz"
     String outprefix="out.megahit"
     String filename_outfile="${outprefix}/final.contigs.fa"
     String filename_outfile_opts="${outprefix}/opts.txt"
     String filename_outfile_log="${outprefix}/log"
     String dollar="$"
     #runtime {backend: "Local"}

     command{
        if [ ${single} == 0 ]
	then
	    cat ${sep = " " reads_files } > ${reads_input}
	else
	    ln -s ${reads_files[0]} ./${reads_input}
	fi

	megahit -t ${dollar}(nproc) --k-list  23,43,63,83,103,123 -m 100000000000 -o ${outprefix} --12 ${reads_input}
	echo " -t ${dollar}(nproc) --k-list  23,43,63,83,103,123 -m 100000000000 -o ${outprefix} --12 ${reads_input}" >| ${filename_outfile_opts}	
     }
     output {
            File out = filename_outfile
	    File opt = filename_outfile_opts
	    File log = filename_outfile_log
     }
}

task create_agp {
    File contigs_in
    String container

    String java="-Xmx48g"
    String prefix="assembly"
    String filename_contigs="${prefix}.contigs.fasta"
    String filename_scaffolds="${prefix}.scaffolds.fasta"
    String filename_agp="${prefix}.agp"
    String filename_legend="${prefix}.scaffolds.legend"
    #runtime {backend: "Local"}

    command{
        fungalrelease.sh ${java} in=${contigs_in} out=${filename_scaffolds} outc=${filename_contigs} agp=${filename_agp} legend=${filename_legend} mincontig=200 minscaf=200 sortscaffolds=t sortcontigs=t overwrite=t
  }
    output{
	File outcontigs = filename_contigs
	File outscaffolds = filename_scaffolds
	File outagp = filename_agp
    	File outlegend = filename_legend
    }
}
