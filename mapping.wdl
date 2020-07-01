workflow mapping {
    Array[File] input_files
    File input_reference
    Boolean nersc = false

    String bbtools_container="bryce911/bbtools:38.86"

    if (length(input_files) == 1 ){
       call mappingtask as single_run {
           input: reads=input_files[0], reference=input_reference, is_nersc=nersc, container=bbtools_container
       }
    }
    if (length(input_files) > 1 ){
       scatter (input_file in input_files) {
           call mappingtask as multi_run {
               input: reads=input_file, reference=input_reference, is_nersc=nersc, container=bbtools_container
           }
       }
    }
    call finalize_bams{
        input: insing=single_run.outbamfile, inmult=multi_run.outbamfile, is_nersc=nersc, container=bbtools_container
    }
    output {
        File final_outbam = finalize_bams.outbam
        File final_outsam = finalize_bams.outsam
        File final_outbamidx = finalize_bams.outbamidx
        File final_outcov = finalize_bams.outcov
        File final_outflagstat = finalize_bams.outflagstat
    }
}

task finalize_bams {
    File? insing
    Array[File]? inmult
    Boolean is_nersc
    String container

    String single = if(defined(insing)) then "1" else "0"
    String run_prefix = if(is_nersc) then "shifter --image=" + container + " -- " else ""    
    String java="-Xmx50g"
    String filename_outsam="pairedMapped.sam.gz"
    String filename_sorted="pairedMapped_sorted.bam"
    String filename_sorted_idx="pairedMapped_sorted.bam.bai"
    String filename_cov="pairedMapped_sorted.bam.cov"
    String filename_flagstat="pairedMapped_sorted.bam.flagstat"    
    String dollar="$"
#    runtime { backend: "Local"}

    command{
        if [ ${single} == "1" ]
        then
                ln -s ${insing} ${filename_sorted}      
        else
                ${run_prefix} samtools merge ${filename_sorted} ${sep=" " inmult}
        fi
        ${run_prefix} samtools index ${filename_sorted}
        ${run_prefix} reformat.sh threads=${dollar}(nproc) ${java} in=${filename_sorted} out=${filename_outsam} overwrite=true
        ${run_prefix} pileup.sh in=${filename_sorted} out=${filename_cov}
        ${run_prefix} samtools flagstat ${filename_sorted} 1>| ${filename_flagstat} 2>| ${filename_flagstat}.e          
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

    Boolean is_nersc = true
    String container

    String run_prefix = if(is_nersc) then "shifter --image=" + container + " -- " else ""
    String filename_unsorted="pairedMapped.bam"
    String filename_sorted="pairedMapped_sorted.bam"
    String dollar="$"
    
#    runtime { backend: "Local"}


    command{
        ${run_prefix} bbmap.sh threads=${dollar}(nproc)  nodisk=true \
        interleaved=true ambiguous=random rgid=filename \
        in=${reads} ref=${reference} out=${filename_unsorted}
        ${run_prefix} samtools sort -m200M -@ ${dollar}(nproc) ${filename_unsorted} -o ${filename_sorted}
  }
  output{
      File outbamfile = filename_sorted
   }
}

