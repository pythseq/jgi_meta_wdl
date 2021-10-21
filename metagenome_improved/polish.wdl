import "index.wdl" as index
import "align_pbmm2.wdl" as align
import "sort.wdl" as sort
import "racon_polish.wdl" as racon

workflow polish {
    String smrtlink_container="bryce911/smrtlink:9.0.0.92188"
    String racon_container="quay.io/biocontainers/racon:1.4.13--he513fc3_0"
    File input_fastq
    File ref

    call index.index as index_round1 {
        input: round= "rd1" + ".mmi", ref=ref
    }
    call align.pbmm2 {
        input: index=index_round1.outidx, container=smrtlink_container, fastq=input_fastq, bam="rd1.bam"
    }

    call sort.sort {
         input: bam=pbmm2.outbam, container=smrtlink_container, filename_sam="round1.filtered.sam"
    }
    call racon.racon {
         input: ref=ref, container=racon_container, input_sam=sort.outsam, input_fastq=input_fastq, filename_polished="polished_assembly.rd1.fasta"
    }
    call index.index as index_round2 {
        input: round= "rd2" + ".mmi", ref=racon.outfasta
    }
    call align.pbmm2  as align_round2 {
        input: index=index_round2.outidx, container=smrtlink_container, fastq=input_fastq, bam="rd2.bam"
    }
        call sort.sort as sort_round2 {
         input: bam=align_round2.outbam, container=smrtlink_container, filename_sam="round2.filtered.sam"
    }
    call racon.racon as racon_round2 {
         input: ref=racon.outfasta, container=racon_container, input_sam=sort_round2.outsam, input_fastq=input_fastq, filename_polished="polished_assembly.fasta"
    }
    output {
        File final_polished=racon_round2.outfasta
    }
}
