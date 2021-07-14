workflow run_singlem {

 String bbtools_container="bryce911/bbtools:38.90"
 String singlem_container="wwood/singlem:0.13.2-dev11.a6cc1b4"
 File in_fastq
 String pf = basename(in_fastq, ".fastq.gz")
 Int threads
 
    call split_reads {
     input: reads=in_fastq, container=bbtools_container
    }
    call singlem {
     input: n_threads=threads, container=singlem_container, r1=split_reads.out1, r2=split_reads.out2, in_name=pf
    }
    call singlem_remove_off_target {
     input: input_csv=singlem.out_csv, container=singlem_container, in_name=pf
    }
    call singlem_cluster {
     input: input_target_csv=singlem_remove_off_target.out_targeted_csv, container=singlem_container,in_name=pf
    } 
    call singlem_krona {
     input: input_clustered_csv=singlem_cluster.out_clustered_csv, container=singlem_container,in_name=pf
    }
}

task split_reads {
   String container
   File reads
   String filename_outfile1="r1.fastq.gz"
   String filename_outfile2="r2.fastq.gz"
   command {
  reformat.sh in=${reads} out=${filename_outfile1} out2=${filename_outfile2}
   }

    runtime {
        docker : container
        }

   output {
            File out1=filename_outfile1
            File out2=filename_outfile2


   }
}

task singlem {
  File r1
  File r2
  String in_name
  Int n_threads
  String container
  String otu_table = "${in_name}.out.csv"
       
   command {
    /singlem/bin/singlem pipe --working-directory-tmpdir --singlem-packages `ls -d /pkgs/*spkg` --assignment-method diamond --diamond-prefilter --diamond-prefilter-performance-parameters '--block-size 0.5 --target-indexed -c1' --diamond-prefilter-db /pkgs/53_db2.0-attempt4.0.60.faa.dmnd --diamond-taxonomy-assignment-performance-parameters '--target-indexed -c1' --threads ${n_threads} --forward ${r1} --reverse ${r2} --otu-table ${otu_table}
   }
    runtime {
        docker : container
        }

   output {
      File out_csv=otu_table

   }
}

task singlem_remove_off_target {

  File input_csv
  String container 
  String in_name
  String otu_table_targeted = "${in_name}.out.targeted.csv"
   command {
     /singlem/bin/singlem summarise --singlem-packages `ls -d /pkgs/*spkg` --input-otu-table  ${input_csv} --output-otu-table ${otu_table_targeted} --exclude-off-target-hits 
   }
    runtime {
        docker : container
        }

   output {
      File out_targeted_csv=otu_table_targeted

   }
}

task singlem_cluster {
  File input_target_csv
  String container
  String in_name
  String otu_table_clustered = "${in_name}.out.targeted.clustered.csv"
   command {
     /singlem/bin/singlem summarise --cluster --cluster-id 0.89 --input-otu-table  ${input_target_csv} --output-otu-table ${otu_table_clustered} 
   }
    runtime {
        docker : container
        }

   output {
      File out_clustered_csv=otu_table_clustered

   }
}

task singlem_krona {
 File input_clustered_csv
 String container
 String in_name
 String plot_name = "${in_name}.out.targeted.clustered.krona.html"
    command {
     /singlem/bin/singlem summarise --input_otu_tables ${input_clustered_csv} --krona ${plot_name}
    }
    runtime {
        docker : container
        }

   output {
      File out_clustered_krona=plot_name

   }
}
