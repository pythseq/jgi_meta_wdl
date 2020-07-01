import "rqcfilter2.wdl" as rqcfilter2
import "metagenome_assy.wdl" as metagenome_assy
import "mapping.wdl" as mapping

workflow metagenome_filtering_assembly_and_alignment {
    Array[File] input_files
    Boolean nersc = false

	call rqcfilter2.rqcfilter2 as filter {
		input: input_files=input_files, nersc=nersc
	}
    call metagenome_assy.metagenome_assy as assy {
        input: input_files=filter.final_filtered, nersc=nersc
    }
    call mapping.mapping {
       input: input_files=filter.final_filtered, input_reference=assy.final_contigs, nersc=nersc
    }
}

