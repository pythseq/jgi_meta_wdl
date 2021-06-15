import "metatranscriptome_assy.wdl" as metatranscriptome_assy
import "mapping.wdl" as mapping
import "krona.wdl" as krona

workflow metatranscriptome_assembly_and_alignment {
    Array[File] input_files
    Array[File] raw_files = []
    Boolean nersc = false

    call metatranscriptome_assy.metatranscriptome_assy as assy {
        input: input_files=input_files, nersc=nersc
    }
    call mapping.mapping {
       input: input_files=input_files, input_reference=assy.final_contigs, nersc=nersc
    }
    if (length(raw_files) > 0 ){    
       call krona.krona {
    	   input: input_files=raw_files, nersc=nersc
    	}
    }
}

