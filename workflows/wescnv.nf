// Check input path parameters to see if they exist
def checkPathParamList = [
    params.outdir,
    params.input,
    params.reference,
    params.target_bed
]

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

// Set input
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { INPUT_CHECK    } from '../subworkflows/input_check'
include { MOSDEPTH       } from '../modules/nf-core/mosdepth'
include { HTSNIMTOOLS_COUNT_READS     } from '../modules/local/hts_nim_tools/count-reads'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WESCNV {

    INPUT_CHECK( ch_input )



    // MOSDEPTH for avaerge coverage per capture region
    INPUT_CHECK.out.alignments.map { meta, aln, index ->
               return [meta, aln, index, file(params.target_bed)]
           }
        .set { ch_mosdepth_in }
    //ch_mosdepth_in.view()
    MOSDEPTH ( ch_mosdepth_in, [[:],[file(params.reference)]])
    HTSNIMTOOLS_COUNT_READS ( ch_mosdepth_in, [[:],[file(params.reference)]])

}
