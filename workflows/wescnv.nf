// Check input path parameters to see if they exist
def checkPathParamList = [
    params.outdir,
    params.input,
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
include { PARTITIONGS    } from '../modules/local/partition'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ACNE {

    INPUT_CHECK( ch_input )



    // big GS files can be partitioned for efficiency
    // 1: partition into batches of size partition_n. (affects PBF and GC model creation steps)
    // 2: also partition how many ind files are passed to the CNV caller.

    // note also that splitting uses python polars - so more efficient using uncompressed input
    if (params.partition) {
        if (!params.partition_n) {
            log.error "must inlude n samples to split large GS project"
            exit 1
        }
        // PARTITIONGS returns uncompressed
        PARTITIONGS(INPUT_CHECK.out.gsfiles, params.partition_n)
            .transpose()
            .map {
                meta, partition ->
                def meta_clone = meta.clone()
                partition_suffix = partition.baseName.split('-').last()
                meta_clone.sampleID = meta_clone.id
                meta_clone.partition = partition_suffix
                meta_clone.id=meta_clone.id+'_'+partition_suffix
                [ meta_clone, partition ]
            }
            .set { ch_pre_split }

        // pfb and gcmodel files, unique to each subbatch instantiated by PARTITIONGS
        MAKEPFB( ch_pre_split )

        PENNCNV_GC( MAKEPFB.out.output, params.gc_model )

        PENNCNV_GC.out.output
            .join(MAKEPFB.out.output)
            .map{
                meta, gc_files, pfb_files ->
                [ meta, gc_files, pfb_files ]
            }
            .set { ch_gc_pfb }

        //
        SPLITGS( ch_pre_split )
            .transpose()
            .set { ch_post_split }

        // now merge ind split, pfb and gc files
        ch_post_split
            .combine( ch_gc_pfb, by: 0 )
            .set{ ch_pre_penn_call }

        PENNCNV_DETECT( ch_pre_penn_call, params.hmm)

        /*
        cleanup of raw calls post ind calls
        can concatenate rawcnv files to make one file per sub-batch/partition.
        NB: Still need batch speific pfb for merging adjacebt CNV calls
        */
        // cat calls
        PENNCNV_DETECT.out.raw_call
            .groupTuple(by: 0)
            .set { ch_pre_cat_calls }

        CONCATENATE_PENN_CALLS( ch_pre_cat_calls )
        // cat logs
        PENNCNV_DETECT.out.raw_log
            .groupTuple(by: 0)
            .set { ch_pre_cat_logs }

        CONCATENATE_PENN_LOGS( ch_pre_cat_logs )

        // merge split calls
        CONCATENATE_PENN_CALLS.out.output
            .combine(MAKEPFB.out.output, by: 0)
            .set{ ch_merge_cat_calls}

        PENNCNV_MERGE( ch_merge_cat_calls )

        PENNCNV_MERGE.out.output
            .combine(CONCATENATE_PENN_LOGS.out.output, by: 0)
            .set{ ch_pre_filter }

        PENNCNV_FILTER( ch_pre_filter )

        PENNCNV_FILTER.out.output
            .map{
                meta, cnv, qcpass, qcsum ->
                def meta_clone = meta.clone()
                meta_clone.id = meta_clone.sampleID
                meta_clone.remove("partition")
                meta_clone.remove("sampleID")
                [ meta_clone, cnv, qcpass, qcsum ]
                }
            .groupTuple(by: 0)
            .view()
            .set{ ch_post_filter }
        CONCATENATE_PARTITIONS( ch_post_filter )

    } else {
        SPLITGS( INPUT_CHECK.out.gsfiles )

    }
    //BATCH_CALL(params.runID, ch_input_sample)
}
