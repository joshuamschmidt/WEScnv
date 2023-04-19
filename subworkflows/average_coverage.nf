/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAKEPFB        } from '../modules/local/make_pfb'
include { SPLITGS        } from '../modules/local/split_gs'
include { PENNCNV_GC     } from '../modules/local/penn_gc'


workflow AVERAGE_COVERAGE {
  take:
    gs_file

    main:
    MAKEPFB(gs_file)
    SPLITGS(gs_file)
    Channel
        .from( MAKEPFB.out )
        .first()
        .set { gc_in_ch }
    gc_in_ch.view()
    //PENNCNV_GC(gc_in_ch, Channel.fromPath(params.gc_model))
}
