#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.cag_batchsize = 10000

// Default options
params.distance_threshold = 0.5
params.distance_metric = "cosine"
params.linkage_type = "average"
params.famli_batchsize = 10000000

include { makeInitialCAGs } from "./cag_util" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round1 } from "./cag_util" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round2 } from "./cag_util" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round3 } from "./cag_util" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round4 } from "./cag_util" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round5 } from "./cag_util" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round6 } from "./cag_util" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round7 } from "./cag_util" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round8 } from "./cag_util" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round9 } from "./cag_util" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round10 } from "./cag_util" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_final } from "./cag_util" params(
    distance_threshold: params.distance_threshold,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)



workflow CAG_contig_oriented_wf {
    //
    //  CAG Stuff starts here
    //    This is for the 'contig oriented'
    //
    take:
        gene_abundances_zarr_tar
        gene_lists

    main: 
    // Group shards of genes into Co-Abundant Gene Groups (CAGs)
    makeInitialCAGs(
        gene_abundances_zarr_tar,
        gene_lists
    )

    // Perform multiple rounds of combining shards to make ever-larger CAGs
    refineCAGs_round1(
        makeInitialCAGs.out[0].toSortedList().flatten().collate(2),
        makeInitialCAGs.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round2(
        refineCAGs_round1.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round1.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round3(
        refineCAGs_round2.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round2.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round4(
        refineCAGs_round3.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round3.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round5(
        refineCAGs_round4.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round4.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round6(
        refineCAGs_round5.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round5.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round7(
        refineCAGs_round6.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round6.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round8(
        refineCAGs_round7.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round7.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round9(
        refineCAGs_round8.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round8.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round10(
        refineCAGs_round9.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round9.out[1].toSortedList().flatten().collate(2),
    )

    // Combine the shards and make a new set of CAGs
    refineCAGs_final(
        refineCAGs_round10.out[0].toSortedList(),
        refineCAGs_round10.out[1].toSortedList(),
    )

    emit:
        cag_csv = refineCAGs_final.out[0]
        cag_abund_feather = refineCAGs_final.out[1]
}