nextflow.enable.dsl=2

include { FQ_TO_SAM_WF } from './modules/fq_to_sam/main.nf' params ( params )
include { FAIDX_WF } from './modules/samtools/faidx/main.nf' params ( params )
include { IMPORT_WF } from './modules/samtools/import/main.nf' params ( params )
include { SORT_WF } from './modules/samtools/sort/main.nf' params ( params )
include { EXPANSIONHUNTER_WF } from './modules/expansion_hunter/main.nf' params ( params )
include { REPORTING_WF } from './modules/reporting/main.nf' params ( params )

workflow STR_WF {
    take:
        ch_reads
        ch_fasta
        ch_variant_catalog
    main:
       
       
        FQ_TO_SAM_WF(ch_fasta,ch_reads)

        FAIDX_WF(ch_fasta)

        ch_sam = FQ_TO_SAM_WF.out.sam

        IMPORT_WF(ch_sam)

        ch_bam = IMPORT_WF.out.bam

        SORT_WF(ch_bam)

        ch_sorted_bam = SORT_WF.out.sorted_bam
        ch_bai = SORT_WF.out.bai

        EXPANSIONHUNTER_WF(ch_sorted_bam, ch_fasta, ch_variant_catalog, ch_bai)

        ch_exph_json = EXPANSIONHUNTER_WF.out.json

        REPORTING_WF(ch_exph_json)
}

workflow {
    // ch_reads = Channel.fromFilePairs(params.fastq_pairs + '/*_R{1,2}.fq.gz', size: -1)
    ch_reads = Channel.fromFilePairs(params.fastq_pairs + '/*_R{1,2}.fastq.gz', size: -1)
    ch_fasta = Channel.fromPath(params.fasta)
    ch_variant_catalog = Channel.fromPath(params.variant_catalog)
    ch_sam = Channel.fromPath(params.output_dir + "/*.sam")

    STR_WF(ch_reads, ch_fasta, ch_variant_catalog)
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}