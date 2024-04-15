nextflow.enable.dsl=2

process EXPANSIONHUNTER {
    tag "${sample_name}"
    publishDir "${params.output_dir}", mode: 'copy'
        
    input:
        tuple val(sample_name), path(bam)
        path(fasta)
        path(variant_catalog)
        tuple val(sample_name), path(bai)
    
    output:
        tuple val(sample_name), path("*.vcf"), emit: vcf
        tuple val(sample_name), path("*_exph.json"), emit: json
    
    script:
        """
        echo "Running ExpansionHunter ....."

        ExpansionHunter --reads ${bam} --reference ${fasta} --variant-catalog ${variant_catalog} --output-prefix ${sample_name}_exph
        
        echo "ExpansionHunter finished!"
        """
}

workflow EXPANSIONHUNTER_WF {
    take: 
        ch_sorted_bam
        ch_fasta
        ch_variant_catalog
        ch_bai
    main:
        EXPANSIONHUNTER(ch_sorted_bam, ch_fasta, ch_variant_catalog, ch_bai)
    emit: 
        vcf = EXPANSIONHUNTER.out.vcf
        json = EXPANSIONHUNTER.out.json
}

workflow {
    ch_sorted_bam = Channel.fromPath(params.output_dir + "/*_sorted.bam")
    ch_fasta = Channel.fromPath(params.fasta)
    ch_variant_catalog = Channel.fromPath(params.variant_catalog)

    EXPANSIONHUNTER_WF(ch_sorted_bam, ch_fasta, ch_variant_catalog, ch_bai)
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}