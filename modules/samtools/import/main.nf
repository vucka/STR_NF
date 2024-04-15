nextflow.enable.dsl=2

process IMPORT {
    tag "${sample_name}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        tuple val(sample_name), path(sam)
    
    output:
        tuple val(sample_name), path("${sample_name}_unsorted.bam"), emit: bam
    
    script:
        """
        echo "Converting SAM to BAM ....."

        samtools view -S -b ${sam} > ${sample_name}_unsorted.bam

        echo "BAM file created!"
        """
}

workflow IMPORT_WF {
    take:
        ch_sam
    main:
        IMPORT(ch_sam)
    emit:
        bam = IMPORT.out.bam
}

workflow {
    ch_sam = Channel.fromPath(params.output_dir + "/*.sam")

    IMPORT_WF(ch_sam)
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}
