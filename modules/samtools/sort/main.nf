nextflow.enable.dsl=2

process SORT {
    tag "${sample_name}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        tuple val(sample_name), path(bam)
    
    output:
        tuple val(sample_name), path("${sample_name}_sorted.bam"), emit: sorted_bam
        tuple val(sample_name), path("${sample_name}_sorted.bam.bai"), emit: bai
    
    script:
        """
        echo "Converting SAM to BAM ....."

        samtools sort ${sample_name}_unsorted.bam -o ${sample_name}_sorted.bam
        samtools index ${sample_name}_sorted.bam 

        echo "BAM file created!"
        """
}

workflow SORT_WF {
    take:
        ch_bam
    main:
        SORT(ch_bam)
    emit:
        sorted_bam = SORT.out.sorted_bam
        bai = SORT.out.bai

}

workflow {
    ch_bam = Channel.fromPath(params.output_dir + "/*_unsorted.bam")

    SORT_WF(ch_bam)
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}
