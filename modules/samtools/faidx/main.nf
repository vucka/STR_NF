nextflow.enable.dsl=2

process FAIDX {
    tag "${sample_name}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path(fasta)
    
    output:
        path("*.fai"), emit: fai
    
    script:
        """
        echo "Generating .fa.fai file ....."

        samtools faidx ${fasta} > genome.fa.fai

        echo "INDEX file created!"
        """
}

workflow FAIDX_WF {
    take:
        ch_fasta
    main:
        FAIDX(ch_fasta)
    emit:
        fai = FAIDX.out.fai
}

workflow {
    ch_fasta = Channel.fromPath(params.fasta)
    ch_fasta.view()
    ch_fai = Channel.fromPath(params.output_dir + "/*.fai")
    ch_fai.view()

    FAIDX_WF(ch_fasta)
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}
