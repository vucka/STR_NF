nextflow.enable.dsl=2

process EXTRACTING {
    tag "${sample_name}"
    publishDir "${params.output_dir}", mode: 'copy'
        
    input:
        tuple val(sample_name), path(json)
    
    output:
        tuple val(sample_name), path("*_reporting.json"), emit: reporting_json
    
    script:
        """
        echo "Extracting information from JSON file ....."

        python3 ~/../bin/exp_hunter_extraction.py ${sample_name}_exph.json ${sample_name}_reporting.json

        echo "Extraction complete ....."

        """
}

process REPORTING {
    tag "${sample_name}"
    publishDir "${params.output_dir}", mode: 'copy'
        
    input:
        tuple val(sample_name), path(json)
    
    output:
        tuple val(sample_name), path("*.pdf"), emit: pdf
    
    script:
        """
        echo "Generating STR report ....."

        python3 ~/../bin/generate_pdf_report.py ${sample_name}_reporting.json ${sample_name}_report.pdf
        
        echo "STR report generated!"
        """
}

workflow REPORTING_WF {
    take:
        ch_exph_json
    main:
        EXTRACTING(ch_exph_json)
        ch_json = EXTRACTING.out.reporting_json
        REPORTING(ch_json)
    emit:
        reporting_json = EXTRACTING.out.reporting_json
        pdf = REPORTING.out.pdf
}

workflow {
    ch_exph_json = Channel.fromPath(params.output_dir + "/*_exph.json")

    REPORTING_WF(ch_exph_json)
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}