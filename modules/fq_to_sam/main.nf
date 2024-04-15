nextflow.enable.dsl=2

process FASTA_INDEX {

  publishDir "${params.output_dir}", mode: 'copy'

  input:
 
    path(fasta)
  
  output:

   path "*.{amb,ann,bwt,pac,sa}", emit: bwa_index
  
  script:
    """

    bwa index ${fasta}


    """
}

process FQ_TO_SAM {
  tag "${sample_name}"
  publishDir "${params.output_dir}", mode: 'copy'

  input:
    tuple val(sample_name), path(reads)
    path(fasta)
    path(fasta_index)

  output:
    
    tuple val(sample_name), path ("${sample_name}.sam"), emit: sam
  
  script:
    """
    echo "Converting FASTQs to SAM ....."

    bwa mem -M -t 4 -R '@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA' ${fasta} ${reads} > ${sample_name}.sam
    
    
    echo "SAM file created!"
    """
}

workflow FQ_TO_SAM_WF{
  take: 
    ch_fasta
    ch_reads
    

  main:
    FASTA_INDEX(ch_fasta)
    ch_fasta_index = FASTA_INDEX.out.bwa_index
    FQ_TO_SAM(ch_reads, ch_fasta,ch_fasta_index)
  
  emit:
   
    sam = FQ_TO_SAM.out.sam
}

workflow {
  ch_reads = Channel.fromFilePairs(params.fastq_pairs + '/*_R{1,2}.fq.gz', size: -1)
  ch_fasta = Channel.fromPath(params.fasta)
  

  FQ_TO_SAM_WF(ch_fasta,ch_reads)
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}
