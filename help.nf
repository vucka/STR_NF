nextflow.enable.dsl=2

def printHeader() {
  
  log.info """\
  CNV   P I P E L I N E
  ===================================
  fastq files    : ${ params.reads }
  \n
  """

}

def helpMessage() {

  yellow = "\033[0;37m"
  blue = "\033[0;35m"
  white = "\033[0m"
  red = "\033[0;31m"

  log.info """\
    cnv-pipeline 
    Usage:
        nextflow run main.nf [options]
    Script Options: see nextflow.config
        ${yellow}
        [required]
        --reads             FILE    Path to fastq files specified as a glob pattern
        OR
       
        [optional]
        
        --publish_dir       DIR     Path to run output directory
                                    DEFAULT: ${params.publish_dir}
     
        --email_on_fail     STR     Email to receive upon failure
                                    DEFAULT: ${params.email_on_fail}
                                    
        --help              BOOL    Display help message
        
    ${yellow}
    """.stripIndent()


}

workflow{
  printHeader()
  helpMessage()
}