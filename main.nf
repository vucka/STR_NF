nextflow.enable.dsl=2

include { printHeader; helpMessage } from './help' params ( params )
include { sayHello_WF } from './modules/test_module/main.nf' params ( params )


if ( params.help ) {
    helpMessage()
    exit 0
}





workflow {
    
    
     ch_x =  Channel.of('Bonjour', 'Ciao', 'Hello', 'Hola')
       
    sayHello_WF( 
                ch_x
                    
             )
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}