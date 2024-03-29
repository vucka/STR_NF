nextflow.enable.dsl=2

process sayHello {
  input: 
    val x
  output:
    stdout
  script:
    """
    echo '$x world!'
    """
}
workflow sayHello_WF {
    take:
        ch_x
    main:

    ch_x =  Channel.of('Bonjour', 'Ciao', 'Hello', 'Hola')
            sayHello(ch_x)
}
workflow {
     ch_x =  Channel.of('Bonjour', 'Ciao', 'Hello', 'Hola')
    
        sayHello_WF(ch_x)
        
}


// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}