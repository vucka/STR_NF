# TEMPLATE NEXTFLOW PIPELINE

This is a template repo that will be used for BioX pipelines

## Structure:
    │── configs/ 
    │     ├── environment.config
    │     ├── genomes.config
    │     ├── modules.config
    │     ├── test.config
    │── containers/
    │     ├── Tool_name/
    │         ├── Dockerfile
    │── modules/
    │     ├── module_name/
    │         ├── main.nf
    |── src/
    │     ├── test.py
    │── tools/
    │     ├── Submodule_name/
    │         ├── containters/
    │         ├── main.nf
    │
    │── .gitignore
    │── README.md
    │── help.nf
    │── main.nf
    │── nextflow.config

- <ins>Configs dir:<ins>

    - contains all config files to be used in the pipeline

- <ins>Containers dir:<ins>
    
    - contains containers of all tools and Dockerimages to be used in the pipeline

- <ins>Modules dir:</ins>
    
    - contains all modules to be used in the pipeline

- <ins>Src dir:</ins>

    - contains external scripts (e.g. R-scripts, python scripts) that will be used in the main pipeline
    - another option is to do that by dockerizing scripts - the scripts then should be placed in <ins>Containers dir<ins>

- <ins>Tools dir:</ins>
    
    - contains all Submodules that will be used in the main pipeline
    - Adding submodules into tools directory:
        In order to add tools that are placed in different repo, go into tools directory and type:

    ```
        git submodule add [link of a tool repo]
    ```

    - When cloning a git repo locally in order to clone submodules placed in tools directroy run this commands:

    ```
        git submodule update --init --recursive
    ```
