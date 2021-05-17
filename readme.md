*This code was developed as part of a project carried out during the course of Microbial Metagenomics 
(Molecular Biology master degree) at the University of Padova. 
The project was supervised by Prof. Stefano Campanaro and Dr. Arianna Basile.*

## Introduction
The purpose of this program is to check for the presence of a 16s gene in bacterial/archaeal genomes, ideally Metagenome-Assembled Genomes. It does so in two steps:
1. search for the 16s gene in the genomes by automating an nhmmer search
2. use the nhmmer hits to extract the identified 16s sequence from the genomes.

## Requirements
This program comes with two hmm models for the 16s gene, one for Archeal genomes (arc.ssu.rnammer.hmm) and one for Bacterial genomes (bac.ssu.rnammer.hmm). It is possible to work with a user-provided hmm model, but the program as not been tested for this.

This program has been designed to work with a specific folder structure (see example directory): 
- .fasta/.fna/.fa files in a subfolder to the current Working Directory. You are to specify the name of the folder as an argument when launching the program from command line.
- The .hmm file for hnmmer has to be in the same directory as the program.

### Programs
- HMMER 3.3 (http://hmmer.org/)
- Python 3.8+

## Command line instructions
### Launching the program
The program is provided as a .py file. It is necessary to allow the file to run as a program by setting the appropriate permission in file properties.
To launch the program, open a terminal session in the directory where the program is located and enter the command `./digimon.py`.
### Arguments
When launching the program, arguments are used to specify both our program’s parameters, and to be passed to the automated nhmmer command.
##### MANDATORY arguments
- genomes folder: you must point to the folder where the microbial genomes are located. These genomes will be searched by nhmmer and our program to retrieve their 16s genes. (genomes folder must be a subfolder to the program’s working directory).\
`-f <folder_name>` or `--folder <folder_name>`.
- Query hmm: you can indicate which .hmm you want nhmmer to use as query when looking for the 16s genes. You can use one of the provided hmm files or one of your pleasing.\
	`-q b` or `--query b`: nhmmer will use “bac.ssu.rnammer.hmm” file (indicated for bacterial genomes). In place of `b`, you can also use one of the following notations: ba, bac, bacteria.\
	`-q a` or `--query a`: nhmmer will use “arc.ssu.rnammer.hmm” file (indicated for archaeal genomes). In place of `a`, you can also use one of the following: a, ar, archaea.\
	`-q <hmm_filename>` or `--query <hmm_filename>`: nhmmer will use your specified hmm file.
    
##### OPTIONAL arguments
- minimum 16s sequence length: tell the program to export 16s gene sequences only if longer than this number of bases. This command only affects retrieval of the 16s gene sequence by our program and the ensuing multifasta file and info file. It does not affect nhmmer hits and nhmmer output and log files. Defaults to 0\
`-l <number>` or `--length <number>`
- E-value: allows you to run the nhmmer commands with this specified E-value. Must be integer or float with “.” as decimal separator. Defaults to 0.0001\
`-e <e-value>` or `--evalue <e-value>`
	
	
## Output
