# DIGIMON16s    
![Senza nome](https://user-images.githubusercontent.com/84088778/118647347-89103f00-b7e1-11eb-87cf-248cc7abf5db.png)

*This code was developed as part of a project carried out during the course of Microbial Metagenomics 
(Molecular Biology master degree) at the University of Padova. 
The project was supervised by Prof. Stefano Campanaro and Dr. Arianna Basile.*

## IDentIfication of 16s Genes In MicrObial geNomes

The taxonomic assignment is one of the fundamental steps for understanding the microbial community.
The 16s genes are highly conserved across species of Archea and Bacteria but this genes also contain nine hypervariable regions that can provide species-specific signature sequence useful for their identification.
Since much more species are isolated the amount of data is incresing;therefore there is the necessity to find a way to extract the 16s sequence and assing a taxonomy with a computational method.
We developed a software that is able to identify the 16s genes inside a reference genome and extract it so afterwords will be possible to assing a taxonomy.

This program works in two steps:
1. search for the 16s gene in the genomes by automating an nhmmer search from the HMMER software;
2. use the nhmmer hits to extract the identified 16s sequence from the genomes.

## Requirements
This program has been designed to work with a specific folder structure (see example directory structure.zip): 
- .fasta/.fna/.fa/.faa files in a **subfolder** to the current Working Directory. You are to specify the name of the folder as an argument when launching the program from command line.
- The .hmm file for nhmmer has to be in the **same** directory as the program.

### Files
The automated nhmmer step makes use of an Hidden Markov Model file.
This program comes with two hmm models for the 16s gene, one for Archeal genomes (arc.ssu.rnammer.hmm) and one for Bacterial genomes (bac.ssu.rnammer.hmm). It is possible to work with a user-provided hmm model, but the program as not been tested for this.

### Programs
- [HMMER (v.3.3.3)](http://hmmer.org/)
- Python 3.8+

## Command line instructions
### Launching the program
The program is provided as a .py file. It is necessary to allow the file to run as a program by setting the appropriate permission in file properties.
To launch the program, open a terminal session in the directory where the program is located and enter the command `./digimon.py`.
### Arguments
When launching the program, arguments are used to specify both our program’s parameters, and to be passed to the automated nhmmer command.
##### MANDATORY arguments
- genomes folder: you must point to the folder where the microbial genomes are located. These genomes will be searched by nhmmer and our program to retrieve their 16s genes. Genomes folder must be a subfolder to the program’s working directory.\
`-f <folder_name>` or `--folder <folder_name>`.
- Query hmm: you can indicate which .hmm you want nhmmer to use as query when looking for the 16s genes. You can use one of the provided hmm files or one of your pleasing.\
	`-q b` or `--query b`: nhmmer will use “bac.ssu.rnammer.hmm” file (indicated for bacterial genomes). In place of `b`, you can also use one of the following notations: ba, bac, bacteria.\
	`-q a` or `--query a`: nhmmer will use “arc.ssu.rnammer.hmm” file (indicated for archaeal genomes). In place of `a`, you can also use one of the following: a, ar, archaea.\
	`-q <hmm_filename>` or `--query <hmm_filename>`: nhmmer will use your specified hmm file.
    
##### OPTIONAL arguments
- minimum 16s sequence length: tell the program to export 16s gene sequences only if longer than this number of bases. This command only affects retrieval of the 16s gene sequence by our program and the ensuing multifasta file and info file. It does not affect nhmmer hits and nhmmer output and log files. Defaults to **0**\
`-l <number>` or `--length <number>`
- E-value: allows you to run the nhmmer commands with this specified E-value. Must be integer or float with “.” as decimal separator. Defaults to **0.0001**\
`-e <e-value>` or `--evalue <e-value>`
	
## Testing setup
Query sequences as hidden markov models(HMM) were obtained from [RNAmmer](http://www.cbs.dtu.dk/services/RNAmmer/).

As a testing dataset, we used 1628 metagenome-assembled genomes originating from multiple anaereobic digesters, imported as FASTA files, taken from _Campanaro S, et al. New insights from the biogas microbiome by comprehensive genome-resolved metagenomics of nearly 1600 species originating from multiple anaerobic digesters. Biotechnol Biofuels. 2020 Feb 24;13:25_

## Output
The software will generate 4 files:
- `My_16S_genes_<folder>_<query>.Fasta`: multi-fasta file containing all the valid 16S sequences obtained. Each sequence is identified by the code of the scaffold in which it has been found (with an incremental `_<integer>` suffix in case of multiple alignments on the same scaffold) and the name of the genome;
- `My_16S_genes_info_<folder>_ <query>.txt`: tabular file which stores the number of 16S genes found per genomes, the length of the sequences, the start and end position of the alignment; 
- `nhmmer_log_<folder>_<query>.txt`: file containing information printed out by the nhmmer command;
- `nhmmer_output_<folder>_<query>.txt`: file containing information on the alignments performed by nhmmer, from which our program will save the alignment positions. 
