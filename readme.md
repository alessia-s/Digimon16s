# DIGIMON16s    
![Senza nome](https://user-images.githubusercontent.com/84088778/118647347-89103f00-b7e1-11eb-87cf-248cc7abf5db.png)

*This code was developed by Anna Corrà, Davide Santinello and Alessia Strapazzon as part of a project carried out during the course of Microbial Metagenomics 
(Molecular Biology master degree) at the University of Padova under the supervision of Prof. Stefano Campanaro and Dr. Arianna Basile.*

## IDentIfication of 16s Genes In MicrObial geNomes

Taxonomic assignment is one of the fundamental steps for understanding the compostion of a microbial community.
The 16S rRNA gene is highly conserved among species of Archea and Bacteria but this genes also contain nine hypervariable regions that can provide a species-specific signature sequence useful for their identification.
Since the number of microbial isolates with a sequenced genome is rapidly increasing in public databases, there is a strong need of methods to efficiently identify and extract the 16S rRNA gene sequence from the genomic sequence. This can allow a rapid taxonomic assignment based on computational methods.
We developed a script which is able to perform a fast identification and extraction of the 16S rRNA gene within a bunch of reference genomes.

This script works in two steps:
1. search for the 16S rRNA gene in the genomes by performing an automatic nhmmer search with HMMER software;
2. use the nhmmer hits to extract the identified 16S sequence from the genomes.

## Requirements
This program has been designed to work with a specific folder structure:
- .fasta/.fna/.fa/.faa/.ffn files in a **subfolder** to the working directory. The user has to provide the name of the folder as an argument when launching the program from the command line.
- The .hmm file for nhmmer has to be in the **same** directory where the software is located.
```
Digimon16s-folder
	|------ <Your genomes folder>
	|	|------ <your genome file>
	|	|------ <your genome file>
	|	|------ ...
	|
	|------ bac.ssu.rnammer.hmm
	|------ arc.ssu.rnammer.hmm
	|------ <your HMM file> (optional)
	|------ digimon.py
```

### Files
The automated nhmmer step makes use of two Hidden Markov Model files, one for Bacteria and one for Archaea.
This program comes with two hmm models for the 16S gene, one for Archaeal genomes (arc.ssu.rnammer.hmm) and one for Bacterial genomes (bac.ssu.rnammer.hmm) obtained from [RNAmmer](http://www.cbs.dtu.dk/services/RNAmmer/). It is possible to work with a user-provided hmm model, but the program as not been tested for this.

### Programs
- [HMMER (v.3.3.3)](http://hmmer.org/)
- Python 3.8+

## Command line instructions
### Launching the program
The program is provided as a .py file. It is needed to provide appropriate permissions to .py file properties to allow the file to run as a program.
To launch the script, open a terminal session in the directory where the program is located and enter the command `./digimon.py`.
### Arguments
When launching the script, arguments are used to specify both the program parameters, and to be passed to the automated nhmmer command.
##### MANDATORY arguments
- genomes folder: you must point to the folder where the microbial genomes are located. These genomes will be searched by nhmmer, while DIGIMON-16S will retrieve the 16S gene sequence from the genome. Genomes folder must be a subfolder to the working directory on the program.\
`-f <folder_name>` or `--folder <folder_name>`.
- Query hmm: user can also define which .hmm you want nhmmer to use as query when looking for the 16S genes. The search can be performed using the hmm files provided by the software or another one selected by the user.\
	`-q b` or `--query b`: nhmmer will use “bac.ssu.rnammer.hmm” file (specific for bacterial genomes). In place of `b`, you can also use one of the following notations: ba, bac, bacteria.\
	`-q a` or `--query a`: nhmmer will use “arc.ssu.rnammer.hmm” file (specific for archaeal genomes). In place of `a`, you can also use one of the following: a, ar, archaea.\
	`-q <hmm_filename>` or `--query <hmm_filename>`: nhmmer will use your specified hmm file.
    
##### OPTIONAL arguments
- minimum 16S sequence length: provide to the script to export 16s gene sequences only if longer than this number of bases. This command only affects retrieval of the 16s gene sequence by our program and the ensuing multifasta file and info file. It does not affect nhmmer hits and nhmmer output and log files. Defaults to **0**\
`-l <number>` or `--length <number>`
- E-value: allows you to run the nhmmer commands with this specified E-value. Must be integer or float with “.” as decimal separator. Defaults to **0.0001**\
`-e <e-value>` or `--evalue <e-value>`

### Example
`$ ./digimon.py -f MAG_bact -q bac -e 0.01 -l 600` \
this command will operate on fasta files contained in the subfolder "MAG_bact", with "bac.ssu.rnammer.hmm" as query, "0.01" as E-value for nhmmer command, minimum length of 600 nucleotides for an identified 16s sequence to be considered as valid.
	
## Testing setup
Query sequences as hidden markov models(HMM) were obtained from [RNAmmer](http://www.cbs.dtu.dk/services/RNAmmer/).

As a testing dataset, we used 1628 metagenome-assembled genomes originating from multiple anaereobic digesters, imported as FASTA files, taken from _Campanaro S, et al. New insights from the biogas microbiome by comprehensive genome-resolved metagenomics of nearly 1600 species originating from multiple anaerobic digesters. Biotechnol Biofuels. 2020 Feb 24;13:25_

## Output
The software will generate 4 files:
- `My_16S_genes_<folder>_<query>.fasta`: multi-fasta file containing all the valid 16S sequences obtained. Each sequence is identified by the code of the scaffold in which it has been found (with an incremental `_<integer>` suffix in case of multiple alignments on the same scaffold) and the name of the genome;
- `My_16S_genes_info_<folder>_<query>.txt`: tabular file which stores the number of 16S genes found per genomes, the length of the sequences, the start and end position of the alignment; 
- `nhmmer_log_<folder>_<query>.txt`: file containing information printed out by the nhmmer command;
- `nhmmer_output_<folder>_<query>.txt`: file containing information on the alignments performed by nhmmer, from which our program will save the alignment positions. 
