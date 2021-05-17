# DIGIMON16s
## Introduction 

## IDentIfication of 16s Genes In MicrObial geNomes

The taxonomic assignment is one of the fundamental steps for understanding the microbial community.
The 16s genes are highly conserved across species of Archea and Bacteria but this genes also contain nine hypervariable regions that can provide species-specific signature sequence useful for their identification.
Since much more species are isolated the amount of data is incresing;therefore there is the necessity to find a way to extract the 16s sequence and assing a taxonomy with a computational method.
We developed a software that is able to identify the 16s genes inside a reference genome and extract it so afterwords will be possible to assing a taxonomy.

## Usage

For the software to work are necessary: 

HMMER (v.3.3.3) that was used for identify the 16s sequences; 

Python (v.3.8) to execute the software; 

Query sequences (since HMMER works also with querys and not only profiles) provided by an hidden markov models(HMM) one for the bacteria called:bac.ssu.rnammer.hmm and one for the archaea called:arc.ssu.rnammer.hmm and were obtained from RNAmmer "http://www.cbs.dtu.dk/services/RNAmmer/"

As a test, to verify if the program was working correctly, we used the 1600 metagenome-assembled genomes originating from multiple anaereobic digesters, imported as FASTA files.



### Output files 
The software will generate 4 file:
- My_16S_genes_( < folder >_<query> ).Fasta: it is a multi-fasta file containing all the valid 16S sequences obtained. Each sequence is identified by the code of the scaffold in which has been found and the name of the genome;
- My_16S_genes_info_(<folder>_<query>).txt: it is a tabular file in which are saved the number of 16S genes found per genomes, the length of the sequences, the start and end position of the alignment; 
- nhmmer_log_(<folder>_<query>).txt: file containing all the printed out information by the nhmmer software;
- nhmmer_output_(<folder>_<query>).txt: file containing information of the alignment performed by nhmmer, from which the software will save the alignment positions. 
