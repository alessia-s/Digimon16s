#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import re #regular expression
from math import floor # floor of the division
import os
import subprocess
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--folder", required=True) #permette di dichiarare la cartella in cui cercare gli fna, come argomento nel terminale
parser.add_argument("-q", "--query", required=True)     #permette di specificare se si vuole usare la query hmm per batteri o archea
parser.add_argument("-l", "--length", required=False)
parser.add_argument("-e", "--evalue", required=False)
args = parser.parse_args()


#*******************************************************************
if args.length is None:
    min_len = 0
else: 
    try:
        min_len = int(args.length)
    except ValueError:
        print("Minimum length must be an integer")
        quit()
        
if args.evalue is None:
    e_value = "0.0001"
else:
    try:
        e_value = float(args.evalue)
        e_value = str(e_value)
    except ValueError:
        print("e-value must be a number. \".\" as decimal separator")
        quit()
#***********************************************************************

if args.query.lower() in ["b", "ba", "bac", "bacteria"]:
    query = "bac.ssu.rnammer.hmm"
elif args.query.lower() in ["a", "ar", "arc", "archea"]:
    query = "arc.ssu.rnammer.hmm"
elif os.path.isfile(args.query) == True:
    query = args.query
else:
    print("Please, choose a valid query argument")
    quit()
    
print("\nlooking for fasta files in folder:", args.folder)
print("Query:", query)
print("searching for hits longer than:", min_len)
print("nhmmer E-value =", e_value)
print()

FNA_list=[]  #crea una lista vuota per i file fna su cui vogliamo lavorare
for files in sorted(os.listdir(args.folder)):   #uno ad uno i file presenti nella cartella dove cercare i fasta
    #aggiungiamo alla lista dei file solo quelli con estensione corretta, nel caso ci siano altri file nella cartella 
    if os.path.splitext(files)[1] in [".fna", ".fasta", ".fa", ".faa", ".ffn"]:
        FNA_list.append(files)
if FNA_list == []:
    print("No suitable fasta file in folder")
    quit()

processed_files = 0
processed_list = []
nhmmer_output_list = []
nhmmer_log=[]

for i in FNA_list:
    print("processing", args.folder + "/%s" % i)
    subprocess.run(["nhmmer", "-E", e_value, "-o", "log_temp.txt", "--dfamtblout=OUTPUTTEMP.txt", query , args.folder +"/%s" % i])
    processed_files += 1
    processed_list.append(i)
    nhmmer_output_list.append("@" + i +"\n" + open("OUTPUTTEMP.txt", "r").read())
    nhmmer_log.append("@" + i +"\n" + open("log_temp.txt", "r").read())

with open("nhmmer_output_%s_%s.txt" % (args.folder, query.replace(".", "")), "w") as f:
    print(*nhmmer_output_list, sep="\n", file=f)
with open("nhmmer_log_%s_%s.txt" % (args.folder, query.replace(".", "")), "w") as f:
    print(*nhmmer_log, sep="\n", file=f)
os.remove("OUTPUTTEMP.txt")
os.remove("log_temp.txt")
print ("\nProcessed", processed_files, "files")
print("Output:", "nhmmer_output_%s_%s.txt" % (args.folder, query.replace(".", "")))
print("Log:", "nhmmer_log_%s_%s.txt" % (args.folder, query.replace(".", "")))
print()


def get_the_gene(start, end, target, file_fasta, strand):
    
    #---------------------------    
    if strand == '-':
        i = end
        end = start
        start= i
    
    x = open(file_fasta,'r')
    file = x.read()
    x.close()

    file= file.replace('\n','')
    seq = file.split('>')
    
    target_s = '^'+ target + '\D'
    for s in seq:
        f = re.search(target_s, s)
        if f is not None:
            my_seq = s

    my_seq= my_seq.replace(target_name,'')

    gene_16S = my_seq[int(start) -1 :int(end)]
    
    if strand == '-':
        gene_16S = rev_complement(gene_16S)
        
    if len(gene_16S) < min_len:
        return 'not long enough'
    
    return gene_16S

def save_genes(gene, target, genome_name, my_16S):
    

    genome_name = genome_name.partition('.')[0]
    genome_name = genome_name.split('/')[1]
    head ='>'+ target + '.*' + genome_name
    
    count = 0
    for values in my_16S:
        f = re.search(head, values)
        if f is not None:
            count = count +1
            
    if count > 0:
        DNA ='>' + target + '_' + str(count) + '+' + genome_name + split_in_70(gene)
    else:
        DNA ='>' + target + '+' + genome_name + split_in_70(gene)
    
    my_16S = my_16S.append(DNA)
    
    return my_16S



def save_in_Fasta(my_16S):
        
    with open('My_16S_genes_%s_%s.fasta' %(args.folder, query.replace(".", "")), 'w') as f: 
        for i in my_16S:
            f.write(i)



def split_in_70(DNA):
    if len(DNA) > 0:
        if len(DNA) < 70:
            DNA = DNA + '\n'
            return '\n' + DNA
        else:
            cycle = floor(len(DNA)/70)
            i = 0
            dna = ''
            temp = DNA
            while i < cycle:
                peace = temp[0:70]
                dna = dna + peace + '\n'
                temp = temp[70:]
                i += 1
            return '\n' + dna + temp + '\n'
    return '\n'



def rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases[::-1]]
    return ''.join(bases)


def get_the_info(start, end, genome_name, info, target):

    genome_name = genome_name.partition('.')[0] 
    match = re.match(r".*?\/(.*)\.*",genome_name)
    genome_name = match.group(1)
    count=0
    t = 0 
    for i in range(len(info)):
        if t == 2:
            break
        if genome_name in info[i]:
            t = 2
            jo = ''.join(info[i:len(info)])
            tag = '^' + target +'.*'
            jo = jo.replace('\t','\n')
            for j in jo.split('\n'):
                f = re.search(tag, j)
                if f is not None:
                    count = count +1
     
    if count > 0:    
        target = target +'_' + str(count)

    
    alignment = target + '\t' + str(abs(int(end)-int(start))) + '\t' + start + '\t'+ end    
   
    head = genome_name +'\n'+ alignment  

    info = info.append(head)

    return info

def add_n_hits_to_info(my_16S, info):
    
    diz = {}
    
    for ele in my_16S:
        match = re.match(r".*?\+(.*)\n.*",ele)
        genome_name = match.group(1)
        count = 0
        for val in my_16S:
            if genome_name+'\n' in val:
                count +=1
        diz[genome_name] = count
        
    for i in range(len(info)):
        for genome_name in diz.keys():
            if genome_name in info[i]:
                
                p = info[i].split('\n')
                if len(p[1]) > 0:
                    info[i] = p[0] +'\t' + str(diz[genome_name]) + '\t' + p[1] + '\n'

    return info


def save_in_info(info):
    with open('My_16S_genes_%s_%s_info.txt' %(args.folder, query.replace(".", "")), 'w') as f:
        header = 'Genome_name' + '\t'+ 'n_16S_genes' +'\t'+ 'scaffold_code'+'\t'+'length'+'\t'+'start'+'\t'+ 'end' + '\n'
        f.write(header)
        for i in info:
            f.write(i)

x = open("nhmmer_output_%s_%s.txt" % (args.folder, query.replace(".", "")),'r')
result = x.read()
x.close()

my_16S = []
info = []

result= result.split('@') 

for i in range(1,len(result)):
    que = result[i]
    lines = que.split('\n')
    print(lines[0])
    if len(lines) < 9:
        print('No match found\n')
    else:
        if len(lines) > 9:
            print(str(len(lines)- 8) + ' matching hits were found.\n')
        n = 1
        for j in range (6,len(lines)-2):
                
            values = lines[j].split()
            
            fna = str(args.folder + '/' + lines[0])

            hmm_start = values[6]
            hmm_end = values[7]
            
            ali_start= values[9]
            ali_end= values[10]
            print('The match number '+ str(n) +' aligns from position: ' + ali_start + ' to position: ' + ali_end)
            n = n + 1
            target_name = values[0]
            print('The sequence I am looking for has ' + target_name + ' as target name.')
            gene_seq = get_the_gene(ali_start, ali_end, target_name, fna, values[8])
            print('The sequence of the 16S gene is:\n' + gene_seq + '\n')
            
            if gene_seq != 'not long enough':
                save_genes(gene_seq, target_name, fna, my_16S)
                get_the_info(ali_start, ali_end, fna, info, target_name)

save_in_Fasta(my_16S)
add_n_hits_to_info(my_16S, info)
save_in_info(info)

