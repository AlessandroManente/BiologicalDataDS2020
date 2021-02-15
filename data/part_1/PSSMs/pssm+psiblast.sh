#!/bin/bash

# For this step you will need:
#- Blast+ in folder "binx"
#	NCBI ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 
#or 
#	wget -P ../binx/ ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
#	tar -xf ../binx/ncbi-blast-2.11.0+-x64-linux.tar.gz -C ../binx/

# Create database
#../binx/ncbi-blast-2.11.0+/bin/makeblastdb -dbtype prot -in ../data/msa/uniprot_sprot.fasta

#-SwissProt Database in folder "data"
#	wget -P ../data/msa/ ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#	gunzip ../data/msa/uniprot_sprot.fasta.gz

echo "Choose MSA Method: (C)offee, Clustal-(O) or (M)uscle."
read msamethod
 
echo "Enter Try number"
read try

echo "Enter database used (swissprot, uniref90, uniref50, uniref100)"
read db

echo "How many psiblast iteration should I do?"
read iterations


#create PSSM model
../../../binx/ncbi-blast-2.11.0+/bin/psiblast -db ../../../data/msa/uniprot_sprot.fasta -in_msa ../MSAs/MSA_${msamethod}_${try}_${db}.fa -out_ascii_pssm ./ascii_pssm_${msamethod}_${try}_${db}.txt -out_pssm ./pssm_${msamethod}_${try}_${db}.pssm
# use PSSM model to make a psiblast query
../../../binx/ncbi-blast-2.11.0+/bin/psiblast -in_pssm ./pssm_${msamethod}_${try}_${db}.pssm -db ../../../data/msa/uniprot_sprot.fasta -outfmt 5 -num_iterations ${iterations} -evalue 0.01 > ./out_psiblast_${msamethod}_${try}_${db}_${iterations}iterations.xml