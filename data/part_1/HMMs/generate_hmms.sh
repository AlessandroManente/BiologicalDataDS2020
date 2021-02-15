#!/bin/bash

#For now, in order for this to run you'll have to have a folder called "data" containing the SwissProt Database, and a folder called "binx" containing the hmmer source code, both in the same folder where the Git Repository is.
# 
#This script generates all the files needed for the HMM part. 
#
#Script Arguments:
# - The method of the MSA we want to take the MSA from: (C)offee, #Clustal-(O) or (M)uscle.
# - The Try number of the MSA.
# 
# In this way any MSA that follows our naming convention will be #uniquely determined.

 
echo "Choose MSA Method: (C)offee, Clustal-(O) or (M)uscle."
read msamethod
 
echo "Enter Try number"
read try

echo "Enter database used (swissprot, uniref90, uniref50, uniref100)"
read db


 # Build HMM model
 ../../../binx/hmmer-3.3.1/src/hmmbuild ./HMM_${msamethod}/hmm_model_${msamethod}_${try}_${db}.hmm ../MSAs/MSA_${msamethod}_${try}_${db}.fa
 
 # Command line commands to search SwissProt, using the created HMM model as query. This commands will output all 3 different formats of output.

../../../binx/hmmer-3.3.1/src/hmmsearch --tblout ./HMM_${msamethod}/hmmsearch_out_${msamethod}_${try}_${db}.tblout --domtblout ./HMM_${msamethod}/hmmsearch_out_${msamethod}_${try}_${db}.domtblout ./HMM_${msamethod}/hmm_model_${msamethod}_${try}_${db}.hmm ../../../data/uniprot_sprot.fasta > ./HMM_${msamethod}/hmmsearch_out_${msamethod}_${try}_${db}.align


