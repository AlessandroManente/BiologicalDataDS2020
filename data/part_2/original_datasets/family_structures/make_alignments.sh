#!/bin/bash

# This Bash script will geneate all the pairwise sequence alignments that we need to do. It will output every file in the ./temp folder
echo "Enter model type (psiblast or hmm)"
read model

echo "If PSIBLAST, how many iterations{"
read iterations


echo "Enter MSA method (C, M or O)"
read msamethod

echo "Enter try number"
read try

echo "Enter database (swissprot, uniref90, uniref50 or uniref100)"
read db

if [ $model == 'psiblast' ]
then

directoryname=pdbs_${model}_${msamethod}_${try}_${db}_${iterations}iterations
else
directoryname=pdbs_${model}_${msamethod}_${try}_${db}
fi

for ent1 in ./${directoryname}/*.ent; do
	for ent2 in ./${directoryname}/*.ent; do
		#echo "${ent1}_${ent2}"
		TMalign ${ent1} ${ent2} > ./temp/$(basename ${ent1})_$(basename ${ent2}).out
	done
done
