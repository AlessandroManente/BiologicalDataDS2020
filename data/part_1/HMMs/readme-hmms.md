# HMM README

This folder contains the hmm models generated from our MSAs. For each MSA try we generate:

- A HMM model;
- The output of our HMM search on Swissprot, using such HMM as input. The output is in three different formats: align, tblout, domtblout.

  Steps to create HMMs and do the search are packed inside the "generate_hmms" bash script.

To use it you will have to download Swissprot:

`wget -P ../data/msa/ ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz gunzip ../data/msa/uniprot_sprot.fasta.gz`

------ Naming Conventions -----

- HMM MODEL: hmm_model_X_Y;
- ALIGN OUTPUT: hmmsearch_out_X_Y.align;
- TBLOUT: hmmsearch_out_X_Y.tblout;
- DOMTBLOUT: hmmsearch_out_X_Y.domtblout;

For all of those, it must be:

- X = Shortname for the MSA method: "C", "O", "M" for T-Coffee, Clustal-Omega and MUSCLE respectively.
- Y = try number.
