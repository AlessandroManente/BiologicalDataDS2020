# MSAs README

This folder contains the Multiple Sequence Alignments from the results of the BLAST search.

MSAs are performed on the following web services

- T-Coffee: <https://www.ebi.ac.uk/Tools/msa/tcoffee/>;
- MUSCLE: <https://www.ebi.ac.uk/Tools/msa/muscle/>;
- Clustal-Omega: <https://www.ebi.ac.uk/Tools/msa/clustalo/>.

Naming convention:

MSA_X_Y_Z.fa, where:

- X = Shortname of the method used: "C", "O", "M" for T-Coffee, Clustal-Omega and MUSCLE respectively.
- Y = Try number (NB: must coincide with the respective BLAST search's try number).
- Database = (swissprot, uniref90, uniref50, uniref100)
