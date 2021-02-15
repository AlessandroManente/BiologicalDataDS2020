# Description of the folder

Here we will keep only the baseline datasets for part 2. (?)

To generate **"family_sequences"** dataset for a given model, go here: https://www.uniprot.org/uploadlists/ and paste in the search bar all the UNIPROT Accession Codes found by the model.
Doing this will allow us to map those UniProt codes to Uniref90 database, which contains "clusterized" versions of those sequences.
Specifically, if 2 or more sequences are very similar, they will be grouped inside a single cluster; that cluster inside Uniref90 will be considered as a single database entry. The name for the entry will be a representative protein of the cluster, that will give the name to the Uniref90 entry (i.e.it will be like "UniRef90_D0V3Y4").

To generate **"family_structures"** dataset, check the notebook!
