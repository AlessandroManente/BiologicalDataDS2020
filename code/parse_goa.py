#!/usr/binx/env python

import gzip
import copy
from code import parse_go_obo
import os


def gen_block(f):
    """
    Genrator function that parses GOA GAF files (https://www.ebi.ac.uk/GOA/downloads)
    The generator yields a block of lines corresponding to the same protein
    UniProtKB       A0A024R1R8      hCG_2014768             GO:0002181      PMID:21873635   IBA     PANTHER:PTN002008372|SGD:S000007246     P       HCG2014768, isoform CRA_a       hCG_2014768     protein taxon:9606      20171102        GO_Central
    UniProtKB       A0A024RBG1      NUDT4B          GO:0003723      GO_REF:0000037  IEA     UniProtKB-KW:KW-0694    F       Diphosphoinositol polyphosphate phosphohydrolase NUDT4B NUDT4B  protein taxon:9606      20191109        UniProt
    UniProtKB       A0A024RBG1      NUDT4B          GO:0005829      GO_REF:0000052  IDA             C       Diphosphoinositol polyphosphate phosphohydrolase NUDT4B NUDT4B  protein taxon:9606      20161204        HPA
    """
    name, old_name = None, None
    chunk = []
    for line in f:
        line = line.decode()
        if line and line[0] != "!":
            _, name, _, _, term, _, ec, _, namespace, protein_name = line.split(
                "\t")[:10]
            if name != old_name and old_name:
                yield (
                    old_name, set(chunk)
                )  # return a set as there can be repetitions, i.e. the same term with different evidence codes
                chunk = []
            old_name = name
            chunk.append(term)
    # Last line
    if old_name:
        yield (old_name, set(chunk))


if __name__ == "__main__":

    # Get ontology data from parse_go_obo lib
    graph = parse_go_obo.parse_obo(os.getcwd() + "\\data\\function\\go.obo")
    ancestors, depth, roots = parse_go_obo.get_ancestors(graph)
    children = parse_go_obo.get_children(ancestors)

    # *** How many proteins are directly annotated with "regulation of kinase activity" (GO:0043549)
    #     and "mitochondrion" (GO:0005739)?
    proteins_counts = {}  # { term : number_of_proteins }
    with gzip.open(os.getcwd() + "\\data\\function\\goa_human.gaf.gz") as f:
        for acc, annotations in gen_block(f):
            for term in annotations:
                proteins_counts.setdefault(term, 0)
                proteins_counts[term] += 1
    print(graph["GO:0043549"]["def"], proteins_counts["GO:0043549"])
    print(graph["GO:0005739"]["def"], proteins_counts["GO:0005739"])

    # *** How many proteins are "regulation of kinase activity" (GO:0043549)?
    #     According to the "true path rule" you have to count direct annotations
    #     plus proteins annotated with children of that term

    proteins = {
    }  # { term : proteins_annotated_with_term } It contains proteins annotated directly with the term or its children
    with gzip.open(os.getcwd() + "\\data\\function\\goa_human.gaf.gz") as f:
        for acc, annotations in gen_block(f):
            # Copy direct annotations
            terms = copy.copy(annotations)
            # Add ancestors
            for term in annotations:
                terms.update(ancestors.get(term, set()))
            # For each term add protein accession to proteins dict
            for term in terms:
                proteins.setdefault(term, set()).add(acc)
    print(graph["GO:0043549"]["def"], len(proteins["GO:0043549"]))
    print(graph["GO:0005739"]["def"], len(proteins["GO:0005739"]))

    # *** Which are the 5 most abundant "biological process" terms in mitochondrial proteins (GO:0005739 mitochondrion)?
    terms_count = {
    }  # { term : count } count within the mitochondrial proteins set
    with gzip.open(os.getcwd() + "\\data\\function\\goa_human.gaf.gz") as f:
        for acc, annotations in gen_block(f):
            if acc in proteins["GO:0005739"]:
                # Copy direct annotations
                terms = copy.copy(annotations)
                # Add ancestors
                for term in annotations:
                    terms.update(ancestors.get(term, set()))
                # For each term add protein accession to proteins dict
                for term in terms:
                    terms_count.setdefault(term, 0)
                    terms_count[term] += 1

    # Sort by count and filter by biological_process namespace
    data = sorted([(k, v) for k, v in terms_count.items()],
                  key=lambda x: x[1],
                  reverse=True)
    for (k, v) in list(
            filter(lambda x: graph[x[0]]["namespace"] == "biological_process",
                   data))[:20]:
        print(k, v, graph[k]["def"])

    # *** Which are the top 20 most enriched terms in mitochondrial proteins (annotated with GO:0005739)?
    #     Measure the ratio between a term in mitochondrial proteins and in the rest of the human proteins
    #     and select those with the higher "fold-increase" (ratio)
    terms_set = {}  # { term : count }  mitochondrial proteins
    terms_rest = {}  #  { term : count }  other proteins
    proteins_set = 0  # number of mitochondrial proteins
    proteins_rest = 0  # number of remaining proteins
    with gzip.open(os.getcwd() + "\\data\\function\\goa_human.gaf.gz") as f:
        for acc, annotations in gen_block(f):
            # Copy direct annotations
            terms = copy.copy(annotations)
            # Add ancestors
            for term in annotations:
                terms.update(ancestors.get(term, set()))
            # For each term add protein accession to proteins dict

            if acc in proteins["GO:0005739"]:
                proteins_set += 1
                for term in terms:
                    terms_set.setdefault(term, 0)
                    terms_set[term] += 1
            else:
                proteins_rest += 1
                for term in terms:
                    terms_rest.setdefault(term, 0)
                    terms_rest[term] += 1

    data = []
    for term in terms_set:
        ratio_set = (terms_set[term] + 1) / proteins_set  # add pseudo count
        ratio_rest = terms_rest.get(term,
                                    1) / proteins_rest  # add pseudo count
        fold_increase = ratio_set / ratio_rest
        data.append(
            (term, terms_set[term], terms_rest.get(term,
                                                   0), ratio_set, ratio_rest,
             fold_increase, graph[term]["namespace"], graph[term]["def"]))
    for ele in sorted(data, key=lambda x: x[5], reverse=True)[:20]:
        print("{} {} {} {:.2g} {:.2g} {:.2g} {} {}".format(*ele))
