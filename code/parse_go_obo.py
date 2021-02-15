#!/usr/bin/env python
import os


def parse_obo(obo_file):
    # Parse the ontology, exclude obsolete terms
    graph = {}  # { term_id : term_object }
    obj = {}  # { id: term_id, name: definition, is_a: list_of_parents, is_obsolete: True, namespace: namespace }
    with open(obo_file) as f:
        for line in f:
            line = line.strip().split(": ")
            if line and len(line) > 1:
                # print(line)
                k, v = line[:2]
                if k == "id" and v.startswith("GO:"):
                    obj["id"] = v
                elif k == "alt_id" and v.startswith("GO:"):
                    obj.setdefault("alt_id",[])
                    obj['alt_id'].append(v)
                elif k == "name":
                    obj["def"] = v
                elif k == "is_a" and v.startswith("GO:"):
                    obj.setdefault("is_a", []).append(v.split()[0])
                elif k == "is_obsolete":
                    obj["is_obsolete"] = True
                elif k == "namespace":
                    obj["namespace"] = v
            else:
                if obj.get("id") and not obj.get("is_obsolete"):
                    if "is_a" not in obj:
                        obj["is_root"] = True
                    graph[obj["id"]] = obj
                    # print(obj)
                obj = {}
    return graph


def get_ancestors(graph):
    """
    Build a dictionary of ancestors
    and calculate term depth (shortest path)
    """

    roots = set()
    for node in graph:
        if graph[node].get("is_root"):
            roots.add(node)

    depth = {}
    ancestors = {}  # { term : list_of_ancestor_terms }
    for node in graph:
        c = 0
        node_ancestors = []
        node_parents = graph[node].get("is_a")

        # Loop parents levels (parents of parents) until no more parents
        while node_parents:
            c += 1

            # Set root
            if node not in depth and roots.intersection(set(node_parents)):
                depth[node] = c

            # Add ancestors
            node_ancestors.extend(node_parents)

            # Update the list of parents (1 level up)
            node_parents = [term for parent in node_parents for term in graph[parent].get("is_a", [])]

        ancestors[node] = set(node_ancestors)
    return ancestors, depth, roots


def get_children(ancestors):
    children = {}  # { node : list_of_children }, leaf terms are not keys
    for node in ancestors:
        for ancestor in ancestors[node]:
            children.setdefault(ancestor, set()).add(node)
    return children


if __name__ == "__main__":

    # Gene ontology OBO file is available here
    # http://geneontology.org/docs/download-ontology/

    graph = parse_obo(os.getcwd() + "\\data\\function\\go.obo")
    ancestors, depth, roots = get_ancestors(graph)

    # *** How many different ancestors and leaf terms?
    ancestors_ids = set([ancestor for node in ancestors for ancestor in ancestors[node]])
    print(len(ancestors_ids))  # number of ancestors
    print(len(set(graph.keys()) - ancestors_ids))  # number of leaf terms

    # *** How many ancestors and leaf terms?
    children = get_children(ancestors)
    print(len(children))  # number of ancestors
    print(len(set(graph.keys()) - set(children.keys())))  # number of leaf terms

    # *** How many terms are children of GO:0016791?
    print(len([node for node in ancestors if "GO:0016791" in ancestors[node]]))  # number of children
    print(len(children["GO:0016791"]), graph["GO:0016791"]["def"])  # number of children

    # *** How many terms for each sub-ontology?
    for root in roots:
        print(root, graph[root], len(children[root]))

    # *** Which is the term with the largest number of ancestors?
    node = sorted([(len(ancestors[node]), node) for node in ancestors])[-1]
    print(node, graph[node[1]])

    # *** How many leaf terms for each sub-ontology?
    for root in roots:
        print(graph[root]["namespace"], len(set(children[root]) - set(children.keys())))  # number of leaf terms

    # *** How many terms at minimum depth 2 (1 node between the root and the term)
    print(len([node for node in depth if depth[node] == 2]))

    # *** How many terms for each branch at depth 1?
    print(len([node for node in depth if depth[node] == 1]))




