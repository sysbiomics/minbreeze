#!/usr/bin/env python

# Convert from BIOM with taxonomy into mothur-like output
# Input:
#    arg[1][tsv]: otu_table
#    arg[2][tsv]: tax_table
#    arg[3][tsv]: output file
#
    
import sys
import copy
import argparse
from operator import itemgetter
from functools import cmp_to_key

import pandas as pd
import numpy as np
from joblib import Parallel, delayed

REL_FLAG = False  # Relative abundance
TRIM_UNCLASS = False #  Trim unclassified.

otu_tab = sys.argv[1]
tax_tab = sys.argv[2]
output  = sys.argv[3]

otu_table = pd.read_csv(otu_tab, sep="\t", index_col=0)
tax_table = pd.read_csv(tax_tab, sep="\t", index_col=0)

print(otu_table)
print(tax_table)
# Make sure that otu_table and tax_table is joinable
#df.index.names = ['Date']

# Add tax column into OTU table (which only has hash)
#otu_table["total"] = otu_table.sum(axis=1)

addon = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]
newDf = []
for index, row in tax_table.iterrows():
    newRow = []
    prev = ""
    for a, item in zip(addon, row):
        if pd.isna(item):
            # Use prev name + unclassified
            newRow.append(prev + "_unclassified")
        else:
            prev = a + item
            newRow.append(prev)
    newDf.append(newRow)
    #tax_table.loc[index,] = newRow

tax_table = pd.DataFrame(newDf, columns=tax_table.columns, index=tax_table.index)

@cmp_to_key
def cmp_taxa(key1, key2):
    if key1 == "k__Unassigned" or key1 == "unclassified":
        return 1
    if key2 == "k__Unassigned" or key2 == "unclassified":
        return -1

    input = [key1, key2]
    # Lex sorted
    return int(sorted(input) == input)

class Tree(object):
    """ Tree where each node has value.
    """

    def __init__(self, name, value, parent):
        self.name = name
        self._children = {}
        self._parent = parent
        self.value = value

    def add(self, key, value):
        """ Add value to children
        """
        children = self._children
        if key in children:
            children[key].value += value
        else:
            children[key] = Tree(key, value, self)

    def insert(self, key, value):
        children = self._children
        if key in children:
            raise ValueError("Key already exists")
        else:
            children[key] = Tree(key, value, self)
        # Create a new tree if it is not available
        return children[key]

    def get(self, key):
        return self._children[key]

    def children(self):
        return list(self._children.values())

    def __getitem__(self, key):
        return self._children[key]

    def keys(self):
        return self._children.keys()


def create_tree():
    return Tree("root", 0, None)


def update_chain(keys, value, tree):
    # Chain update, when update value in leaf, or middle node, also
    # update every ancestral node.
    tree.value += value

    nextNode = tree
    for key in keys:
        # Update value of current node
        nextNode.add(key, value)
        nextNode = nextNode[key]

def normalize_tree(tree):
    """ Normalize value in tree to percentage
    """
    val = tree.value
    tree.value = 100 * tree.value / val
    for child in tree.children():
        t_normalize_tree(child, val)

def t_normalize_tree(node, denom):
    """
    """
    node.value = 100 * node.value / denom
    for child in node.children():
        t_normalize_tree(child, denom)


# TODO: Refactor this one

def tree_to_list(tree):
    return tree_to_list_recur(tree, rankID=[], output=[])

def tree_to_list_recur(node, rankID, output):
    """ Convert tree to list like the output from mothur. Deep first and append to output
    Args:
        node
    """
    # Collect taxa name and value at this stage
    if (TRIM_UNCLASS) and node.name.endswith("_unclassified"): # Don't collect unclassified tip,
        return

    taxlevel = len(rankID)
    rankIDs = ".".join(rankID)
    line_value = [(taxlevel, rankIDs, node.name), node.value]
    output.append(line_value)

    if len(node._children) == 0: # Stop
        return

    sortedKeys = sorted(node._children.keys())
    for i, key in enumerate(sortedKeys):
        nextRankID = rankID + [str(i)]
        tree_to_list_recur(node[key], nextRankID, output)

    return output

# Build a tree for each sample.
trees = {}
samples = otu_table.columns.values

# Join table for use
fulltab = otu_table.join(tax_table)
fulltab["taxa"] = fulltab.reindex(columns=["kingdom", "phylum", "class", "order", "family", "genus", "species"]).values.tolist()
fulltab = fulltab.drop(columns=["kingdom", "phylum", "class", "order", "family", "genus", "species"])

def create_tree_from_sample(fulltab, sample):
    """ Create taxonomy tree with expected abundance
    """
    tree = create_tree()
    sampledf = fulltab.reindex(columns= [sample,"taxa"] )
    for index, (abundance, taxa_lst) in sampledf.iterrows():
        update_chain(taxa_lst, abundance, tree)
    if REL_FLAG is True:
        normalize_tree(tree)
    return tree

# Good thing that joblib return thing in order
treevals = Parallel(n_jobs=10, backend="multiprocessing")(delayed(create_tree_from_sample)(fulltab, sample) for sample in samples)

tree = {}
for treename, treeval in zip(samples, treevals):
    trees[treename] = treeval
#for sample in samples:
#    print(sample)
#    tree = create_tree()
#    sampledf = fulltab.reindex(columns= [sample,"taxa"] )
#    for index, (abundance, taxa_lst) in sampledf.iterrows():
#        update_chain(taxa_lst, abundance, tree)
#    if REL_FLAG is True:
#        normalize_tree(tree)
#    trees[sample] = tree

# Convert list into series.
series = []
for name, tree in trees.items():
    # First three item are keys()
    lists = tree_to_list(tree)  # [(indices), list]
    indice = [i[0] for i in lists]
    values = [i[1] for i in lists]

    sIndex = pd.MultiIndex.from_tuples(indice)
    if REL_FLAG == True:
        taxaOtu = pd.Series(values, sIndex, dtype=np.float64)
    else:
        taxaOtu = pd.Series(values, sIndex, dtype=np.int64)
    taxaOtu.name = name
    series.append(taxaOtu)

df = pd.concat(series, axis=1)
df.index.set_names(["level", "taxtree", "tip"], inplace=True)

df.to_csv(output, sep="\t", header=True, index=True)


#
# DEBUG
#

"""
col = "total"
arr = otu_table[col]

tree = create_tree()
for otu, abundance in arr.iteritems():
    # Look for taxnomy
    # make sure that index is
    taxa = list(tax_table.loc[otu,])
    update_chain(taxa, abundance, tree)

tree_to_list(tree)

"""