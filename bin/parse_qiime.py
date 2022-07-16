#!/usr/bin/env python

""" Parse QIIME2's taxonomy output into tab-delimited taxonomy table.
"""
import sys
import pandas as pd

def gg_parse(s):
    """ Parse taxonomy string in GG format. Return 7 levels of taxonomy.
       Args:
          s: taxonomy string. 
    """
    
    #TODO: Make it all lowercase
    abbr_dct = {"k": "kingdom", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}
    taxa_dct = {"kingdom": "", "phylum": "", "class": "", "order": "",
                "family": "", "genus": "", "species": ""}  # Because groupby exclude None value.
    items = s.split("; ")
   
    # Check
    if not s.startswith("k__"):
        # Probalbly cannot identify at all.
        return taxa_dct
    
    #if len(items) != 7:
    #    raise TaxaStringError()
        
    for token in items:
        abbrv, taxa = token.split("__")
        taxa_lvl = abbr_dct[abbrv]
        taxa = taxa if taxa else ""  # If empty, leave it as empty string
        # If it is bracket, then remove it
        if len(taxa) > 0 and taxa[0] == "[" and taxa[-1] == "]":
            taxa = taxa[1:-1]
        
        taxa_dct[taxa_lvl] = taxa
    
    return taxa_dct

def rdp_parse(s):
    """ Parse RDP taxonomy string with 7 level format (SILVA uses it.)
        D_0__Bacteria;D_1__Epsilonbacteraeota;D_2__Campylobacteria;D_3__Campylobacterales;D_4__Thiovulaceae;D_5__Sulfuricurvum;D_6__Sulfuricurvum sp. EW1
        The ambiguous_taxa will be convert to empty string.
        Args:
          s: String of taxnomy
    """
    abbr_dct = {"D_0": "kingdom", "D_1": "phylum", "D_2": "class", "D_3": "order",
                "D_4": "family", "D_5": "genus", "D_6": "species"}
    taxa_dct = {"kingdom": "", "phylum": "", "class": "", "order": "",
                "family": "", "genus": "", "species": ""}
    tokens = s.split(";")
    #if tokens[0] == "Ambiguous_taxa" or tokens[0] == "Unassigned":
    #    return taxa_dct
    for token in tokens: # D_0__Bacteria, or Ambiguous_taxa
        if token == "Ambiguous_taxa" or token == "Unassigned":
            return taxa_dct
        taxLv, taxName = token.split("__")
        # Make the output behave like GG parse
        taxLv = abbr_dct[taxLv]
        taxa_dct[taxLv] = taxName
        
    return taxa_dct

def qiime2_parse(s):
    """ Parse taxonomy string in qiime2 (later) format. Return 7 levels of taxonomy.
       Args:
          s: taxonomy string. 
    """
    
    #TODO: Make it all lowercase
    abbr_dct = {"d": "kingdom", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}
    taxa_dct = {"kingdom": "", "phylum": "", "class": "", "order": "",
                "family": "", "genus": "", "species": ""}  # Because groupby exclude None value.
    items = s.split("; ")
   
    # Check
    if not s.startswith("d__"):
        # Probalbly cannot identify at all.
        return taxa_dct
    
    #if len(items) != 7:
    #    raise TaxaStringError()
        
    for token in items:
        abbrv, taxa = token.split("__")
        taxa_lvl = abbr_dct[abbrv]
        taxa = taxa if taxa else ""  # If empty, leave it as empty string
        # If it is bracket, then remove it
        if len(taxa) > 0 and taxa[0] == "[" and taxa[-1] == "]":
            taxa = taxa[1:-1]
        
        taxa_dct[taxa_lvl] = taxa
    
    return taxa_dct


inp = sys.argv[1]
oup = sys.argv[2]
parser = sys.argv[3]

if parser == "silva":
    parser_fn = rdp_parse
elif parser == "gg":
    parser_fn = gg_parse
elif parser == "qiime2":
    parser_fn = qiime2_parse
else:
    raise ValueError("Parser need to be either silva or gg")

alignment = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]  

taxa_table = pd.read_csv(inp, sep="\t")
taxa_str = pd.DataFrame.from_records(taxa_table.Taxon.apply(parser_fn)).reindex(columns=alignment)
# First and last column are
#taxa_tab = pd.concat([taxa_table["Feature ID"], taxa_str, taxa_table["Confidence"]], axis=1)
taxa_tab = pd.concat([taxa_table.iloc[:, 0], taxa_str, taxa_table.iloc[:,-1]], axis=1)
taxa_tab = taxa_tab.set_index("Feature ID")

taxa_tab.to_csv(oup, sep="\t")
