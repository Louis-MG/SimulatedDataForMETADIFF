#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = '1.0.0'
__author__ = 'Alban Mathieu'
__date__ = '01-12-2020'

import sys
import argparse
import textwrap
#import pandas as pd
#import requests as rq
from pathlib import Path
import re
from collections import defaultdict
#from bs4 import BeautifulSoup as bs
#from Bio import SeqIO
#####################################
def warn(msg):
    """Utils: print a warning
    """
    print(f"\033[93m{msg}\033[0m")
##########################################################################
tree = lambda: defaultdict(tree)
hash_all = tree()

hash_all2 = tree()
hash = defaultdict(lambda :defaultdict(defaultdict))
node2parent = defaultdict(lambda :defaultdict(defaultdict))
node2type = defaultdict(lambda :defaultdict(defaultdict))
tax2name = defaultdict(lambda :defaultdict(defaultdict))
name2tax = defaultdict(lambda :defaultdict(defaultdict))
merged2real = defaultdict(lambda :defaultdict(defaultdict))
##########################################################################
def get_nodes(nodes = None):
    if nodes is None:
        raise Exception("nodes is empty.")
    if not nodes.is_file():
        raise Exception(f"nodes.dmp file '{nodes}' does not exist.")
    with nodes.open('rt') as f:
        for line in f:
            node = line.split("\t")[0]
            parent = line.split("\t")[2]
            type = line.split("\t")[4]
            node2parent[node] = parent
            node2type[node] = type

def get_names(names = None):
    if names is None:
        raise Exception("names is empty.")
    if not names.is_file():
        raise Exception(f"names.dmp file '{names}' does not exist.")
    with names.open('rt') as f:
        for line in f:
            l2 = line.strip()
            tax = l2.split("\t")[0]
            name = l2.split("\t")[2]
            name2tax[name] = tax
            if l2.split("\t")[6] == "scientific name":
                tax2name[tax] = name

def get_merged(merged = None):
    if merged is None:
        raise Exception("merged is empty.")
    if not merged.is_file():
        raise Exception(f"merged.dmp file '{merged}' does not exist.")
    with merged.open('rt') as f:
        for line in f:
            l2 = line.strip()
            merge = l2.split("\t")[0]
            real = l2.split("\t")[2]
            merged2real[merge] = real

##########################################################################
def read_list(list = None):
    """Read list file.
    """
    report = []
    if list is None:
        raise Exception("list is empty.")
    if not list.is_file():
        raise Exception(f"list file '{list}' does not exist.")
    with list.open('rt') as f:
        for line in f:
            report.append(line.strip())
    return report

##########################################################################
sample_list = []
taxonomy_list = []
def read_sample_data(report = None):
    """Read input sample file.
    """
    report2 = None
    if report is None:
        raise Exception("report is empty.")
    for line in report:
        report2 = Path(line)
        if not report2.is_file():
            raise Exception(f"report2 file '{report2}' does not exist.")
        psamlev = report2.name #permet de retrouver le nom du fichier sans le path
        samlev = re.findall(r"(^.*)_bracken_(S)_report.txt", psamlev)
        sample = samlev[0][0]
        print(sample)
        #level = samlev[0][1]
        if sample not in sample_list:
            sample_list.append(sample)
        with report2.open('rb') as f:
            for li in f:
                data = li.decode("utf-8").strip()
                val = data.split("\t")[1]
                lev2 = data.split("\t")[3]
                tax = data.split("\t")[4]
                if lev2 == "S":
                    #hash_all[tax][sample][lev2] = int(val)
                    taxonomy = get_tax(taxi=tax, taxonomy_list = taxonomy_list)
                    sk = taxonomy[0]
                    try:
                        hash["SK"][sk][sample] += int(val)
                    except KeyError:
                        hash["SK"][sk][sample] = int(val)
                    phy = ";".join(taxonomy[:2])
                    try:
                        hash["P"][phy][sample] += int(val)
                    except KeyError:
                        hash["P"][phy][sample] = int(val)
                    cla = ";".join(taxonomy[:3])
                    try:
                        hash["C"][cla][sample] += int(val)
                    except KeyError:
                        hash["C"][cla][sample] = int(val)
                    ord = ";".join(taxonomy[:4])
                    try:
                        hash["O"][ord][sample] += int(val)
                    except KeyError:
                        hash["O"][ord][sample] = int(val)
                    fam = ";".join(taxonomy[:5])
                    try:
                        hash["F"][fam][sample] += int(val)
                    except KeyError:
                        hash["F"][fam][sample] = int(val)
                    gen = ";".join(taxonomy[:6])
                    try:
                        hash["G"][gen][sample] += int(val)
                    except KeyError:
                        hash["G"][gen][sample] = int(val)
                    spe = ";".join(taxonomy)
                    try:
                        hash["S"][spe][sample] += int(val)
                    except KeyError:
                        hash["S"][spe][sample] = int(val)
    return hash

##########################################################################

def write_result(otu = None,spe = None,gen = None,cla = None,fam = None,ord = None,phy = None,kin = None, hash = None):
    """Write end results to file
    """
    if otu is None:
        raise Exception("otu is empty.")
    if spe is None:
        raise Exception("spe is empty.")
    if gen is None:
        raise Exception("gen is empty.")
    if fam is None:
        raise Exception("fam is empty.")
    if cla is None:
        raise Exception("cla is empty.")
    if ord is None:
        raise Exception("ord is empty.")
    if phy is None:
        raise Exception("phy is empty.")
    if kin is None:
        raise Exception("kin is empty.")
    if hash is None:
        raise Exception("hash is empty.")
    #if hash_all2 is None:
    #    raise Exception("hash_all2 is empty.")
    #out_all = Path(f"{output}")
    #with out_all.open("w") as oua:
    #    sam = "\t".join(sample_list)
    #    oua.write(f"all_level\t{sam}\n")
    #    for tax, hash_allb in hash_all2.items():
    #        tabline = []
    #        tabline.append(tax)
    #        for i in sample_list:
    #            if i in hash_allb.keys():
    #                tabline.append(hash_allb[i])
    #            else:
    #                tabline.append(0)
    #        if sum(tabline[1:]) > 0:
    #            string_line = [str(do) for do in tabline]
    #            linef = "\t".join(string_line) #une autre idée \t".join(str(do) for do in tabline)
    #            oua.write(f"{linef}\n")
    for lev2, uhash in hash.items():
        out2 = None
        out3 = None
        if lev2 == "S":
            out2 = Path(f"{otu}")
            out3 = Path(f"{spe}")
        elif lev2 == "G":
            out2 = Path(f"{gen}")
        elif lev2 == "F":
            out2 = Path(f"{fam}")
        elif lev2 == "O":
            out2 = Path(f"{ord}")
        elif lev2 == "C":
            out2 = Path(f"{cla}")
        elif lev2 == "P":
            out2 = Path(f"{phy}")
        elif lev2 == "SK":
            out2 = Path(f"{kin}")
        else:
            print(f"error in parsing {lev2}")
        with out2.open("w") as ou:
            sam = "\t".join(sample_list)
            ou.write(f"{lev2}\t{sam}\n")
            for tax, uhash2 in uhash.items():
                ou.write(f"{tax}")
                for i in sample_list:
                    try:
                        ou.write(f"\t{uhash2[i]}")
                    except KeyError:
                        ou.write(f"\t0")
                ou.write(f"\n")
        if lev2 == "S":
            with out3.open("w") as ou3:
                sam = "\t".join(sample_list)
                ou3.write(f"{lev2}\t{sam}\n")
                for tax, uhash2 in uhash.items():
                    tax2 = re.sub(';\d+$', '', tax)
                    ou3.write(f"{tax2}")
                    for i in sample_list:
                        try:
                            ou3.write(f"\t{uhash2[i]}")
                        except KeyError:
                            ou3.write(f"\t0")
                    ou3.write(f"\n")

##########################################################################


def get_tax(taxi = None, taxonomy_list = None):
    """get all taxonomy using names and nodes files
    """
    if taxi is None:
        raise Exception("taxi is empty.")
    parent = None
    #taxonomy_list = listi
    if taxi != "1":
        if node2parent[taxi]:
            parent = node2parent[taxi]
        else:
            print(f"No parent or no node {taxi}\n")
        if node2type[taxi] == "species":
            taxonomy = tax2name[taxi]
            taxonomy_list = ["other","other","other","other","other","other","other"]
            taxonomy_list.append(taxi)
            taxonomy_list[6] = "s__" + taxonomy
        elif node2type[taxi] == "genus":
            taxonomy = tax2name[taxi]
            #print(tax2name[taxi])
            #print(taxi)
            #print(taxonomy_list)
            taxonomy_list[5] = "g__" + taxonomy
        elif node2type[taxi] == "family":
            taxonomy = tax2name[taxi]
            taxonomy_list[4] = "f__" + taxonomy
        elif node2type[taxi] == "order":
            taxonomy = tax2name[taxi]
            taxonomy_list[3] = "o__" + taxonomy
        elif node2type[taxi] == "class":
            taxonomy = tax2name[taxi]
            taxonomy_list[2] = "c__" + taxonomy
        elif node2type[taxi] == "phylum":
            taxonomy = tax2name[taxi]
            taxonomy_list[1] = "p__" + taxonomy
        elif node2type[taxi] == "superkingdom":
            taxonomy = tax2name[taxi]
            taxonomy_list[0] = "k__" + taxonomy
        else:
            taxonomy_list = taxonomy_list
        get_tax(parent, taxonomy_list)
    return taxonomy_list


##########################################################################
def get_argparser():
    """Create the argument parser
    """
    parser = argparse.ArgumentParser(prog = "Parse and filter sam alignement",
                                     add_help=True,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""\
    extract ID samples occurences and report number of sequences mapped
    """))
    # attr
    parser.add_argument("-l", "--list", type=str, help="A path to the report list", required=True)
    parser.add_argument("-o", "--otu",type=str, help="A path to otu level results file", required=True)
    parser.add_argument("-S", "--species",type=str, help="A path to species level results file", required=True)
    parser.add_argument("-G", "--genus",type=str, help="A path to genus level results file", required=True)
    parser.add_argument("-F", "--family",type=str, help="A path to family level results file", required=True)
    parser.add_argument("-C", "--classes",type=str, help="A path to class level results file", required=True)
    parser.add_argument("-O", "--order",type=str, help="A path to order level results file", required=True)
    parser.add_argument("-P", "--phylum",type=str, help="A path to phylum level results file", required=True)
    parser.add_argument("-K", "--kingdom",type=str, help="A path to kingdom level results file", required=True)
    parser.add_argument("-na", "--names",type=str, help="A path to names.dmp", required= True)
    parser.add_argument("-no", "--nodes",type=str, help="A path to nodes.dmp", required=True)
    parser.add_argument("-me", "--merged",type=str, help="A path to merged.dmp", required=True)
    return parser




def main(*args, **kwargs):
    """Main function
    """
    # create result dir
    # res_dir = write_result(kwargs["output"])
    # exec
    get_names(names=kwargs["names"])
    get_nodes(nodes=kwargs["nodes"])
    get_merged(merged=kwargs["merged"])
    queryf = read_list(list=kwargs["list"])
    hash = read_sample_data(report=queryf)
    #hash_all2 = create_all_res(hash=hash)
    write_result(  hash= hash,
                   otu=kwargs["otu"],
                   spe=kwargs["species"],
                   gen=kwargs["genus"],
                   fam=kwargs["family"],
                   cla=kwargs["classes"],
                   ord=kwargs["order"],
                   phy=kwargs["phylum"],
                   kin=kwargs["kingdom"])
    # print(f"fasta selection '{kwargs['query']}' done.")

if __name__ == '__main__':
    # build arg parser, top level parser
    parser = get_argparser()

    # parse arguments
    parsed_args, _ = parser.parse_known_args()

    # args and kwargs
    kw = {
        "list": Path(parsed_args.list) if parsed_args.list is not None else None,
        "otu": Path(parsed_args.otu) if parsed_args.otu is not None else None,
        "species": Path(parsed_args.species) if parsed_args.species is not None else None,
        "genus": Path(parsed_args.genus) if parsed_args.genus is not None else None,
        "family": Path(parsed_args.family) if parsed_args.family is not None else None,
        "classes": Path(parsed_args.classes) if parsed_args.classes is not None else None,
        "order": Path(parsed_args.order) if parsed_args.order is not None else None,
        "phylum": Path(parsed_args.phylum) if parsed_args.phylum is not None else None,
        "kingdom": Path(parsed_args.kingdom) if parsed_args.kingdom is not None else None,
        "names": Path(parsed_args.names) if parsed_args.names is not None else None,
        "nodes": Path(parsed_args.nodes) if parsed_args.nodes is not None else None,
        "merged": Path(parsed_args.merged) if parsed_args.merged is not None else None,
    }
    # exec
    main(**kw)
