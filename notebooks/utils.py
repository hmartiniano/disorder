#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import re
import itertools
import pathlib
import logging
from itertools import zip_longest

import requests
import requests_cache
requests_cache.install_cache('variant_cache')

import pandas as pd
import numpy as np

from pyfaidx import Fasta

#logging.basicConfig(level=logging.DEBUG)


AA_CODES = [['A', 'Ala', 'Alanine'],
            ['C', 'Cys', 'Cysteine'],
            ['D', 'Asp', 'Aspartic Acid'],
            ['E', 'Glu', 'Glutamic Acid'],
            ['F', 'Phe', 'Phenylalanine'],
            ['G', 'Gly', 'Glycine'],
            ['H', 'His', 'Histidine'],
            ['I', 'Ile', 'Isoleucine'],
            ['K', 'Lys', 'Lysine'],
            ['L', 'Leu', 'Leucine'],
            ['M', 'Met', 'Methionine'],
            ['N', 'Asn', 'Asparagine'],
            ['P', 'Pro', 'Proline'],
            ['Q', 'Gln', 'Glutamine'],
            ['R', 'Arg', 'Arginine'],
            ['S', 'Ser', 'Serine'],
            ['T', 'Thr', 'Threonine'],
            ['V', 'Val', 'Valine'],
            ['W', 'Trp', 'Tryptophan'],
            ['Y', 'Tyr', 'Tyrosine']]

aa_codes = pd.DataFrame(AA_CODES, 
                        columns=["IUPAC_aa_code", "aa_code", "aminoacid"],
                       ).set_index("IUPAC_aa_code")

conv = aa_codes["aa_code"].to_dict() #["aa_code"]

inv_conv = {v:k for k, v in conv.items()}


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def convert_to_hgvsp(protein_id, mutation):
    old, pos, new = re.findall(r"([A-Z])(\d+)([A-Z])", mutation)[0] 
    new_representation = protein_id + ":p." + conv[old] + pos + conv[new]
    logging.info(f"Converting {mutation} to {new_representation}") 
    return new_representation


def get_gene_id(symbol):

    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/human/{symbol}"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    return r.json()[0]["id"]


def lookup_gene_id(gene_id):

    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{gene_id}"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    return r.json()


def get_protein_sequence(gene_id):

    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{gene_id}?type=protein;species=human"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    return r.json()



def get_gene(symbol):
    """Query ensembl for gene and retrieve protein sequence for the canonical transcript."""
    gene_id = get_gene_id(symbol)
    gene_data = lookup_gene_id(gene_id)
    seq = get_protein_sequence(gene_data["canonical_transcript"].split(".")[0])
    gene_data["protein"] = seq
    return gene_data


def get_variants(gene, biotype="protein_coding"):
    """deprecated"""
    gene_id = gene["id"]
    server = "https://rest.ensembl.org"
    ext = f"/overlap/id/{gene_id}?feature=variation&species=homo_sapiens&so_term=SO:0001583"
    if biotype is not None:
        ext += f"&biotype={biotype}"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    print(r.url)
    return pd.DataFrame(r.json())


def get_variants(gene):
    """ Retrieve missense variants overlaping a protein coding region"""
    protein_id = gene["protein"]["id"]
    server = "https://rest.ensembl.org"
    ext = f"/overlap/translation/{protein_id}?feature=transcript_variation;species=homo_sapiens&type=missense_variant" #&so_term=SO:0001583"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    print(r.url)
    return pd.DataFrame(r.json())


def recode(variant):
    """Recode variants for the annotation step"""
    server = "https://rest.ensembl.org"
    ext = f"/variant_recoder/human/{variant}?"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    print(repr(decoded))
    alleles = list(decoded[0].keys())
    logging.info(f"Alleles: {alleles}")
    return decoded[0][alleles[0]]['hgvsc'][0]


def process_consequences(data):
    df = pd.DataFrame(data)
    df = df.dropna(subset=["colocated_variants", "transcript_consequences"])
    df = df.explode("colocated_variants")
    df = df.explode("transcript_consequences")
    df2 = df.set_index("input").transcript_consequences.apply(pd.Series)
    df = pd.merge(df.drop(columns=["transcript_consequences"]), df2, left_on="input", right_index=True) 
    df = df[df.biotype == "protein_coding"]
    df = df[~df.canonical.isna()]
    df2 = df.set_index("input").colocated_variants.apply(pd.Series)
    df = df[df.biotype == "protein_coding"]
    df = df[~df.canonical.isna()]
    df = pd.merge(df.drop(columns=["colocated_variants"]), df2, left_on="input", right_index=True) 
    df = df.drop_duplicates(subset=["start_x", "end_x", "allele_string_x"])
    #df = df.drop_duplicates()
    return df


def get_consequences(variant):

    server = "https://rest.ensembl.org"
    ext = f"/vep/human/hgvs/{variant}?CADD=1&LoF=1&protein=1&uniprot=1&hgvs=1&canonical=1"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    df = process_consequences(r.json())
    return df


def get_consequences_bulk(variants):

    server = "https://rest.ensembl.org"
    ext = "/vep/human/id?CADD=1&LoF=1&protein=1&uniprot=1&hgvs=1&canonical=1"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    result = []
    for grp in grouper(variants, 100):
        grp = ", ".join(['"{0}"'.format(x) for x in grp if x is not None])
        print(grp)
        r = requests.post(server+ext, headers=headers, data=f'{{ "ids" : [ {grp} ] }}')

        print(r.url)
        if not r.ok:
          r.raise_for_status()
          sys.exit()

        df = process_consequences(r.json())
        result.append(df)
    result = pd.concat(result)
    return result


def write_fasta(name, seq, directory="."):
    if not pathlib.Path(directory).exists():
        os.mkdir(directory)
    fname = f"{directory}/{name}.fasta"
    with open(fname, "w") as f:
        f.write(f">{name}\n")
        for group in grouper(seq, 60): 
            seq = "".join([i for i in group if i is not None])
            f.write(f"{seq}\n")
    return fname



def convert_from_hgvsp(mutation):
    try:
        old, pos, new = re.findall(r"ENSP[0-9]+\.[0-9]+\:p\.([a-zA-Z]+)(\d+)([a-zA-Z]+)", mutation)[0] 
    except:
        print(f"could not convert {mutation}")
        return None
    old = inv_conv[old]
    new = inv_conv[new]
    return old + pos+ new



def annotate_variants(variants, additional_variants, hgvsp=False):
    """"""
    protein_id = variants.seq_region_name.iloc[0]
    variant_data = []
    csq = get_consequences_bulk(variants["id"].tolist())
    variant_data.append(csq)
    for variant in additional_variants:
        r = convert_to_hgvsp(protein_id, variant)
        print(r)
        r = recode(r)
        print(r)
        print(30 * "#")
        csq = get_consequences(r)
        variant_data.append(csq)
    variant_data = pd.concat(variant_data)
    variant_data = variant_data[~variant_data.hgvsp.str.contains("=")]
    variant_data["mutation"] = variant_data.hgvsp.apply(convert_from_hgvsp)
    variant_data = variant_data.set_index("mutation")
    variant_data = variant_data.drop_duplicates(subset=["hgvsp"])
    return variant_data


def mutate(sequence, mutation, hgvsp=False):
    if hgvsp:
        old, pos, new = re.findall(r"ENSP[0-9]+\.[0-9]+\:p\.([a-zA-Z]+)(\d+)([a-zA-Z]+)", mutation)[0] 
        old = inv_conv[old]
        new = inv_conv[new]
    else:
        old, pos, new = re.findall(r"([A-Z])(\d+)([A-Z])", mutation)[0] 
    pos = int(pos)
    ndx = pos - 1
    logging.info(f"{mutation} : mutating {old} on position {pos} to {new}")
    assert sequence[ndx] == old
    mutated = list(sequence)
    mutated[ndx] = new
    mutated = "".join(mutated)
    logging.info(sequence)
    logging.info(" " * ndx + "*")
    logging.info(mutated)
    return mutated



