import networkx as nx
import matplotlib.pyplot as plt
import random as rnd
import rdkit
from rdkit import Chem
import logging
import time
from networkx.algorithms.connectivity import EdgeComponentAuxGraph
import random
import logging
from itertools import permutations, product, chain
import pymongo
from pebble import concurrent as pebble_concurrent
from pebble import ProcessExpired
import hashlib
import itertools
import more_itertools as mit
import copy
from rdkit.Chem import AllChem
import os
import tqdm
import pandas as pd
import concurrent.futures
import argparse
import sys
from util import *

#get client
client = get_db() 

# Substitute for your db of choice.
mol_db = client.cutandjoin

# collection names can be altered, but should be altered across all files
# to avoid missing data.
frag_col = mol_db.fragments
rule_col = mol_db.rules
mol_col = mol_db.molecules

argparser = argparse.ArgumentParser(description="Open SMILES file")
argparser.add_argument("filepath", help="Path to the SMILES file")

args = argparser.parse_args()

# read SMILES input file
try:
    with open(args.filepath, 'r') as inputfile:
        smiles = inputfile.read().splitlines()
except FileNotFoundError:
    print(f"File not found")
    sys.exit()
except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit()

for cur_smiles in smiles:

    #setting up database entries
    G_orig = smiles_to_nx(cur_smiles)

    insert_dic_mol = {}
    insert_dic_mol["_id"] = hashlib.md5(cur_smiles.encode()).hexdigest()
    insert_dic_mol["smiles"] = cur_smiles
    insert_dic_mol["inchi"] = Chem.MolToInchi(Chem.MolFromSmiles(cur_smiles), options="-SNon")
    insert_dic_mol["generation"] = 0
    insert_dic_mol["cuts"] = []
    insert_dic_mol["synth"] = []

    insert_frags = []

    copy_G = G_orig.copy()

    bridges = list(nx.bridges(G_orig))
    single_cuts = [[x] for x in bridges]
    # bridge removal
    G_orig.remove_edges_from(bridges)
    connected_components = nx.connected_components(G_orig)

    cutsetsG = []
    cutsetsG += single_cuts

    # get all cuts of remaining scaffolds
    for con in connected_components:
        subg = nx.induced_subgraph(G_orig, con)
        cutsetsG += list_all_cutsets(subg)

    cutsetsG = sorted(cutsetsG, key=lambda x: len(x))

    insert_dic_mol["cuts"] = [c for c in cutsetsG if len(c) > 0]

    G_orig = copy_G.copy()

    # try getting all resonance structures of supplied SMILES
    rdk0 = Chem.ResonanceMolSupplier(
        Chem.MolFromSmiles(cur_smiles), Chem.KEKULE_ALL
    )

    kekule_set_0 = list()

    future = iter_kekul(rdk0)

    try:
        kekule_set_0 = future.result()
    except TimeoutError as error:
        print("Kekulization took longer than %d seconds" % error.args[1])
        kekule_set_0 = fallback_kekule(cur_smiles)
    except ProcessExpired as error:
        print("%s. Exit code: %d" % (error, error.exitcode))
        kekule_set_0 = fallback_kekule(cur_smiles)
    except Exception as error:
        print("Kekulization raised %s" % error)
        print(error.traceback)
        kekule_set_0 = fallback_kekule(cur_smiles)

    for x in cutsetsG:

        G_orig = copy_G.copy()
        if x == []:
            continue

        for k in kekule_set_0:
            # print("GOING OVER KEKULIZATION")
            G_orig = mol_to_nx(k)

            if G_orig == None:
                continue

            cut = []

            # apply cuts to molecule graph, get fragments with placeholder node(s)
            for nodes, attributes in G_orig.nodes(data=True):
                attributes["is_aromatic"] = False

            placeholders = [i + len(G_orig.nodes()) + 1 for i in range(2 * len(x))]

            attr = {
                "atomic_num": 92,
                "formal_charge": 0,
                "chiral_tag": rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
                "hybridization": rdkit.Chem.rdchem.HybridizationType.SP2,
                "num_explicit_hs": 0,
                "is_aromatic": False,
            }

            jk = 0

            bonds = []

            for cutty in range(len(x)):
                G_orig.add_edge(
                    x[cutty][0],
                    placeholders[jk],
                    order=G_orig.get_edge_data(x[cutty][0], x[cutty][1])["order"],
                )
                G_orig.add_edge(
                    x[cutty][1],
                    placeholders[jk + 1],
                    order=G_orig.get_edge_data(x[cutty][0], x[cutty][1])["order"],
                )

                bonds += [G_orig.get_edge_data(x[cutty][0], x[cutty][1])["order"]]

                G_orig.nodes[placeholders[jk]].update(attr)
                G_orig.nodes[placeholders[jk + 1]].update(attr)

                G_orig.remove_edge(x[cutty][0], x[cutty][1])
                jk += 2

            frags = [
                G_orig.subgraph(c).copy() for c in nx.connected_components(G_orig)
            ]

            # try to convert fragment to SMILES
            for frag in frags:
                try:
                    s = getFragments(frag)
                except rdkit.Chem.rdchem.KekulizeException as e:
                    with open("./debug_db.txt", mode="a") as dgo:
                        print("failed to kekulize")
                        dgo.write(str(x) + " " + cur_smiles + "\n")
                        dgo.write(str(e) + "\n\n")
                frag_dic = {
                    "_id": hashlib.md5(s.encode()).hexdigest(),
                    "smiles": s,
                    "origin": [{insert_dic_mol["_id"]: [list(y) for y in x]}],
                    # "target": [target],
                    "ends": len(x),
                    "bonds": bonds,
                }
                insert_frags += [frag_dic]
                # print("FRAGMENT ", s)

    print("INSERT STEP.")

    # try saving fragment to DB
    try:
        mol_col.insert_one(insert_dic_mol)
    except pymongo.errors.DuplicateKeyError:
        old_mol = mol_col.find({"_id": insert_dic_mol["_id"]})[0]
        update_target_mol = []

    for f in insert_frags:
        try:
            frag_col.insert_one(f)
        except pymongo.errors.DuplicateKeyError:
            old_frag = frag_col.find({"_id": f["_id"]})[0]

            update_origin = []
            update_target = []

            if f["origin"][0] in old_frag["origin"]:
                pass
            else:
                update_origin += old_frag["origin"]
                update_origin += f["origin"]

            if update_target == [] and update_origin == []:
                continue
            elif update_target == [] and update_origin != []:
                frag_col.update_one({"_id": f["_id"]}, {"$set": {"origin": update_origin}})
                frag_col.update_one({"_id": f["_id"]},
                                    {"$set": {"target": update_target, "origin": update_origin}})



