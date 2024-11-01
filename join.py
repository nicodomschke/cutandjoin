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
from pebble import concurrent, ProcessExpired
import hashlib
import itertools
import more_itertools as mit
import copy
from rdkit.Chem import AllChem
import os
import tqdm
import pandas as pd
import gurobipy as gp
from gurobipy import GRB
from util import *


client = get_db()

mol_db = client.cutandjoin

frag_col = mol_db.fragments
rule_col = mol_db.rules
mol_col = mol_db.molecules

order_dict = {
    rdkit.Chem.rdchem.BondType.SINGLE : 1,
    rdkit.Chem.rdchem.BondType.DOUBLE : 2,
    rdkit.Chem.rdchem.BondType.TRIPLE : 3,
}


solution_dict = {
    1.0 : rdkit.Chem.rdchem.BondType.SINGLE,
    2.0 : rdkit.Chem.rdchem.BondType.DOUBLE,
    3.0 : rdkit.Chem.rdchem.BondType.TRIPLE,
}


#get all reaction rules that apply to the set of generated fragments 
rules = rule_col.find({"$and":[ {"sat_left_len":{"$gt": 0}}, {"sat_right_len":{"$gt": 0}}]})

rulecnt = 0
for r in rules:
    rulecnt += 1
    print(rulecnt)

    id = r["_id"]
    print(id)

    smarts = r["rule"]
    l_sat_smiles = []
    r_sat_smiles = []

    # get all matching fragments
    for l in r["sat_left"]["target"]:
        for f in frag_col.find({"_id": l}):
            l_sat_smiles.append(f["smiles"])
            
    for l2 in r["sat_right"]["target"]:
        for f in frag_col.find({"_id": l2}):
            r_sat_smiles.append(f["smiles"])
    
    for _, x in tqdm.tqdm(enumerate(itertools.product(l_sat_smiles,r_sat_smiles)), total=len(l_sat_smiles)*len(r_sat_smiles)):
        # get both fragments as networkx graph
        v_graph, w_graph = smiles_to_nx(x[0]), smiles_to_nx(x[1])

        # Prepare variables to store data
        v_dict = {}
        w_dict = {}

        # Retrieve information about the cut-nodes from the cut-set of the first graph
        # in particular, store in a dictionary with 
        for i in v_graph.nodes(data=True):
            if i[1]["atomic_num"] == 92:                                                                            # if the atom is u - a placer holder for forming new edges (practically a free electron)
                neigh = list(v_graph.neighbors(i[0]))                                                               # Find all atoms belonging to this electron
                if len(neigh) > 1:                                                                                  # Sanity check. One electron - one atom
                    print("SOMETHINGS WRONG!!!!")                                       
                    input()
                if neigh[0] not in v_dict:                                                                          # If cut-node has no entry v_dict 
                    v_dict[neigh[0]] = [(i[0] ,order_dict[v_graph.get_edge_data(neigh[0], i[0])["order"]])]         # Add new entry with key = nodekey of cute node, value = list of tuples(u-id, bond order)
                else:                                                                                               # If cut-node has an entry
                    v_dict[neigh[0]] += [(i[0] ,order_dict[v_graph.get_edge_data(neigh[0], i[0])["order"]])]        # Add a new tuple to the existing list with (u-id, bond order)
            #v.remove_node(i[0])
        
        # Retrieve information about the nodes from the cut-set of the first graph - same as above
        for i in w_graph.nodes(data=True):
            if i[1]["atomic_num"] == 92:
                neigh = list(w_graph.neighbors(i[0]))
                if len(neigh) > 1:
                    print("SOMETHINGS WRONG!!!!")
                    input()
                if neigh[0] not in w_dict:
                    w_dict[neigh[0]] = [(i[0] ,order_dict[w_graph.get_edge_data(neigh[0], i[0])["order"]])]
                else:
                    w_dict[neigh[0]] += [(i[0] ,order_dict[w_graph.get_edge_data(neigh[0], i[0])["order"]])]
            #w.remove_node(i[0])
        
        # Extract a list with all cut-nodes C[V], C[W]
        vkeys = sorted(list(v_dict.keys()))
        wkeys = sorted(list(w_dict.keys()))

        # Prepare two new variables to store all tuples(u-id, bond-order) 
        # for each cut-set in one list
        v_list = []
        w_list = []
        for i in range(len(vkeys)):
            v = vkeys[i]
            v_list = v_list + v_dict[v]
        for i in range(len(wkeys)):
            w = wkeys[i]
            w_list = w_list + w_dict[w]
        
        # Sort the lists
        v_list = sorted(v_list)
        w_list = sorted(w_list)
        
        # Save a list with all u-ids
        v_id_list = [i[0] for i in v_list]
        w_id_list = [j[0] for j in w_list]

        # print("v")
        # print(v_dict)
        # print(v_list)
        # print(v_id_list)
        # print("------------------")
        # print("w")
        # print(w_dict)
        # print(w_list)
        # print(w_id_list)
        # input()

        # 1. Create Model
        no_solutions = GRB.MAXINT # Set limit for number of optimal solutions the solver provides
        model = gp.Model() # Create a new gurobi model
        # 2. Add variables, i. e. add matrix A for a_{ef}
        A = model.addVars(v_id_list, w_id_list, vtype=GRB.INTEGER, name="A")
        # 3. Add Constraints
        # 3.1 One Electron/electron pair only with one other electron/electron pair - constraint 1 (equation 2 in the paper)
        for c in v_list:
            model.addConstr(gp.quicksum(A[(c[0],d[0])] for d in w_list) == 1)
        for d in w_list:
            model.addConstr(gp.quicksum(A[(c[0],d[0])] for c in v_list) == 1)

        # 3.2 Labeling constraint - constraint 2 (quation 3 in the paper)
        # Remark: To avoid exceeding the number of permitted constraints, we substitued the 
        # constrained by taking the absolute value of the difference of the labels, in particular: \sum_{f€C_H} \sum_{f€G_G} |l_G(e)-l_H(e)|
        model.addConstr(gp.quicksum((abs(c[1]-d[1])) * A[(c[0],d[0])]  for d in w_list for c in v_list)==0)

        # 3.3 Inhibit multibindings - constraint 3 (equation 4 in the paper)
        for vkey, vvalue in v_dict.items():
            for wkey, wvalue in w_dict.items():
                model.addConstr(gp.quicksum(A[(cv[0], cw[0])] for cv in vvalue for cw in wvalue)<=1)

        # 4. Set Objective function (Equation 5)
        # Remark: We are using here the free electrons as 
        model.setObjective(gp.quicksum(gp.quicksum(gp.quicksum(gp.quicksum(d[1]*A[(c[0],d[0])] for d in w_dict[w]) for w in w_dict.keys()) for c in v_dict[v]) for v in v_dict.keys()), GRB.MAXIMIZE)

        # 5. Specify that multiple optimal solutions are provided by the solver
        model.Params.PoolSearchMode = 2

        # 6. Set number of provided solution to above assigned number 
        # (I guess you could restrict this number here in each of the joins to the number of joins that are maximally possible - the n! we talked about)
        model.Params.PoolSolutions = no_solutions

        # 7. Since gurobi is forgetting on our variables we need to tell him to remember
        model._vars = A

        # 8. Optimize the model

        model.optimize()
        #print("Number of solutions", model.SolCount)

        # 9. Extract the matchings
        # 9.1 Make a union of the two disjoints graph to add the newly formed edges into it
        union_copy = nx.union(v_graph, w_graph, rename=("v", "w"))
        # 9.2 Search each solution
        for i in range(model.SolCount):
            model.setParam(GRB.Param.SolutionNumber,i)                                                                                              # Set the solution number to get the exact solution
            union = union_copy.copy()                                                                                                               # Make a copy of the union to get a unique graph for each solution
            for v in model.getVars():                                                                                                               # Extract the matchings of the nodes
                ctuple = (v.varName.split("[")[1].split("]")[0].split(",")[0],v.varName.split("[")[1].split("]")[0].split(",")[1])                  # First: Extract matched nodes an write into new tuple
                matched = v.Xn                                                                                                                      # Second: Retrieve if atoms were matched or not 
                #print("ctuple", ctuple)
                #print(list(v_graph.neighbors(int(ctuple[0])))[0], list(w_graph.neighbors(int(ctuple[1])))[0], order)
                if abs(matched) == 0.0:                                                                                                             # If no matching, check next atom atom combination
                    continue
                
                getorder = ((int(ctuple[0]), list(v_graph.neighbors(int(ctuple[0])))[0]))                                                           # Get information on original tuple(cut-node, u-id)  
                order = v_graph[getorder[0]][getorder[1]]["order"]                                                                                  # Extract bond order
                union.add_edge(f"v{list(v_graph.neighbors(int(ctuple[0])))[0]}", f"w{list(w_graph.neighbors(int(ctuple[1])))[0]}", order = order)   # Add new edge for matched nodes with bond order

               

            # remove placeholders
            toremove = []
            for k in union.nodes(data=True):
                if k[1]["atomic_num"] == 92:
                    toremove.append(k[0])
            
            union.remove_nodes_from(toremove)
            # generate molecule SMILES
            smiles = getFragments(union)

            # save molecule to database
            insert_dic_mol = {}
            insert_dic_mol["_id"] = hashlib.md5(smiles.encode()).hexdigest()
            insert_dic_mol["smiles"] = smiles
            insert_dic_mol["inchi"] = Chem.MolToInchi(Chem.MolFromSmiles(smiles), options="-SNon")
            insert_dic_mol["generation"] = 1
            insert_dic_mol["cuts"] = []
            insert_dic_mol["synth"] = []

            try:
                mol_col.insert_one(insert_dic_mol)
            except pymongo.errors.DuplicateKeyError:
                continue
