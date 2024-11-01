from rdkit import Chem
import pymongo
import hashlib
import itertools
import more_itertools as mit
import copy
from util import get_db, smiles_to_nx



client = get_db()


# Substitute for your db of choice.
mol_db = client.cutandjoin

# collection names can be altered, but should be altered across all files
# to avoid missing data.
frag_col = mol_db.fragments
rule_col = mol_db.rules
mol_col = mol_db.molecules

# Get all entries of the fragment collection
res = frag_col.find()

# Counter to output number of molecules while they are being read from the fatabase
counter = 0

# List where smiles that produce errors of some sort are collected.
debuglist =  []

# Mapping between rdkit bond types and an integer representation
order_dic = {
    Chem.rdchem.BondType.SINGLE:1,
    Chem.rdchem.BondType.DOUBLE:2,
    Chem.rdchem.BondType.TRIPLE:3,
}

for r in res:
    if counter % 1000 == 0:
        print("FINISHED UPDATING", counter, "MOLECULES")
    counter += 1

    # Parse smiles from the database to networkx graph
    smiles = r["smiles"]
    graph = smiles_to_nx(smiles)

    # Extract the cut-nodes that represent open bonds of a fragment
    target_nodes = [x[0] for x in graph.nodes(data=True) if x[1]["atomic_num"] == 92]

    # Dictionary that collects the open bonds for each node. The open bonds for each atom
    # are derived by the number of uranium atoms it is connected to. The uranium atoms are a place-holder
    # for open bonds.
    ends = {}

    # For each uranium we check which atoms it is connected to. We keep track of the nodes it is connected to
    # in the ends dictionary. We also keep track of the corresponding bond, i.e. edge weight in order to produce
    # proper natural joins.
    for u in target_nodes:

        # iterate over the neighbors of u. val stores the correponding edge with attributes
        for n, val in graph[u].items():

            # If have already encountered n, we add a corresponding value to the dictionary (additional bond)
            if n in ends.keys():
                try:
                    ends[n] += [order_dic[val["order"]]]
                except:
                    print(f"failed: {smiles}")
                    debuglist.append(smiles)

            # Else we initialize a new entry for n.
            else:
                try:
                    ends[n] = [order_dic[val["order"]]]
                except KeyError:
                    print(f"failed: {smiles}")
                    debuglist.append(smiles)

    # Parse the dictionary into an array for saving it to the db. For example [(1,1,1), (1,2)] represents a fragment
    # that resulted from cutting three single bonds (1,1,1) at one atom and a single bond and a double bond (1,2) at
    # another atom. We sort this array in order to identify fragments with the same open-end-configuration easily.
    # (sorting a list of tuples in python can be done canonically via the sorted function).
    bond_array = sorted([tuple(ends[x]) for x in ends.keys()])

    # Add additional database attribute that contains the sorted bonds, bond configuration ("open bonds on cut-vertices)
    # and a string representation of the bond config (for human-readable data when browsing the mongodb).
    frag_col.update_one({"_id": r["_id"]}, {"$set":{"sorted_bonds":sorted(r["bonds"]), "bond_config":bond_array, "bond_config_string": str(bond_array)}})


# Get all (updated) entries of the fragment collection
res = frag_col.find()

# Output any potentially faulty entries, or entries that cause rdkit errors
with open("./add_ligate_info_debug.txt", mode="w") as of:
    for i in debuglist:
        of.write(i + "\n")

def make_sidechains(rest, bonds, docking_partners):
    '''
    Function to create a SMARTS reaction left side for a fragment. Rest is any pre-existing reaction string. the bonds
    list describes the "open" bonds each fragment has. Docking partners is the molecular placeholder for a fragment.
    Parameters
    ----------
    rest any pre-existing reaction string. Can be left empty
    bonds list that contains the "open" bonds of a fragment
    docking_partners list or corresponding placeholder atoms that are used in a fragment.

    Returns
    -------
    String of SMARTS reaction
    '''

    # Carry over the preexisting reaction string
    sidechain = rest

    # For each bond add a "requirement" to the left side. After finish, left side for a fragment with two
    # open single bonds a a single open double bond would be: (-U)(-U)(=U).
    for b in bonds:
        if b == 1:
            sidechain += "(-" + docking_partners[bonds.index(b)] + ")"
        if b == 2:
            sidechain += "(=" + docking_partners[bonds.index(b)] + ")"
        if b == 3:
            sidechain += "(#" + docking_partners[bonds.index(b)] + ")"
    return sidechain


def get_valid_mappings(map_array_0, map_array_1):
    '''
    Returns all combinations of map_array_0 - map_array_1  tuples without replacement.
    :param map_array_0: list of indices
    :param map_array_1: list of indices
    :return: list of tuples that contains all combinations without replacement from map_array_0 x map_array_1
    '''

    if len(map_array_0) >= len(map_array_1):
        comb_list = [list(zip(x, map_array_1)) for x in itertools.permutations(map_array_0, len(map_array_1))]
    else:
        comb_list = [list(zip(x, map_array_0)) for x in itertools.permutations(map_array_1, len(map_array_0))]
        comb_list = [[(t[1], t[0]) for t in m] for m in comb_list]

    clist = set([tuple(sorted(t)) for t in comb_list])
    comb_list = clist

    return comb_list

def bond_representation(combination, index_to_bond_0, index_to_bond_1):
    '''
    Translates a tuple from index to corresponding bond-config. Allows for checking for equivalent mappings
    :param combination: set of tuples of indices that correspond to keys in index_to_bond_0 and _1
    :param index_to_bond_0: dictionary mapping indexto bond-config
    :param index_to_bond_1: dictionary mapping indexto bond-config
    :return: translated tuple
    '''
    new_rep = tuple((index_to_bond_0[c[0]], index_to_bond_1[c[1]]) for c in combination)
    return new_rep

def bond_representation_with_bonds(combination, index_to_bond_0, index_to_bond_1):
    '''
    Translates a tuple from index to corresponding bond-config. Also contains bond-type.
    :param combination: set of tuples of indices that correspond to keys in index_to_bond_0 and _1
    :param index_to_bond_0: dictionary mapping indexto bond-config
    :param index_to_bond_1: dictionary mapping indexto bond-config
    :return: translated tuple
    '''
    new_rep = tuple((c[0], (index_to_bond_0[c[1][0]], index_to_bond_1[c[1][1]])) for c in combination)
    return new_rep

def check_mapping(add, index_to_bond_0, index_to_bond_1, bond_num, mapping_check=None):
    '''
    Check if a given mapping is valid
    :param add: mapping between bonds
    :param index_to_bond_0: dictionary translating index to bond_type
    :param index_to_bond_1: dictionary translating index to bond_type
    :param bond_num: type of bond that is matching
    :param mapping_check: dictionary to stepwisely check mapping
    :return: True if mapping is valid, False else
    '''
    add = [t[1] for t in add]

    if mapping_check == None:
        mapping_check = {}
        mapping_check["l"] = {}
        mapping_check["r"] = {}
        for t in add:
            mapping_check["l"][t[0]] = {"cap": {n: index_to_bond_0[t[0]].count(n) for n in set(index_to_bond_0[t[0]])}, "map": {n:[] for n in set(index_to_bond_0[t[0]])}}
            mapping_check["r"][t[1]] = {"cap": {n: index_to_bond_1[t[1]].count(n) for n in set(index_to_bond_1[t[1]])}, "map": {n: [] for n in set(index_to_bond_1[t[1]])}}
    else:
        for t in add:
            if t[0] not in mapping_check["l"].keys():
                mapping_check["l"][t[0]] = {"cap": {n: index_to_bond_0[t[0]].count(n) for n in set(index_to_bond_0[t[0]])},
                                            "map": {n: [] for n in set(index_to_bond_0[t[0]])}}
            if t[1] not in mapping_check["r"].keys():
                mapping_check["r"][t[1]] = {"cap": {n: index_to_bond_1[t[1]].count(n) for n in set(index_to_bond_1[t[1]])},
                                            "map": {n: [] for n in set(index_to_bond_1[t[1]])}}


    for t in add:
        mapping_check["l"][t[0]]["map"][bond_num] += [t[1]]
        mapping_check["r"][t[1]]["map"][bond_num] += [t[0]]

        if len(mapping_check["l"][t[0]]["map"][bond_num]) > mapping_check["l"][t[0]]["cap"][bond_num] \
                or len(mapping_check["l"][t[0]]["map"][bond_num]) > len(set(mapping_check["l"][t[0]]["map"][bond_num])):
            return False

        if len(mapping_check["r"][t[1]]["map"][bond_num]) > mapping_check["r"][t[1]]["cap"][bond_num] \
                or len(mapping_check["r"][t[1]]["map"][bond_num]) > len(set(mapping_check["r"][t[1]]["map"][bond_num])):

            return False

        if len([True for n in mapping_check["l"][t[0]]["map"] if t[1] in mapping_check["l"][t[0]]["map"][n]]) > 1:
            return False

        if len([True for n in mapping_check["r"][t[1]]["map"] if t[0] in mapping_check["r"][t[1]]["map"][n]]) > 1:
            return False

    return mapping_check

# Select all distinct bond configurations in the database 
res = frag_col.aggregate([{"$group":{"_id":{"bonds":"$sorted_bonds"},"counts":{"$sum":1}}}])

bonds_list = []

for r in res:
    if r["_id"]["bonds"] != None:
        bonds_list += [r["_id"]["bonds"]]

# Generate all possible combinations of bond-configuration, e.g. [(1,), (1,), (2,), (2,)], [(1,2), (1,2)], ...
# for every different bond array [1,1,2,2] we find in the database. Since this covers all combinations of open-bonds
# that are compatible, we generate all these matching configurations of these and create their reaction smarts.
for bonds_list_element in bonds_list:

    print()
    print("********** CURRENT BOND TYPES: ", bonds_list_element, " **********")
    print()
    n = len(bonds_list_element)

    # Get the bond array
    bond_array_original = []
    bond_array_original = sorted(bonds_list_element)

    # Ids for left side mapping for SMARTS reaction string
    rests_0 = ["[*:" + str(i) + "]" for i in range(1, n+1)]
    rests_1 = ["[*:" + str(i) + "]" for i in range(n+1, 2*n + 1)]


    final_combinations = []

    for i in range(1, n + 1):

        # First create all combinations of bond-configs for the current bond types
        combinations = []

        # All "open bonds" are on distinct atoms
        if i == 1:
            added_combinations = [[(k,) for k in bond_array_original]]
            final_combinations += added_combinations
            continue
        # All bonds are on the same atom
        elif i == n:
            added_combinations =[[tuple(bond_array_original)]]
            final_combinations += added_combinations
            continue
        # All other possibilities to distribute the open bonds in bond-array onto i different atoms
        else:
            elems = [k for k in range(1, n + 1 - (i - 1))]
            combinations = [c for c in itertools.combinations_with_replacement(elems, i) if sum(c) == n]

        # Make sure added_combinations is empmty
        added_combinations = []

        # Post processing for other cases
        for c in combinations:
            bond_array = bond_array_original

            combination_collection = [[] for k in range(len(c))]

            combination_identity = 0

            for num_bp in c:
                thiscomb = mit.distinct_combinations(bond_array, num_bp)

                combination_collection[combination_identity] += thiscomb
                combination_identity += 1


            product = list(itertools.product(*combination_collection))

            duplicate_list = []

            for p in product:
                flattened = sorted([item for sublist in p for item in sublist])
                if flattened != bond_array_original:
                    continue
                else:
                    if set(p) not in duplicate_list:
                        duplicate_list += [set(p)]
                        # print(p)
                        added_combinations += [list(p)]

        # Collect the combinations
        final_combinations += added_combinations

    # Loop for checking if a natural join between two combinations of bonds in final_combinations
    # is possible. If a natural join is possible we calculate every distinct possibility of joining
    # the two fragments. Note that we only compute the distinct possibilities based on bond-type.
    # For example, (1,)-(1,) (2,)-(2,) (2,)-(2,) only has one distinct natural join since the different
    # (2,) matchings are the same since they only refer to the rests in the smart-strings which
    # will match the different 2 bonds when plugged into a SMARTS-reaction library like rdkit's
    for i in range(len(final_combinations)):
        bonding0 = final_combinations[i]

        rxn_str0 = "("

        # Generate left side SMARTS for mol0 by using sidechains
        for k in range(len(bonding0)):
            s0 = make_sidechains(rests_0[k], bonding0[k], ["[U]" for l in range(len(bonding0[k]))])
            rxn_str0 += s0 + "."

        rxn_str0 = rxn_str0[0:-1] + ")"

        # find longest bond configuration for simple compatibility check
        longest_0 = len(max(bonding0, key=lambda x: len(x)))

        for j in range(i, len(final_combinations)):
            bonding1 = final_combinations[j]

            rxn_str1 = "("

            # Generate left side SMARTS for mol1 by using sidechains
            for k in range(len(bonding1)):
                s1 = make_sidechains(rests_1[k], bonding1[k], ["[U]" for l in range(len(bonding1[k]))])
                rxn_str1 += s1 + "."

            rxn_str1 = rxn_str1[0:-1] + ")"

            longest_1 = len(max(bonding1, key=lambda x: len(x)))

            # Baseline filter, if we have one rest with n open ends, we need at least n different atoms on the other side
            if longest_0 > len(bonding1) or longest_1 > len(bonding0):
                continue

            bond_array_uniq_sort = sorted(list(set(bond_array_original)))

            # Derive mapping subcomponents. Tuples containing ones can only be mapped to a tuple containing a one.
            # Bond mapping array collects those components for both partners.
            bond_mapping_array_0 = [[t for t in bonding0 if e in t] for e in set(bond_array_original)]
            bond_mapping_array_1 = [[t for t in bonding1 if e in t] for e in set(bond_array_original)]

            # Array that contains the corresponding bond type that is mapped at each index of bond mapping array
            bond_indices = list(set(bond_array_original))

            for b in range(len(bond_mapping_array_0)):

                # Only contains the 1s/2s/... of the corresponding bond mapping array.
                # We precheck here if a mapping is actually possible.
                clean_bond_mapping_array_0 = [[ee for ee in e if ee == bond_array_uniq_sort[b]] for e in bond_mapping_array_0[b]]
                clean_bond_mapping_array_1 = [[ee for ee in e if ee == bond_array_uniq_sort[b]] for e in bond_mapping_array_1[b]]

                # Check what the longest tuple is for a bond type. For example, we cannot naturally join (1,1,1,1,1) if the other
                # array contains only three tuples.
                longest_b_0 = len(max(clean_bond_mapping_array_0, key=lambda x: len(x)))
                longest_b_1 = len(max(clean_bond_mapping_array_1, key=lambda x: len(x)))

                no_mapping = False

                # Simple compatibility checks
                if longest_b_0 > len(bond_mapping_array_1[b]):
                    no_mapping = True
                    break
                if longest_b_1 > len(bond_mapping_array_0[b]):
                    no_mapping = True
                    break

            if no_mapping:
                continue

            # Mappings for bond->id, id->rest and index->bond for dealing with the different levels of mapping.
            # Index in bond-array to string of smarts identifier, i.e. 0:[1:*], 1:[2:*] etc for each molecule
            index_to_rest_0 = {}
            index_to_rest_1 = {}

            # Bond-config type to indices of the bonding0/1 arrays
            bond_to_index_0 = {}
            bond_to_index_1 = {}

            # smarts identifiers to the bond-cofig they refer to
            index_to_bond_1 = {}
            index_to_bond_0 = {}

            # Simply building the mappings
            for b in range(len(bonding0)):
                index_to_rest_0[b] = rests_0[b]
                try:
                    bond_to_index_0[bonding0[b]] += [b]
                except KeyError:
                    bond_to_index_0[bonding0[b]] = [b]
                index_to_bond_0[b] = bonding0[b]

            for b in range(len(bonding1)):
                index_to_rest_1[b] = rests_1[b]
                try:
                    bond_to_index_1[bonding1[b]] += [b]
                except KeyError:
                    bond_to_index_1[bonding1[b]] = [b]
                index_to_bond_1[b] = bonding1[b]

            # In bond combinations we will store every mol0 bond-config with valid mappings to bonds in mol1
            bond_combinations = {}
            for b in range(len(bond_indices)):
                bond_combinations[bond_indices[b]] = {bond:[] for bond in bond_mapping_array_0[b]}

            # bond_mapping_array contains for each bond-type i.e. 1,2,... a list of the bonds-configuration that match
            # this type. Aditionally, each bond-configuration appears as many times as bonds this configuration
            # can match to. For example, (1,1,2) would appear twice in the first element of bond_mapping_array since it
            # has to accept 2 new 1-bonds.
            for b in range(len(bond_mapping_array_0)):

                # Bonds
                bond_set_0 = list(set(bond_mapping_array_0[b]))
                bond_set_1 = list(set(bond_mapping_array_1[b]))

                # Count of bond-type
                bond_counts_0 = [bond_mapping_array_0[b].count(t) for t in bond_set_0]
                bond_counts_1 = [bond_mapping_array_1[b].count(t) for t in bond_set_1]

                # Tuples conf and count of conf
                counted_bond_set_0 = [(t, bond_mapping_array_0[b].count(t)) for t in bond_set_0]
                counted_bond_set_1 = [(t, bond_mapping_array_1[b].count(t)) for t in bond_set_1]

                # Mapping conf count of conf
                counted_bond_map_0 = {t: bond_mapping_array_0[b].count(t) for t in bond_set_0}
                counted_bond_map_1 = {t: bond_mapping_array_1[b].count(t) for t in bond_set_1}

                for m in range(len(bond_set_0)):

                    lpivot = bond_set_0[m]
                    left = []
                    right = []

                    # We collect in left all the indices corresponding to the current bond-config.
                    # Indices appear multiple times if bond-config can be matched multiple times.
                    for e in bond_to_index_0[lpivot]:
                        for c in range(lpivot.count(bond_indices[b])):
                            left += [e]

                    # We collect in right all the indices of bonds that match the current bond-config.
                    # Indices appear multiple times if bond-config can be matched multiple times.
                    for r in set(bond_mapping_array_1[b]):
                        for e in bond_to_index_1[r]:
                            for x in range(r.count(bond_indices[b])):
                                right += [e]

                    right = sorted(right)

                    # Comb gives us all combination of left-right tuples
                    # without replacement. Combinations contatins a set
                    # of tuples that are a mapping from left to right
                    comb = get_valid_mappings(left, right)

                    remove_combs = set()

                    # Remove illegal combinations, i.e. mappings where we map two atoms twice
                    # like (1,2) -= (1,2), hence forming a double and a single bond between a
                    # single atom
                    for bond1 in set(bond_mapping_array_1[b]):
                        for c in comb:
                            indices = bond_to_index_1[bond1]
                            for index in indices:
                                matches = [t[0] for t in c if t[1] == index]

                                # If we have double matches, we remove them.
                                if len(set(matches)) != len(matches):

                                    # We dont remove them immediately since we iterate over combs
                                    remove_combs.add(c)
                                    break

                    # Remove illegal mappings after checking
                    for c in remove_combs:
                        comb.remove(c)

                    # Store the mappings
                    bond_combinations[bond_indices[b]][index_to_bond_0[left[0]]] = list(comb)

            # We now have computed for each different bond-config in mol0 all individually valid mappings to bonds in
            # mol1. We now need to plug the individually valid mappings together such that the resulting mapping
            # corresponds to a natural join.

            # bond_combinations still contains a lot of mappings that are the same if you interpret every (2,) (2,) as
            # equivalent. Since when using SMARTS with rdkit it matches every fitting matching molecular components to these
            # (2,)-structures corresponding to (=U) in the final SMARTS string we can consider each mapping between
            # multiple of the same bond-config as equal. This loop removes them by first translating the
            # indices in bond_combinations to bond_types and then checking them for duplicates.
            for b in bond_combinations.keys():
                for bond in bond_combinations[b].keys():

                    if len(bond) != 1:
                        continue

                    cleaned_mapping_array = copy.copy(bond_combinations[b][bond])
                    already_removed = []

                    for o1 in range(0, len(bond_combinations[b][bond])):

                        if o1 in already_removed:
                            continue
                        e = bond_combinations[b][bond][o1]
                        rep_e = sorted(bond_representation(e, index_to_bond_0, index_to_bond_1))


                        for o2 in range(o1+1, len(bond_combinations[b][bond])):

                            if o2 in already_removed:
                                continue

                            e2 = bond_combinations[b][bond][o2]
                            rep_e2 = sorted(bond_representation(e2, index_to_bond_0, index_to_bond_1))

                            if rep_e2 == rep_e:
                                cleaned_mapping_array.remove(e2)
                                already_removed += [o2]

                    bond_combinations[b][bond] = cleaned_mapping_array


            # list-indices for all mappings for each bond type. A valid combination is an element of each list in
            # decision array such that the bond-mapping it corresponds to is a natural join.
            decision_array = [[b for b in range(len(bond_combinations[t][bond_type]))] for t in bond_combinations.keys() for bond_type in bond_combinations[t]]

            # bond type for which each bond-config is matched
            decision_bonds = [(t, bond_type) for t in bond_combinations.keys() for bond_type in bond_combinations[t]]

            # Sorting the decision array by length of the combinations.
            decision_indices = [x for x in range(len(decision_array))]
            decision_indices = sorted(decision_indices, key=lambda x: len(decision_array[x]))
            decision_array = [decision_array[x] for x in decision_indices]
            decision_bonds = [decision_bonds[x] for x in decision_indices]

            # All possible combinations to plug the bond-configs together
            decision_tuples = list(itertools.product(*decision_array))

            valid_tuples = {}

            final_combination_collection = []

            # Infer which mappings of bonds is actually a natural join. We save in valid_tuples
            # which combination of mappings we already checked at some point to reduce some overhead.
            for x in range(len(decision_tuples)):

                current_tuple = []
                current_mapping = []
                mapping_check = None

                # stepwise add a additional mapping
                for y in range(len(decision_tuples[x])):
                    current_tuple += [decision_tuples[x][y]]

                    # check if the mapping combination is checked already
                    try:
                        if valid_tuples[tuple(current_tuple)] == False:
                            current_mapping = []
                            mapping_check = None
                            break
                        else:
                            mapping_check = copy.deepcopy(valid_tuples[tuple(current_tuple)])
                            current_mapping += [(decision_bonds[y][0], t) for t in bond_combinations[decision_bonds[y][0]][decision_bonds[y][1]][decision_tuples[x][y]]]
                            continue
                    except KeyError:
                        pass

                    # Check if mapping is valid.
                    res = check_mapping([(decision_bonds[y][0], t) for t in bond_combinations[decision_bonds[y][0]][decision_bonds[y][1]][decision_tuples[x][y]]],
                                  index_to_bond_0, index_to_bond_1, decision_bonds[y][0], mapping_check=mapping_check)

                    # If not valid then save false it in valid tuples
                    if res == False:
                        valid_tuples[tuple(current_tuple)] = False
                        current_mapping = []
                        mapping_check = None
                        break
                    # If valid then we extend the current mapping and keep going
                    else:
                        valid_tuples[tuple(current_tuple)] = copy.deepcopy(res)
                        mapping_check = res
                        current_mapping += [(decision_bonds[y][0], t) for t in bond_combinations[decision_bonds[y][0]][decision_bonds[y][1]][decision_tuples[x][y]]]

                # If final mapping is valid we keep it.
                if current_mapping != []:
                    final_combination_collection += [tuple(current_mapping)]


            cleaned_final_combination_collection = copy.deepcopy(final_combination_collection)

            already_removed = []

            # Duplicate removal for equal mappings. Can unfortunately occur again after adding all combinations.
            for x in range(len(final_combination_collection)):

                if x in already_removed:
                    continue

                c0 = final_combination_collection[x]
                rep0 = sorted(bond_representation_with_bonds(c0, index_to_bond_0, index_to_bond_1))

                for y in range(x+1, len(final_combination_collection)):

                    if y in already_removed:
                        continue

                    c1 = final_combination_collection[y]
                    rep1 = sorted(bond_representation_with_bonds(c1, index_to_bond_0, index_to_bond_1))

                    if rep0 == rep1:
                        cleaned_final_combination_collection.remove(c1)
                        already_removed += [y]

                        continue

            final_combination_collection = cleaned_final_combination_collection


            # Create rdkit mol object in order to parse it into the right side SMARTS.
            # The idea is to create a atom-mapped molecule where the atom-map indices correspond to
            # the indices of the left side of the SMARTS for mol0 and mol1. Bonds are then introduced to
            # the molecule such that they correspond to the natural join edges between mol0 and mol1. This atom-
            # mapped molecule is then parsed into a SMARTS string which models the natural join.
            for e in final_combination_collection:

                # Output the mapping so that a human can read it
                prnt_str = "MAPPING (type, mol0, mol1):\n"
                for t in e:
                    prnt_str += "(|" + str(t[0]) + "|, " + str(index_to_bond_0[t[1][0]]) + ", " + str(
                        index_to_bond_1[t[1][1]]) + "),\n "
                print(prnt_str)

                mol_graph = Chem.RWMol()

                # Keep track of rdkit internal indices
                index_to_atom_0 = dict()
                index_to_atom_1 = dict()

                for t in e:
                    # Add atoms for matches in bond0
                    if index_to_rest_0[t[1][0]] not in index_to_atom_0.keys():
                        atom = Chem.Atom("S")
                        atom.SetAtomMapNum(int(index_to_rest_0[t[1][0]].split(":")[-1][0:-1]))
                        index = mol_graph.AddAtom(atom)
                        index_to_atom_0[index_to_rest_0[t[1][0]]] = index
                    # Add atoms for matches in bon1
                    if index_to_rest_1[t[1][1]] not in index_to_atom_1.keys():
                        atom = Chem.Atom("S")
                        atom.SetAtomMapNum(int(index_to_rest_1[t[1][1]].split(":")[-1][0:-1]))
                        index = mol_graph.AddAtom(atom)
                        index_to_atom_1[index_to_rest_1[t[1][1]]] = index

                    # Add bond with corresponding type
                    if t[0] == 1:
                        mol_graph.AddBond(index_to_atom_0[index_to_rest_0[t[1][0]]], index_to_atom_1[index_to_rest_1[t[1][1]]], order=Chem.rdchem.BondType.SINGLE)
                    elif t[0] == 2:
                        mol_graph.AddBond(index_to_atom_0[index_to_rest_0[t[1][0]]], index_to_atom_1[index_to_rest_1[t[1][1]]], order=Chem.rdchem.BondType.DOUBLE)
                    elif t[0] == 3:
                        mol_graph.AddBond(index_to_atom_0[index_to_rest_0[t[1][0]]], index_to_atom_1[index_to_rest_1[t[1][1]]], order=Chem.rdchem.BondType.TRIPLE)
                    elif t[0] == 12:
                        mol_graph.AddBond(index_to_atom_0[index_to_rest_0[t[1][0]]],
                                          index_to_atom_1[index_to_rest_1[t[1][1]]], order=Chem.rdchem.BondType.AROMATIC)
                    else:
                        continue

                # Parse right side to SMARTS
                output = Chem.MolToSmarts(mol_graph.GetMol()).replace("#16","*")

                # Optional output
                print("SMARTS RIGHT     ", output)
                print("RXN SMARTS       ", rxn_str0, ".", rxn_str1, ">>", output, sep="")
                print()

                # Create a proper SMARTS Rule
                rule_string = rxn_str0 + "." + rxn_str1 + ">>(" + output + ")"

                # Compute hash for insert into db
                rule_hash = hashlib.md5(rule_string.encode()).hexdigest()

                # Create db representation
                insert_dic = {
                    "_id":rule_hash,
                    "rule":rule_string,
                    "left":rxn_str0,
                    "right":rxn_str1,
                    "left_bonds":sorted(bonding0),
                    "left_bonds_str":str(sorted(bonding0)),
                    "right_bonds":sorted(bonding1),
                    "right_bonds_str":str(sorted(bonding1)),
                    "product":output,
                    "sat_left":{},
                    "sat_right":{}
                }

                # Try to insert rule. DuplicateKeyError means we already inserted this rule at some point
                # since then the hashes clash.
                try:
                    rule_col.insert_one(insert_dic)
                except pymongo.errors.DuplicateKeyError:
                    pass

                # Separator for output
                print("\n---------------------------------\n")



rules = rule_col.find()

# For each rule we find matching fragments (left and right side)
for rule in rules:

    print("PROCESSING RULE ", rule["_id"])
    print("PROCESSING RULE ", rule["rule"])

    bonds_left = rule["left_bonds_str"]
    bonds_right = rule["right_bonds_str"]

    print("TRYING LEFT ", bonds_left, " - ", type(bonds_left))
    print("TRYING RIGHT ", bonds_right, " - ", type(bonds_right))

    # Find fragments with bond configurations that match the left and right side of the rule
    left_matches = frag_col.find({"bond_config_string":bonds_left})
    right_matches = frag_col.find({"bond_config_string":bonds_right})

    # Initialize potentially already existing sat-left and sast-right
    left_sat = rule["sat_left"]
    left_sat.update({"target": []})
    right_sat = rule["sat_right"]
    right_sat.update({"target": []})

    added_left = False
    added_right = False

    # Try to add left-side satisfying fragments for mol0 individually, only if they are not contained in the list yet.
    for l in left_matches:

        if l["_id"] in left_sat["target"]:
            continue
        else:
            added_left = True
            left_sat["target"] += [l["_id"]]

    # Try to add left-side satisfying fragments for mol1 individually, only if they are not contained in the list yet.
    for r in right_matches:

        if r["_id"] in right_sat["target"]:
            continue
        else:
            added_right = True
            right_sat["target"] += [r["_id"]]

    update_dic = {}

    left_len = sum([len(left_sat[k]) for k in left_sat.keys()])
    right_len = sum([len(right_sat[k]) for k in right_sat.keys()])

    print("LEN LEFT ", left_len)
    print("LEN RIGHT ", right_len)

    # Check if we updated anything
    if added_left:
        update_dic["sat_left"] = left_sat
        update_dic["sat_left_len"] = left_len
    else:
        update_dic["sat_left_len"] = left_len

    # Check if we updated anything
    if added_right:
        update_dic["sat_right"] = right_sat
        update_dic["sat_right_len"] = right_len
    else:
        update_dic["sat_right_len"] = right_len

    # Update db
    if added_left or added_right:
        rule_col.update_one({"_id": rule["_id"]}, {"$set":update_dic})

    print("\n--------------------\n")
