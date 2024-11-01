import networkx as nx
import rdkit
from rdkit import Chem
import logging
from pebble import concurrent as pebble_concurrent
from pebble import ProcessExpired
from rdkit.Chem import AllChem
import pymongo

def getFragments(G):
    '''

    Parameters
    ----------
    G networkx graph

    Returns
    -------
    smiles parsed from the networkx graph
    '''

    # Initialize rdkit molecule.
    mol = Chem.RWMol()

    # Parse atom attributes.
    atomic_nums = nx.get_node_attributes(G, "atomic_num")
    chiral_tags = nx.get_node_attributes(G, "chiral_tag")
    formal_charges = nx.get_node_attributes(G, "formal_charge")
    node_is_aromatics = nx.get_node_attributes(G, "is_aromatic")
    node_hybridizations = nx.get_node_attributes(G, "hybridization")
    num_explicit_hss = nx.get_node_attributes(G, "num_explicit_hs")
    node_to_idx = {}

    # Initialize and add each atom with corresponding attributes to molecule.
    for node in G.nodes():
        a = Chem.Atom(atomic_nums[node])
        a.SetChiralTag(chiral_tags[node])
        a.SetFormalCharge(formal_charges[node])
        a.SetIsAromatic(formal_charges[node])
        a.SetHybridization(node_hybridizations[node])
        a.SetNumExplicitHs(num_explicit_hss[node])
        idx = mol.AddAtom(a)

        # Mapping to keep track of the atom indices graph->atom.
        node_to_idx[node] = idx

    # Get all the edged, i.e. bonds from the networkx graph.
    bond_types = nx.get_edge_attributes(G, "order")

    # Add bonds to molecule.
    for edge in G.edges():

        # Source, target.
        first, second = edge

        # Index parsing from graph->molecule
        ifirst = node_to_idx[first]
        isecond = node_to_idx[second]

        # Get bond type.
        bond_type = bond_types[first, second]

        # Add bond.
        mol.AddBond(ifirst, isecond, bond_type)

    # Explicitly remove aromaticity.
    for atoms in mol.GetAtoms():
        atoms.SetIsAromatic(False)

    # Parse molecule object to smiles.
    smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)

    # Finally return the fragment smiles.
    return smiles


def smiles_to_nx(smiles):
    '''
    Function to parse a smiles to a networkx graph. Node labels and edge labels are parsed as well.
    Labels for nodes are atomic_num, formal_charge, chiral_tag, hybridization, num_explicit_hs, is_aromatic; label for
    edges is order for the type of bond.
    Parameters
    ----------
    smiles string in smiles format that represents a molecule

    Returns
    -------
    networkx graph that is equivalent to the molecule that the supplied smiles encodes.
    '''
    # Create rdkit mol object from smiles
    mol = Chem.MolFromSmiles(smiles)

    # Initialize networkx graph.
    G = nx.Graph()

    # Add node to G for each atom found in molecule. Keep basic atom attributes as node-attributes
    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            formal_charge=atom.GetFormalCharge(),
            chiral_tag=atom.GetChiralTag(),
            hybridization=atom.GetHybridization(),
            num_explicit_hs=atom.GetNumExplicitHs(),
            is_aromatic=atom.GetIsAromatic(),
        )

    # Add all bonds
    for bond in mol.GetBonds():

        # Add the bond to G. Edge attribute order reflects bond type (1, 2, 3, ...)
        G.add_edge(
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), order=bond.GetBondType()
        )
    return G


def mol_to_nx(mol):
    '''
    Parses a rdkit molecule to a networkx graph. Node labels and edge labels are parsed as well.
    Labels for nodes are atomic_num, formal_charge, chiral_tag, hybridization, num_explicit_hs, is_aromatic; label for
    edges is order corresponding to the type of bond.
    Molecule supplied should be kekulized, i.e. should not contain aromatic bonds.
    Parameters
    ----------
    mol rdkit mol object with no aromatic bonds

    Returns networkx graph representing mol
    -------

    '''

    # Initialize networkx graph
    G = nx.Graph()

    # Add all atoms
    for atom in mol.GetAtoms():

        # Add node to G for each atom found in molecule. Keep basic atom attributes as node-attributes
        G.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            formal_charge=atom.GetFormalCharge(),
            chiral_tag=atom.GetChiralTag(),
            hybridization=atom.GetHybridization(),
            num_explicit_hs=atom.GetNumExplicitHs(),
            is_aromatic=atom.GetIsAromatic(),
        )

    # Add all bonds
    for bond in mol.GetBonds():

        # We do not allow aromatic bonds and require kekulization prior to parsing.
        if bond.GetBondType() == 12:
            return None

        # Add the bond to G. Edge attribute order reflects bond type (1, 2, 3, ...)
        G.add_edge(
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), order=bond.GetBondType()
        )
    return G

# This ensures that kekulization finishes in a timely manner. Some molecules lead to very long kekulization
# processes (a lot of kekule forms). Change the timeout parameter to allow for longer time for kekulization.
# Currently, 3 seconds is the cutoff for kekulization of a single molecule.
@pebble_concurrent.process(timeout=3)
def iter_kekul(molset):
    '''
    Function to iterate and de-aromatisizes a list of rdkit mol objects.
    Parameters
    ----------
    molset iterable of rdkit mol objects

    Returns
    -------
    list of non-aromatic forms of supplied molecules. (Requires prior proper kekulization)
    '''
    kekule_set = list()
    for m in molset:
        for a in m.GetAtoms():
            a.SetIsAromatic(False)
        kekule_set.append(m)
    return kekule_set

# Fallback for kekulization. Only one kekule form is returned.
def fallback_kekule(smiles):
    '''
    Returns a single (arbitrary) kekule rdkit mol object parsed from supplied smiles
    Parameters
    ----------
    smiles string that represents a molecule.

    Returns
    -------
    rdkit mol object that represents the supplied smiles.
    '''
    s = Chem.MolFromSmiles(smiles)
    for a in s.GetAtoms():
        a.SetIsAromatic(False)
    return [s]

def rec_cutset(G, S, T, visited, cutsets):
    '''
    Recursive cut iteration. Supplied is a graph G, a current subset (S), nodes to still deal with (T),
    already visited node-subsets (visited) and list where all cuts are collected (cutsets). The idea of the algorithm
    is to check for every possible combination of nodes in G whether there exists a cut that separates S from all
    other nodes in G. Valid in this sense means that removing all edges between S and E(G)\S yields a graph
    that has exactly two connected components.
    Parameters
    ----------
    G networkx graph
    S current subset of nodes (to check if a 2-cut exists for them)
    T nodes that still have to be dealt with (S union T should yield the set of all nodes of G)
    visited list of already visited node sets
    cutsets list where valid cuts are collected

    Returns
    -------
    Nothing. Valid cuts are collected in the supplied cutset.
    '''

    # If the current subset has already been checked we return.
    fS = frozenset(S)
    if fS in visited:
        return

    # Copy of graph G to work non-destructive
    gc = G.copy()

    # Remove all nodes in current subset S from copy of G. If S constitutes a valid fragment, gc
    # will contain exactly one connected component after removal of S.
    gc.remove_nodes_from(S)

    # Check if S is a valid fragment. If not, we return.
    if nx.number_connected_components(gc) > 1:
        return

    # Add S to the already visited list.
    visited.add(fS)

    # Collect all the edges that separate S from the rest of nodes in gc/G.
    cut = [(s, n) for s in S for n in G.neighbors(s) if n not in S]

    # Save the cut in cutset
    cutsets.append(cut)

    # Prepare the next recursive steps of the algorithm. Every neighbor in G of any vertex in S is added into a
    # new set SNew. Then, this algorithm is executed again for the new set SNew which ultimately leads to all
    # combinations of nodes in G being checked (for a valid cut that separated them from all other nodes in G).
    for s in S:
        for n in G.neighbors(s):
            if n in S:
                continue
            SNew = S.union({n})
            TNew = T.difference({n})

            # recursive call
            rec_cutset(G, SNew, TNew, visited, cutsets)


def list_all_cutsets(G):
    '''
    This is in essence a wrapper function for the recursive cutset algorithm implemented in rec_cutset. For explaination
    of the algorithm see the docstring of rec_cutset. This functions initializes S with each individual node in G.
    Calling rec_cutsets with each individual node in G leads to all combinations of nodes in V(G) being checked for
    whether there is a valid cut between these nodes and the rest of nodes in G. Ultimately, list_all_cutsets returns
    a list of all valid cuts in a graph G.
    Parameters
    ----------
    G networkx graph

    Returns
    -------
    list of all cuts in G that yield two connected components.
    '''
    visited = set()
    all_nodes = set(G.nodes)
    cutsets = []
    for n in all_nodes:

        # Initialize S with each individual node.
        S = set({n})
        T = all_nodes.difference(S)

        # Call recursive cutset algorithm
        rec_cutset(G, S, T, visited, cutsets)
    return cutsets

def get_db():
    '''
    Function that returns a predefined mongodb connection.
    Editing of the connection parameters should be done here.
    Returns
    -------
    mongodb connection (pymongo.MongoClient). Different collections can then be accessed via client.collection.
    '''
    # localhost connection string (no user and no password need to be specified)
    client = pymongo.MongoClient("mongodb://localhost:27017")

    # custom connection to mongodb host - replace <USER>,<HOST>,<PASSWORD> with your credentials
    # client = pymongo.MongoClient("mongodb://<USER>:<PASSWORD>@<HOST>:<PORT>/?authMechanism=DEFAULT")

    return client
