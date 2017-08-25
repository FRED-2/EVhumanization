from collections import defaultdict
from tempfile import mkstemp
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '../../..'))

import numpy
from os import system

_NUM_DIGITS = 3
_LARGE_DISTANCE = 100000.0

"""Authors: TAH, CPS"""


def calculate_distance_list(pdb_file, chains, sifts_file, model_number=0):
    """
    Calculates a list of minimum residue distances for all residue in one or
    multiple PDB chains

    :param pdf_file: filename of pdb file
    :type pdb_file: str
    :param chains: chain names in the pdb to include in distance calculation
    :type chains: list of str
    :param sifts_file: file name of SIFTS file
    :type sifts_file: str
    :param model_number: number of model to use from pdb structure
    :type model_number: int
    :returns: dictionary with list of residue indeces and pairwise distances
    for each chain combination
    :rtype: dict
    """

    from xml.dom import minidom
    from Bio.PDB.PDBParser import PDBParser
    import Bio.PDB.Polypeptide
    from scipy.spatial.distance import cdist


    pdb_to_up, index_to_res_name = defaultdict(lambda: defaultdict(int)), defaultdict(lambda: defaultdict(int))

    if sifts_file is not None:
        # read SIFTS mapping
        xmldoc = minidom.parse(sifts_file)
        itemlist = xmldoc.getElementsByTagName('residue')

        for s in itemlist:
            cross_refs = s.getElementsByTagName('crossRefDb')
            pdb_res_num, pdb_res_name, pdb_chain, up_res_num, up_res_name = None, None, None, None, None

            for r in cross_refs:
                if r.attributes['dbSource'].value == "PDB":
                    pdb_res_num = r.attributes['dbResNum'].value
                    pdb_res_name = r.attributes['dbResName'].value
                    pdb_chain = r.attributes['dbChainId'].value
                elif r.attributes['dbSource'].value == "UniProt":
                    up_res_num = r.attributes['dbResNum'].value
                    up_res_name = r.attributes['dbResName'].value

            if (pdb_chain in chains) and pdb_res_num != None and up_res_num != None and pdb_res_num.isdigit():
                pdb_to_up[pdb_chain][int(pdb_res_num)] = int(up_res_num)
                index_to_res_name[pdb_chain][int(pdb_res_num)] = up_res_name

    # extract residues from 3D structure
    parser = PDBParser()
    structure = parser.get_structure("somestructure", pdb_file)

    residue_lists = {}
    pair_dist_lists = {}
    atom_coordinates_chains = defaultdict(lambda: defaultdict())

    for chain in chains:
        model = structure[model_number][chain]
        residue_list_unfiltered = [res.get_id() for res in model]

        # filter heteroatoms and other oddities
        residue_list = sorted([resseq for (het, resseq, icode) in residue_list_unfiltered
                               if het == " " and icode == " "])

        residue_lists[chain] = residue_list

        # compile atom coordinates into one NumPy matrix per residue for cdist calculation
        atom_coordinates = defaultdict()
        for res in residue_list:
            atom_coordinates[res] = numpy.zeros((len(model[res]), 3))
            if sifts_file is None:
                print model[res].get_resname()
                index_to_res_name[chain][res] = Bio.PDB.Polypeptide.three_to_one(model[res].get_resname())

            for i, atom in enumerate(model[res]):
                atom_coordinates[res][i] = atom.get_coord()
        atom_coordinates_chains[chain] = atom_coordinates

    # calculate minimum atom distance for each residue pair and store - inter and intra
    # print chains
    # print pdb_to_up
    for chain1 in chains:
        for chain2 in chains:
            pair_dist_list = []
            for r1 in residue_lists[chain1]:
                for r2 in residue_lists[chain2]:
                    if chain1 != chain2 or r1 < r2:
                        dist = numpy.min(cdist(atom_coordinates_chains[chain1][r1], atom_coordinates_chains[chain2][r2], 'euclidean'))
                        # print r1, r1, "-", r2, r2, "-", round(dist, _NUM_DIGITS)
                        if sifts_file is not None and r1 in pdb_to_up[chain1] and r2 in pdb_to_up[chain2]:
                            pair_dist_list.append(
                                (pdb_to_up[chain1][r1], r1, index_to_res_name[chain1][r1],
                                 pdb_to_up[chain2][r2], r2, index_to_res_name[chain2][r2],
                                 round(dist, _NUM_DIGITS))
                                )
                        else:
                            pair_dist_list.append(
                                    (r1, r1, index_to_res_name[chain1][r1],
                                     r2, r2, index_to_res_name[chain2][r2],
                                     round(dist, _NUM_DIGITS))
                                    )
            pair_dist_lists[(chain1, chain2)] = pair_dist_list

    # sifts_mapping_used is true if sifts mapping exists for every chain
    sifts_mapping_used = True
    for chain in chains:
        r = residue_lists[chain][0]
        if sifts_file is None or r1 not in pdb_to_up[chain]:
            sifts_mapping_used = False

    return pair_dist_lists, sifts_mapping_used


def read_contact_map(cm_file_name):
    """ Read residue distances from a file
    :param cm_file_name: File name containing residues
    :returns: dictionary with list of residue indeces and pairwise distances
    for each chain combination
    :rtype: dict
    """
    pair_dist_lists = defaultdict(list)

    with open(cm_file_name, 'r') as f:
        for line in f:
            chain1, chain2, up_res1, pdb_res1, res_name1, up_res2, pdb_res2, res_name2, dist = line.split(' ')
            pair_dist_lists[(chain1, chain2)].append([int(up_res1), int(pdb_res1), res_name1, int(up_res2), int(pdb_res2), res_name2, float(dist)])
    return pair_dist_lists


def make_contact_map(pdb_id, chains, pdb_file_name=None, sifts_file_name=None,
                     use_sifts_mapping=True, model_number=0, verbose=False, out_file=None):
    """
    Calculates a distance map by downloading the ingredients from PDB and SIFTS,
    or using local files if available

    :param pdb_id: PDB database identifier
    :type pdb_id: str
    :param chains: chain names in the pdb to include in distance calculation
    :type chains: list of str
    :param pdf_file_name: filename of local pdb file
    :type pdb_file_name: str
    :param sifts_file_name: file name of local SIFTS file
    :type sifts_file: str
    :paramn model_number: model in PDB structure to use
    :type model_number: int
    :param verbose: show verbose output
    :type verbose: boolean
    :param out_file: filename for contactfile
    :type out_file: str
    :returns: dictionary with list of residue indeces and pairwise distances
    for each chain combination
    :rtype: dict
    """
    if pdb_file_name is None:
        pdb_file_handle, pdb_file_name = mkstemp()
        system("curl -s http://www.rcsb.org/pdb/files/" + pdb_id + ".pdb > " + pdb_file_name)
        if verbose:
            print >> sys.stderr, "PDB file:", pdb_file_name

    if sifts_file_name is None and use_sifts_mapping is True:
        sifts_file_handle, sifts_file_name = mkstemp()
        system("curl -s ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/" + pdb_id.lower() + ".xml.gz | zcat > " + sifts_file_name)
        if verbose:
            print >> sys.stderr, "SIFTS file:", sifts_file_name

    pair_dist_list, sifts_mapping_used = calculate_distance_list(pdb_file_name, chains, sifts_file_name, model_number)
    # write distances to textfile
    if out_file is not None:
        with open(out_file, 'w') as of:
            for chain_pair in pair_dist_list:
                for dist in pair_dist_list[chain_pair]:
                    # print " ".join(chain_pair) + " " + " ".join(map(str, dist))
                    of.write(" ".join(chain_pair) + " " + " ".join(map(str, dist)) + "\n")
        print "Distances written to", out_file
    return pair_dist_list, sifts_mapping_used


def usage():
    print "usage: python {} <pdb> <chains> <use_sifts flag>".format(sys.argv[0])
    exit(1)

def download_pdb_file(pdb_id, pdb_file_name):
    if not os.path.exists(pdb_file_name) or\
            os.path.getsize(pdb_file_name) == 0:
        # print >> sys.stderr, "Downloading", pdb_file_name
        system("curl -sL http://www.rcsb.org/pdb/files/" + pdb_id +
               ".pdb > " + pdb_file_name)

if __name__ == '__main__':
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print usage()

    pdb_id = sys.argv[1]
    chains = []
    if len(sys.argv) > 2:
        # chains = sys.argv[2].split()
        chains = list(sys.argv[2])

    use_sifts_mapping = True
    if len(sys.argv) > 3:
        use_sifts = sys.argv[3]
        if use_sifts == "False" or use_sifts == "false" or "no" in "use_sifts" or "exclude" in "use_sifts":
            use_sifts_mapping = False


    pdb_file_name = pdb_id + ".pdb"
    download_pdb_file(pdb_id, pdb_file_name)

    outfile = pdb_id + "_" + "".join(sorted(chains)) + "_contact_map.txt"
    make_contact_map(pdb_id, chains, pdb_file_name=pdb_file_name,
                     use_sifts_mapping=use_sifts_mapping, out_file=outfile)
