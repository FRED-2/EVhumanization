import sys
import argparse
import os
import yaml
from collections import OrderedDict
import pandas as pd
import numpy as np
import pickle
from Bio import AlignIO

from ev_humanization.utils import utils
from ev_humanization.utils.ev_couplings_v4 import EVcouplings

ALPHABET_PROTEIN_NOGAP = "ACDEFGHIKLMNPQRSTVWY"
AA_FREQUS_FILE = "../data/aa_freqs_{}_mus.pcl"
MODEL_FILE_IMM = "../data/continuous_immuno_hemiltonian_imm_wconst.mod"
MODEL_FILE_EN = "../data/continuous_immuno_hemiltonian_en_wconst.mod"


def main():
    args = parse_args()

    # make relative paths absolute
    pwd = os.path.abspath(os.path.dirname(__file__))
    freq_file = os.path.join(pwd, AA_FREQUS_FILE.format(args.chain))
    model_file_imm = os.path.join(pwd, MODEL_FILE_IMM)
    model_file_en = os.path.join(pwd, MODEL_FILE_EN)

    # read in config
    with open(args.config) as f:
        config = yaml.safe_load(f)

    if config["sets"]["ignore_pos"] is None:
        config["sets"]["ignore_pos"] = []

    if config["sets"]["exclude_pos"] is None:
        config["sets"]["exclude_pos"] = []

    # read in alignment
    ali = AlignIO.read(args.alignment, "fasta")
    wild_type = str(ali[0].seq)

    # read in EV couplings binary file
    ev_couplings = EVcouplings(
        args.model,
        file_format=args.ev_file_format
    )

    # compute internal offset of the sequence's numbering
    n_lower = 0
    for aa in wild_type:
        if aa.isupper():
            break
        n_lower += 1
    offset = min(ev_couplings.index_list) - n_lower

    # read in alleles
    alleles = pd.read_csv(
        config["sets"]["allele_file"],
        names=["name", "pssm_thresh", "p"]
    )

    # read in PSSMs
    pssms = OrderedDict()
    pssm_thresh_norm = []
    for name, pssm_thresh in zip(alleles.name, alleles.pssm_thresh):
        raw_pssm = getattr(__import__("ev_humanization.utils.tepitopepan_matrices", fromlist=[name]), name)

        # flatten PSSM
        pssm = {}
        for epitope_position, d in raw_pssm.items():
            for aa, value in d.items():
                pssm[(epitope_position, aa)] = value

        for epitope_position in range(config["parameters"]["epi_len"]):
            pssm[(epitope_position, "C")] = 0

        # normalize PSSM
        pssm_mean = np.mean(list(pssm.values()))
        pssm_std = np.std(list(pssm.values()))
        pssm = {
            k: (v - pssm_mean) / pssm_std
            for k, v in pssm.items()
        }

        pssm_thresh_norm.append(
            (pssm_thresh - config["parameters"]["epi_len"] * pssm_mean) / pssm_std
        )

        pssms[name] = pssm

    alleles["pssm_thresh_norm"] = pssm_thresh_norm

    data_file = args.out_basename + ".data"
    with open(data_file, "w") as f:
        # amino acid alphabet, alleles
        f.write(f"set SIGMA := {' '.join(ALPHABET_PROTEIN_NOGAP)};\n")
        f.write(f"set A := {' '.join(alleles.name)};\n\n")

        # wild type length, epitope length, number of mutations
        f.write(f"param N := {len(wild_type)};\n")
        f.write(f"param eN := {config['parameters']['epi_len']};\n")
        f.write(f"param k := {config['parameters']['k']};\n\n")

        f.write("param pssm_thresh :=\n")
        for name, _, _, pssm_thresh_norm in alleles.itertuples(index=False):
            f.write(f"\t{name}\t{pssm_thresh_norm}\n")
        f.write(";\n\n")

        f.write("param p :=\n")
        for name, _, p, _ in alleles.itertuples(index=False):
            f.write(f"\t{name}\t{p:.1f}\n")
        f.write(";\n\n")

        # PSSMs
        f.write("param pssm :=\n")
        for name, pssm in pssms.items():
            f.write(f"[{name},*,*]: {' '.join(map(str, range(1, config['parameters']['epi_len'] + 1)))} :=\n")
            for aa in ALPHABET_PROTEIN_NOGAP:
                f.write(aa + "\t")
                for i in range(config['parameters']['epi_len']):
                    f.write(str(pssm[(i, aa)]) + " ")
                f.write("\n")
        f.write(";\n\n")

        # valid indices of the model
        model_inds = [
            ind for ind in ev_couplings.index_list
            if ind not in config["sets"]["ignore_pos"]
        ]

        # valid pair indices of the model
        eij_inds = sorted(set(
            [(i, j) if i < j else (j, i)
             for i in model_inds
             for j in model_inds
             if i != j and (i < j or j < i)]
        ))

        # pair indices
        eij_inds_str = "\t".join(map(lambda t: f"{t[0] - offset + 1} {t[1] - offset + 1}", eij_inds))
        f.write(f"set Eij := {eij_inds_str};\n")

        # single indices
        hi_inds_str = "\t".join(map(lambda ind: str(ind - offset + 1), model_inds))
        f.write(f"set E := {hi_inds_str};\n\n")

        # wild type sequence
        for i, aa in enumerate(wild_type, start=1):
            f.write(f"set WT[{i}] := {aa.upper()};\n")
        f.write("\n")

        # Kabat numbering of the wild type
        kabat_numbering = utils.retrieve_kabat_numbering(
            wild_type.upper()
        )

        # amino acid frequencies in antibodies
        aa_freqs = pickle.load(open(freq_file, "rb"))

        # compute mutations based on aa frequencies
        mutations = [set(res.upper()) for res in wild_type]
        for i, (kabat, muts) in enumerate(zip(kabat_numbering, mutations)):
            if i + 1 not in config["sets"]["ignore_pos"] + config["sets"]["exclude_pos"]:
                kabat = args.chain[0].upper() + kabat[0]
                try:
                    freqs = aa_freqs[kabat]
                    above_thresh = [aa for aa, _, f in freqs if f >= args.freq_thresh * 100]
                    muts.update(above_thresh)
                except KeyError:
                    print(
                        "Mutations not set for Kabat position", kabat,
                        "since amino acid frequency information not available",
                        file=sys.stderr
                    )

        # allowed mutations per position
        for i, muts in enumerate(mutations, start=1):
            f.write(f"set M[{i}] := {' '.join(muts)};\n")
        f.write("\n")

        # single site EV parameters
        hi = {}
        for i in model_inds:
            for aa in list(ALPHABET_PROTEIN_NOGAP):
                hi[(i, aa)] = -1.0 * ev_couplings.h_i[
                    ev_couplings.index_map[i],
                    ev_couplings.alphabet_map[aa]
                ]
        f.write(f"param h: {' '.join(list(ALPHABET_PROTEIN_NOGAP))} :=\n")
        for i in model_inds:
            f.write(str(i - offset + 1) + "\t" + " ".join(str(hi[i, aa]) for aa in ALPHABET_PROTEIN_NOGAP) + "\n")
        f.write(";\n\n")

        # evolutionary couplings
        eij = {}
        for i, j in eij_inds:
            for ai in mutations[i - offset]:
                for aj in mutations[j - offset]:
                    eij[(i, j, ai, aj)] = -1.0 * ev_couplings.J_ij[
                        ev_couplings.index_map[i],
                        ev_couplings.index_map[j],
                        ev_couplings.alphabet_map[ai],
                        ev_couplings.alphabet_map[aj]
                    ]
        f.write("param eij :=\n")
        for i, j in eij_inds:
            mut_list = list(mutations[j - offset])
            f.write(f"[{i - offset + 1},{j - offset + 1},*,*]: {' '.join(mut_list)} :=\n")
            for ai in mutations[i - offset]:
                f.write(ai + "\t" + " ".join(str(eij[i, j, ai, aj]) for aj in mut_list) + "\n")
        f.write(";\n\n")
        f.write("end;")

    base_path, _ = os.path.splitext(args.out_basename)

    # generate immunogenicity and energy lp files
    for model_file, lp_file in [
        (model_file_imm, base_path + "_imm.lp"),
        (model_file_en, base_path + "_en.lp")
    ]:
        command = f"glpsol -m {model_file} -d {data_file} --check --wcpxlp {lp_file}"
        print(f"RUN COMMAND: {command}")
        os.system(command)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate specification files necessary to run EVdeimmunization"
    )

    # required arguments
    parser.add_argument("--alignment", "-a", required=True,
                        help="alignment in fasta format that served as input to EVcouplings")
    parser.add_argument("--model", "-m", required=True,
                        help="binary eij couplings file with the CDR grafted sequence as target")
    parser.add_argument("--config", "-c", required=True, help="config file in YAML format")
    parser.add_argument("--chain", "-v", required=True, choices=["heavy", "light"],
                        help="heavy/light variable antibody chain")
    parser.add_argument("--out_basename", "-o", required=True,
                        help="basename of out file (multiple files will be created with appropriate extensions added)")

    # optional arguments
    parser.add_argument("--freq_thresh", "-t", type=float, default=0.01,
                        help="amino acid frequency threshold used to determine allowed mutations")
    parser.add_argument("--ev_file_format", "-f", choices=["plmc_v1", "plmc_v2"], default="plmc_v2",
                        help="file format of EVcouplings model file (default: 'plmc_v2')")

    return parser.parse_args()


if __name__ == '__main__':
    main()
