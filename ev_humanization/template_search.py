import argparse
import pandas as pd
import numpy as np
from scipy.spatial import distance
from Bio import AlignIO, SeqIO

from ev_humanization.utils.ev_couplings_v4 import EVcouplings
from ev_humanization.utils import utils


def extract_seq_eijs(seq, model):
    aa_map = model.alphabet_map

    seq_eij = np.zeros((len(seq), len(seq)))
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            seq_eij[i, j] = seq_eij[j, i] = abs(
                model.J_ij[i, j, aa_map[seq[i]], aa_map[seq[j]]]
            )

    return seq_eij


def compute_ev_profile(seq, model):
    seq_eij = extract_seq_eijs(seq, model)

    profile = np.zeros(len(seq))
    for i, aa_i in enumerate(seq):
        profile[i] = (
            model.h_i[i, model.alphabet_map[aa_i]] +
            0.5 * np.sum(seq_eij[i])
        )

    return profile


def main():
    args = parse_args()

    # read in files
    ali = AlignIO.read(args.alignment, "fasta")
    model = EVcouplings(args.model, file_format=args.ev_file_format)
    query = SeqIO.read(args.query, "fasta")
    human_seqs = list(SeqIO.parse(args.human_seqs, "fasta"))

    # query sequence
    query_seq = str(query.seq).upper()

    # first sequence in alignment
    top_ali_seq = str(ali[0].seq)

    # map query sequence to model
    mapped_query_seq = utils.align_to_model(query_seq, top_ali_seq)

    # map all human sequences to the model
    mapped_human_seqs = [
        (seq.id, utils.align_to_model(str(seq.seq), top_ali_seq))
        for seq in human_seqs
    ]

    # compute profile of the query sequence
    query_profile = compute_ev_profile(mapped_query_seq, model)

    # compute profile of all human sequences
    human_profiles = [
        (seq_id, compute_ev_profile(mapped_seq, model))
        for seq_id, mapped_seq in mapped_human_seqs
    ]

    # compare query to each human profile
    dists = [
        (seq_id, distance.euclidean(profile, query_profile))
        for seq_id, profile in human_profiles
    ]

    # sort according to profile distances
    dists.sort(key=lambda x: x[1], reverse=False)

    # sequences closest to the query
    if args.top is not None:
        dists = dists[: args.top]

    # write sequence ids and distances to file
    df = pd.DataFrame(dists, columns=["seq_id", "dist"])
    df.to_csv(args.out, index=False)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Human template search by pairwise comparison of evolutionary profiles"
    )

    parser.add_argument("alignment", help="alignment in fasta format that served as input to EVcouplings")
    parser.add_argument("model", help="binary eij couplings file")
    parser.add_argument("query", help="(non-human) query sequence in fasta format")
    parser.add_argument("human_seqs", help="human sequences in fasta format")
    parser.add_argument("out", help="path to out file (*.csv)")

    parser.add_argument("--top", "-t", type=int, help="number of top sequences written to file "
                                                      "(by default, all sequences are written to file)")
    parser.add_argument("--ev_file_format", "-f", choices=["plmc_v1", "plmc_v2"], default="plmc_v2",
                        help="file format of EVcouplings model file (default: 'plmc_v2')")

    return parser.parse_args()


if __name__ == '__main__':
    main()
