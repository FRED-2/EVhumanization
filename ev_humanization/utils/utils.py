import sys
import requests
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo

ABNUM_URL = "http://www.bioinf.org.uk/cgi-bin/abnum/abnum.pl?plain=1&aaseq={}&scheme=-k"


def retrieve_kabat_numbering(seq):
    # get Kabat numbering from AbNum server
    response = requests.get(ABNUM_URL.format(seq))

    if not response.ok:
        sys.exit(
            "The Abnum server used for Kabat numbering could not be reached.\n"
            "You might want to check your internet connection."
        )

    if not response.text.strip():
        sys.exit(
            "Kabat numbering could not be retrieved from the Abnum server.\n"
            "Please consider providing the Kabat nunbering yourself."
        )

    return read_kabat_numbering(response.text.strip())


def read_kabat_numbering(kabat_str):
    numbering = []
    for line in kabat_str.split("\n"):
        split = line.split()
        aa = split[-1]
        kabat_num = "".join(split[: -1])
        numbering.append((kabat_num[1:], aa))

    return numbering


def align_to_model(seq, model_focus_seq):
    # align sequences
    aln_seq, aln_focus_seq, _, _, _ = pairwise2.align.globalds(
        seq.upper(), model_focus_seq.upper(),
        MatrixInfo.blosum62, -11, -1
    )[0]

    # ignore positions where the focus seq is gapped
    mapped_aln_seq = [
        aa_seq
        for aa_seq, aa_focus in zip(aln_seq, aln_focus_seq)
        if aa_focus != "-"
    ]

    # ignore positions that are excluded from the model
    mapped_aln_seq = [
        aa_seq
        for aa_seq, aa_focus in zip(mapped_aln_seq, model_focus_seq)
        if aa_focus.isupper()
    ]

    return "".join(mapped_aln_seq)
