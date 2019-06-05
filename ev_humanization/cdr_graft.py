import argparse
from Bio import SeqIO

from ev_humanization.utils import utils


ALPHABET = list(" ABCDEFGHIJKLMNOPQRSTUVWXYZ")
KABAT_LIST = [
    (str(pos) + alpha).strip()
    for pos in range(1, 1000)
    for alpha in ALPHABET
]
KABAT_INDEX = {k: i for i, k in enumerate(KABAT_LIST)}

CDR_LIGHT = [("24", "34"), ("50", "56"), ("89", "97")]
CDR_HEAVY = [("31", "35B"), ("50", "65"), ("95", "102")]

map_to_inds = lambda t: (KABAT_INDEX[t[0]], KABAT_INDEX[t[1]])
CDR_LIGHT_INDS = list(map(map_to_inds, CDR_LIGHT))
CDR_HEAVY_INDS = list(map(map_to_inds, CDR_HEAVY))


def find_cdrs(seq, chain):
    # get Kabat numbering of the sequence
    kabat_numbering = utils.retrieve_kabat_numbering(seq)

    # CDRs differ for light and heavy chain
    chain = chain.lower()
    if chain == "heavy":
        CDR_INDS = CDR_HEAVY_INDS
    elif chain == "light":
        CDR_INDS = CDR_LIGHT_INDS

    # compute CDR positions
    cdrs = []
    for kabat, aa in kabat_numbering:
        if aa != "-":
            kabat_ind = KABAT_INDEX[kabat]
            in_cdr = False
            for cdr_ind, (start_ind, end_ind) in enumerate(CDR_INDS):
                in_cdr = start_ind <= kabat_ind <= end_ind
                if in_cdr:
                    break
            cdrs.append((aa, in_cdr))

    return cdrs


def main():
    args = parse_args()

    # read in files
    query = SeqIO.read(args.query, "fasta")
    template = SeqIO.read(args.template, "fasta")

    query_seq = str(query.seq)
    template_seq = str(template.seq)

    # CDRs of the query sequence
    query_cdrs = find_cdrs(query_seq, args.chain)
    query_cdr_seq = "".join([
        aa if in_cdr and aa != "-" else " "
        for aa, in_cdr in query_cdrs
    ]).split()

    # framework regions of the template
    template_cdrs = find_cdrs(template_seq, args.chain)
    template_fr_seq = "".join([
        aa if not in_cdr and aa != "-" else " "
        for aa, in_cdr in template_cdrs
    ]).split()

    # graft query CDRs onto FR of the template
    graft = (
        template_fr_seq[0].lower() + query_cdr_seq[0].upper() +
        template_fr_seq[1].lower() + query_cdr_seq[1].upper() +
        template_fr_seq[2].lower() + query_cdr_seq[2].upper() +
        template_fr_seq[0].lower()
    )

    # write CDR grafted sequence to file
    with open(args.out, "w") as f:
        f.write(f">{query.id}[CDR]+{template.id}[FR]\n{graft}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Graft CDRs of the query sequence onto the framework regions of the template"
    )

    # required arguments
    parser.add_argument("--query", "-q", required=True, help="query sequence")
    parser.add_argument("--template", "-t", required=True, help="(human) template sequence")
    parser.add_argument("--chain", "-v", required=True, choices=["heavy", "light"],
                        help="heavy/light variable antibody chain")
    parser.add_argument("--out", "-o", required=True, help="path to out file (*.fasta)")

    return parser.parse_args()


if __name__ == '__main__':
    main()
