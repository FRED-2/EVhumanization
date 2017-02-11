#!/usr/bin/env python
# encoding: utf-8
"""Module provides the functionality to graft as part of an antibody humanization process.

Residues of the sequence to be humanized are grafted onto a provided human
template sequence.

Three different modes of grafting are supported. (1) CDR grafting, i.e. the
Kabat defined CDRs of the non-human sequence are grafted. (2) SDR Grafting,
i.e. based on a provided three dimensional structure of the non-human antibody
sequence, the residues in contact with the anigen are grafted. (3) Residues
to be grafted can be specified manually.

"""

import sys
import os
from abc import ABCMeta, abstractmethod
import operator
from argparse import ArgumentParser
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import SingleLetterAlphabet
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

from kabat_numbering import KabatNumbering
from contacts_from_structure import download_pdb_file, make_contact_map


class Grafting(object):

    """Abstract class for grafting a subset of residues of the antibody sequence
    to be humanized onto a human template sequence.

    Parameters
    ----------
    sequence : `SeqRecord`
        The non-human antibody sequence to be humanized.
    template : `SeqRecord`
        The human template sequence used to humanize the antibody `sequence`.

    Attributes
    ----------
    sequence : `SeqRecord`
        The non-human antibody sequence to be humanized.
    template : `SeqRecord`
        The human template sequence used to humanize the antibody `sequence`.
    mapping : dict (key=int, value=int)
        Mapping from `sequence` to `template`, as dictionary:
        key=residue number of `sequence`, value=residue number of `template`.
        Residues are not zero-based numbered, i.e. starting from 1.
    residues_to_graft : list of int
        List of residues numbers which correspond to the residues of `sequence`
        to be grafted onto `template`. This variable needs to be set by any
        sub class.
    humanized_seq : `SeqRecord`
        The humanized sequence, i.e. the resulting sequence after grafting
        `residues_to_graft` onto `template`.

    Methods
    -------
    find_residues_to_graft()
        Extract residues to be grafted from the antibody sequence.

    """

    __metaclass__ = ABCMeta

    def __init__(self, sequence, template):
        self.sequence = sequence
        self.template = template
        self.map_sequence_to_template()

    def __str__(self):
        return (('CDRGrafting [sequence=%s, template=%s, mapping=%s, '
                 'residues_to_graft=%s, humanized_seq=%s]')
                % (self.sequence, self.template, self.mapping,
                   self.residues_to_graft, self.humanized_seq))

    @abstractmethod
    def find_residues_to_graft(self):
        """Extract residues to be grafted from the antibody sequence.

        In this method, the field `residues_to_graft` should be filled. Any
        sub class needs to implement this method.

        Raises
        ------
        NotImplementedError
            If the method is not implemented.

        """
        raise NotImplementedError

    def map_sequence_to_template(self, gap_open=-11, gap_extend=-1):
        """Map the antibody sequence to be humanized to the human template sequence.

        Parameters
        ----------
        gap_open : int, optional (default: -11)
            Penalty for introducing a gap during the alignment calculation.
        gap_extend : int, optional (default: -1)
            Penalty for extending a gapped section during the alignment calculation.

        """
        aln_sequence, aln_template, _, _, _ = pairwise2.align.globalds(
            self.sequence.seq.upper(), self.template.seq.upper(),
            matlist.blosum62, gap_open, gap_extend)[0]
        self.mapping = defaultdict(lambda: -1)
        seq_res_num, tem_res_num = 1, 1
        for seq_res, tem_res in zip(aln_sequence, aln_template):
            if seq_res != '-' and tem_res != '-':
                self.mapping[seq_res_num] = tem_res_num
                seq_res_num += 1
                tem_res_num += 1
            elif seq_res == '-':
                tem_res_num += 1
            elif tem_res == '-':
                seq_res_num += 1

    def graft(self):
        """Graft a subset of residues of the non-human sequence onto the human sequence."""
        mod_seq = list(self.template.seq)
        for res_num in self.residues_to_graft:
            tem_res_num = self.mapping[res_num]
            if tem_res_num != -1:
                mod_seq[tem_res_num - 1] = self.sequence.seq[res_num - 1]
        self.humanized_seq = SeqRecord(
            seq=Seq(''.join(mod_seq), SingleLetterAlphabet()),
            id='humanized_%s' % self.sequence.id,
            name='humanized %s' % self.sequence.id,
            description='%s humanized by %s' % (self.sequence.id, self.template.id)
        )
        residues_to_graft_mapped = [self.mapping[res_num]
                                    for res_num in self.residues_to_graft]
        print >> sys.stderr, 'Template indices of residues grafted (CDR residues):'
        print ','.join(map(lambda x: str(x - 1),
                           [x for x in residues_to_graft_mapped if not x == -1]))

    def to_file(self, output):
        """Write the humanized sequence to file in fasta format.

        Parameters
        ----------
        output : str
            Name of the output file.

        """
        with open(output, 'w') as out_file:
            descr = 'Humanization by %s: ' % self.__class__.__name__
            descr += 'fusion of non-human %s and human %s'\
                     % (self.sequence.id, self.template.id)
            out_file.write('>%s %s\n' % (self.humanized_seq.id, descr))
            out_file.write(str(self.humanized_seq.seq))
        print >> sys.stderr, 'Humanized sequence written to %s' % output


class CDRGrafting(Grafting):

    """Humanize an antibody sequence by CDR Grafting.

    Parameters
    ----------
    sequence : `SeqRecord`
        The non-human antibody sequence to be humanized.
    template : `SeqRecord`
        The human template sequence used to humanize the antibody `sequence`.

    Attributes
    ----------
    sequence : `SeqRecord`
        The non-human antibody sequence to be humanized.
    template : `SeqRecord`
        The human template sequence used to humanize the antibody `sequence`.
    mapping : dict (key=int, value=int)
        Mapping from `sequence` to `template`, as dictionary:
        key=residue number of `sequence`, value=residue number of `template`.
        Residues are not zero-based numbered, i.e. starting from 1.
    residues_to_graft : list of int
        List of residues numbers which correspond to the residues of `sequence`
        to be grafted onto `template`.
    humanized_seq : `SeqRecord`
        The humanized sequence, i.e. the resulting sequence after grafting
        `residues_to_graft` onto `template`.
    kabat_seq : `KabatNumbering`
        Kabat numbering of the non-human antibody sequence to be humanized.

    """

    # CDRs of light chain defined by Kabat
    L_1 = ((24, ''), (34, ''))
    L_2 = ((50, ''), (56, ''))
    L_3 = ((89, ''), (97, ''))

    # CDRs of heavy chain defined by Kabat
    H_1 = ((31, ''), (35, 'B'))
    H_2 = ((50, ''), (65, ''))
    H_3 = ((95, ''), (102, ''))

    _ALPHABET = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')

    def __init__(self, sequence, template):
        super(CDRGrafting, self).__init__(sequence, template)
        self.kabat_seq = KabatNumbering(sequence)
        self.find_residues_to_graft()
        self.graft()

    def find_residues_to_graft(self):
        """Extract CDR residues of the antibody sequence to be humanized."""
        mode = self.kabat_seq.kabat_list[0][0][0]
        loop_regions = [CDRGrafting.L_1, CDRGrafting.L_2, CDRGrafting.L_3]\
            if mode == 'L' else [CDRGrafting.H_1, CDRGrafting.H_2, CDRGrafting.H_3]
        self.residues_to_graft = []
        for b, e in loop_regions:
            for i in xrange(b[0], e[0] + 1):
                kabat = (i, '')
                if kabat in self.kabat_seq.kabat_dict:
                    self.residues_to_graft.append(self.kabat_seq.kabat_dict[kabat])
                search_alpha = CDRGrafting._ALPHABET
                if i == e[0]:
                    search_alpha = CDRGrafting._ALPHABET[:operator.indexOf(CDRGrafting._ALPHABET, e[1]) + 1]\
                        if e[1] != '' else []
                for letter in search_alpha:
                    kabat = (i, letter)
                    if kabat in self.kabat_seq.kabat_dict:
                        self.residues_to_graft.append(self.kabat_seq.kabat_dict[kabat])


class SDRGrafting(Grafting):

    """Humanize an antibody sequence by SDR Grafting, i.e. using structural information.

    Parameters
    ----------
    sequence : `SeqRecord`
        The non-human antibody sequence to be humanized.
    template : `SeqRecord`
        The human template sequence used to humanize the antibody `sequence`.
    pdb_id : str
        Unique PDB IDcode.
    chains : list of str
        Chaind IDs of the relevant part of the PDB structure,
        as list: [0] Chain ID of the chain containing the antibody
        sequence of intereset, [1] Chain ID of the chain containing the antigen.
    contact_threshold : float, optional (default: 6)
        Distance threshold defining protein contacts.

    Attributes
    ----------
    sequence : `SeqRecord`
        The non-human antibody sequence to be humanized. If the actual sequence
        to be humanized is only a portion of the sequence of the corresponding
        chain in the PDB file provided, the starting end ending residue numbers
        must be annotated in the header of the fasta file, separated from the
        sequence ID by a forward slash (e.g. >some_id/20-140).
    template : `SeqRecord`
        The human template sequence used to humanize the antibody `sequence`.
    mapping : dict (key=int, value=int)
        Mapping from `sequence` to `template`, as dictionary:
        key=residue number of `sequence`, value=residue number of `template`.
        Residues are not zero-based numbered, i.e. starting from 1.
    residues_to_graft : list of int
        List of residues numbers which correspond to the residues of `sequence`
        to be grafted onto `template`.
    humanized_seq : `SeqRecord`
        The humanized sequence, i.e. the resulting sequence after grafting
        `residues_to_graft` onto `template`.
    pair_dist_list : list of tuple of length 3 (int, int, float)
        List of pairwise distances between antibody and antigen residues.
        Each entry is a tuple: (antibody residue number, antigen residue number,
        distance between the specified residues)

    """

    DEFAULT_CONTACT_THRESHOLD = 6

    def __init__(self, sequence, template, pdb_id, chains,
                 contact_threshold=DEFAULT_CONTACT_THRESHOLD):
        super(SDRGrafting, self).__init__(sequence, template)
        self.set_seq_location()
        self.make_contact_map_wrapper(pdb_id, chains)
        self.find_residues_to_graft(contact_threshold)
        self.graft()

    def set_seq_location(self):
        """Add the location of the non-human antibody sequence provided in the
        fasta header to its record."""
        try:
            bound = self.sequence.id.split('/')[-1]
            loc_feature = SeqFeature(FeatureLocation(int(bound.split('-')[0]),
                                                     int(bound.split('-')[1])),
                                     type='domain')
            self.sequence.features.append(loc_feature)
        except Exception:
            exit('Please provide the location of the non-human antibody domain '
                 'on the sequence of the corresponding PDB chain, '
                 'e.g. >some_id/20-140.')

    def make_contact_map_wrapper(self, pdb_id, chains):
        """Calculate paiwise distances between the antibody and its cognate antigen.

        The method downloads the PDB file, calls a function to calculate all
        pairwise distances of the PDB file ansd then, extracts all relevant
        information.

        Parameters
        ----------
        pdb_id : str
            Unique PDB IDcode.
        chains : list of str
            Chaind IDs of the relevant part of the PDB structure,
            as list: [0] Chain ID of the chain containing the antibody
            sequence of intereset, [1] Chain ID of the chain containing the antigen.

        """
        pdb_file_name = pdb_id + '.pdb'
        download_pdb_file(pdb_id, pdb_file_name)
        pair_dist_lists, sifts_mapping_used = make_contact_map(
            pdb_id, chains, pdb_file_name=pdb_file_name, use_sifts_mapping=True
        )
        os.system('rm ' + pdb_file_name)
        self.pair_dist_list = [
            (up_r1, up_r2, dist)
            for up_r1, _, _, up_r2, _, _, dist in pair_dist_lists[tuple(chains)]
            if up_r1 >= int(self.sequence.features[0].location.start) and
            up_r1 <= int(self.sequence.features[0].location.end)
        ]
        if not sifts_mapping_used:
            if not SDRGrafting.numbering_is_equal(
                self.sequence.features[0].location,
                [r for r, _, _ in self.pair_dist_list]
            ):
                exit('Sifts Mapping not available. SDR Grafting can\'t be applied. '
                     'Try another approach.')

    @staticmethod
    def numbering_is_equal(location, res_list):
        """Test if the provided sequence nunbering and the PDB numbering is equal.

        Parameters
        ----------
        location : `FeatureLocation`
            Location of the non-human antibody sequence on the corresponding
            PDB chain.
        res_list : list of int
            List of residue numbers as given in the PDB file.

        Returns
        -------
        bool
            True, if the numbering is equal. Else, false.

        """
        return range(location.start, location.end + 1) == sorted(list(set(res_list)))

    def find_residues_to_graft(self, contact_threshold=DEFAULT_CONTACT_THRESHOLD):
        """Extract the residues in contact with the antigen.

        Parameters
        ----------
        contact_threshold : float, optional (default: 6)
            Distance threshold defining protein contacts.

        """
        res_list = [r1 for r1, _, dist in self.pair_dist_list
                    if dist <= contact_threshold]
        self.residues_to_graft = sorted(list(set(res_list)))


class ManualGrafting(Grafting):

    """Humanize sequence by grafting of manually specified residues.

    Parameters
    ----------
    sequence : `SeqRecord`
        The non-human antibody sequence to be humanized.
    template : `SeqRecord`
        The human template sequence used to humanize the antibody `sequence`.
    given_residues : list of int
        List of residues numbers which correspond to the residues of `sequence`
        to be grafted onto `template`.

    Attributes
    ----------
    sequence : `SeqRecord`
        The non-human antibody sequence to be humanized. If the actual sequence
        to be humanized is only a portion of the sequence of the corresponding
        chain in the PDB file provided, the starting end ending residue numbers
        must be annotated in the header of the fasta file, separated from the
        sequence ID by a forward slash (e.g. >some_id/20-140).
    template : `SeqRecord`
        The human template sequence used to humanize the antibody `sequence`.
    mapping : dict (key=int, value=int)
        Mapping from `sequence` to `template`, as dictionary:
        key=residue number of `sequence`, value=residue number of `template`.
        Residues are not zero-based numbered, i.e. starting from 1.
    residues_to_graft : list of int
        List of residues numbers which correspond to the residues of `sequence`
        to be grafted onto `template`.

    """

    def __init__(self, sequence, template, given_residues):
        super(ManualGrafting, self).__init__(sequence, template)
        self.find_residues_to_graft(given_residues)
        self.graft()

    def find_residues_to_graft(self, given_residues):
        """Set the residues to graft to the residues provided."""
        self.residues_to_graft = given_residues


ANGSTROM_SYMBOL = r'$\mathrm{\AA}$'  # FIXME : LaTeX not renderd


def init_params(args):
    """Initialize arguments parsed from the command line."""
    required_together = [getattr(args, x) for x in ['pdb_id', 'chains']]
    if not all(required_together) and not required_together == 2 * [None]:
        sys.exit('If structural information is used, both pdb ID (-p) and '
                 'chains (-c) must be specified.')

    args.sequence = SeqIO.read(args.sequence, 'fasta')
    args.template = SeqIO.read(args.template, 'fasta')
    return args


def command_line():
    """Define parser to read arguments from the command line."""
    parser = ArgumentParser(description="""Humanization of a non-human antibody
        sequence, guided by a human template sequence.""")
    parser.add_argument('sequence', help="""The antibody sequence to be humanized,
                        provided in fasta format. If you intend to apply SDR Grafting,
                        please provide the location of the antibody domain
                        on the sequence of the corresponding PDB chain in the
                        header of the fasta entry, separated from the sequence id
                        by a forward slash (e.g. >some_id/20-140).""")
    parser.add_argument('template', help="""The human template sequence used to
                        humanize the non-human antibody sequence, provided in
                        fasta format.""")
    parser.add_argument('out', help="""Name of the output fasta file containing
                        the resulting humanized sequence.""")
    parser.add_argument('--pdb_id', '-p', required=False, help="""PDB ID of the
                        structure to be used to determine antigen contacts of
                        the antibody sequence to be humanized.""")
    parser.add_argument('--chains', '-c', required=False, nargs=2, help="""Chain
                        IDs of the chain containing the antibody sequence (first)
                        and of the chain containing its cognate antigen (second).""")
    parser.add_argument('--contact_threshold', '-ct', required=False, type=float,
                        default=SDRGrafting.DEFAULT_CONTACT_THRESHOLD,
                        help="""Distance threshold defining protein contacts.
                        (default: %.1f%s)""" % (SDRGrafting.DEFAULT_CONTACT_THRESHOLD,
                                                ANGSTROM_SYMBOL))
    parser.add_argument('--manual', '-m', required=False, nargs='+', type=int,
                        help="""The non-human residues to be grafted onto the
                        human template sequence.""")
    return init_params(parser.parse_args())


def main(args):
    """Initialize grafting and write the resulting humanized sequence to file."""
    if args.manual is not None:
        grafting = ManualGrafting(args.sequence, args.template, args.manual)
    elif args.pdb_id is not None:
        grafting = SDRGrafting(args.sequence, args.template, args.pdb_id,
                               args.chains, args.contact_threshold)
    else:
        grafting = CDRGrafting(args.sequence, args.template)

    grafting.to_file(args.out)


if __name__ == '__main__':
    args = command_line()
    main(args)
