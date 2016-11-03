#!/usr/bin/env python
# encoding: utf-8

import sys
from sys import stderr
from os import system
from abc import ABCMeta, abstractmethod
import operator
import argparse
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
    """Abstract class for grafting a subset of residues of the sequence to be
    humanized onto the human template sequence.

    Attributes:
        sequence -- SeqRecord of the sequence to be humanized
        template -- SeqRecord of the human sequence serving as a template
                    in the humanization process
        mapping -- mapping of sequence to template as dictionary:
                   key = residue number of sequence,
                   value = residue number of template
                   (residues not zero-based numbered, i.e. starting from 1)
        residues_to_graft -- list of residues of the sequence to be grafted
                             onto the template
                             (needs to be set in inheriting class)
        humanized_seq -- SeqRecord of the humanized sequence,
                         i.e. resulting sequence by grafting residues_to_graft
                         onto template
    """

    __metaclass__ = ABCMeta

    def __init__(self, sequence, template):
        self.sequence = sequence
        self.template = template
        self.map_sequence_to_template()

    def __str__(self):
        return 'CDRGrafting [sequence=%s, template=%s, mapping=%s, residues_to_graft=%s, humanized_seq=%s]'\
            % (self.sequence, self.template, self.mapping, self.residues_to_graft, self.humanized_seq)

    @abstractmethod
    def find_residues_to_graft():
        """
            Find residues of sequence to graft onto template and fill attribute
            residues_to_graft accordingly.
        """
        pass

    def map_sequence_to_template(self, gap_open=-11, gap_extend=-1):
        """Map sequence to be humanized to human template sequence.

        Args:
            gap_open: penalty for introducing a gap in alignment calculation
                      (default: -11)
            gap_extend: penalty for extending a gapped section in alignment
                        calculation (default: -1)
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
        """
            Graft residues_to_graft onto template and save SeqRecord of the
            resulting humanized sequence.
        """
        mod_seq = list(self.template.seq)
        for res_num in self.residues_to_graft:
            tem_res_num = self.mapping[res_num]
            if tem_res_num != -1:
                mod_seq[tem_res_num-1] = self.sequence.seq[res_num-1]
        self.humanized_seq = SeqRecord(
            seq=Seq(''.join(mod_seq), SingleLetterAlphabet()),
            id='humanized_%s' % self.sequence.id,
            name='humanized %s' % self.sequence.id,
            description='%s humanized by %s' % (self.sequence.id, self.template.id))

    def to_file(self, output):
        with open(output, 'w') as out_file:
            descr = 'Humanization by %s: ' % self.__class__.__name__
            descr += 'fusion of non-human %s and human %s'\
                     % (self.sequence.id, self.template.id)
            out_file.write('>%s %s\n' % (self.humanized_seq.id, descr))
            out_file.write(str(self.humanized_seq.seq))
        print 'Humanized sequence written to %s' % output


class CDRGrafting(Grafting):
    """Humanize sequence by CDR Grafting.

    Additional attributes:
        kabat_seq -- Kabat numbering of the sequence to be humanized
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
        """
            Find residues at the CDRs of the sequence to be humanized and save
            the found residue numbers in residues_to_graft.
        """
        mode = self.kabat_seq.kabat_list[0][0][0]
        loop_regions = [CDRGrafting.L_1, CDRGrafting.L_2, CDRGrafting.L_3]\
            if mode == 'L' else [CDRGrafting.H_1, CDRGrafting.H_2, CDRGrafting.H_3]
        self.residues_to_graft = []
        for b, e in loop_regions:
            for i in xrange(b[0], e[0]+1):
                kabat = (i, '')
                if kabat in self.kabat_seq.kabat_dict:
                    self.residues_to_graft.append(self.kabat_seq.kabat_dict[kabat])
                search_alpha = CDRGrafting._ALPHABET
                if i == e[0]:
                    search_alpha = CDRGrafting._ALPHABET[:operator.indexOf(CDRGrafting._ALPHABET, e[1])+1]\
                        if e[1] != '' else []
                for letter in search_alpha:
                    kabat = (i, letter)
                    if kabat in self.kabat_seq.kabat_dict:
                        self.residues_to_graft.append(self.kabat_seq.kabat_dict[kabat])


class SDRGrafting(Grafting):
    """Humanize sequence by SDR Grafting, i.e. using structural information."""

    def __init__(self, sequence, template, pdb_id, chains, contact_threshold=6):
        super(SDRGrafting, self).__init__(sequence, template)
        self.set_seq_location()
        self.pair_dist_list = self.make_contact_map_wrapper(pdb_id, chains)
        self.find_residues_to_graft(contact_threshold)
        self.graft()

    def set_seq_location(self):
        try:
            bound = self.sequence.id.split('/')[-1]
            loc_feature = SeqFeature(FeatureLocation(int(bound.split('-')[0]),
                                                     int(bound.split('-')[1])),
                                     type='domain')
            self.sequence.features.append(loc_feature)
        except Exception:
            print >> stderr, 'Please provide the location of the domain on the sequence.'
            exit(1)

    def make_contact_map_wrapper(self, pdb_id, chains):
        pdb_file_name = pdb_id + '.pdb'
        download_pdb_file(pdb_id, pdb_file_name)
        pair_dist_lists, sifts_mapping_used = make_contact_map(pdb_id, list(chains),
                                                               pdb_file_name=pdb_file_name,
                                                               use_sifts_mapping=True)
        system('rm ' + pdb_file_name)
        pair_dist_list = [(up_r1, up_r2, dist)
                          for up_r1, _, _, up_r2, _, _, dist in pair_dist_lists[tuple(chains)]
                          if up_r1 >= int(self.sequence.features[0].location.start) and
                          up_r1 <= int(self.sequence.features[0].location.end)]
        if not sifts_mapping_used:
            if not SDRGrafting.numbering_is_equal(self.sequence.features[0].location,
                                                  [r for r, _, _ in pair_dist_list]):
                print >> stderr, 'Sifts Mapping not available.'
                print >> stderr, 'SDR Grafting can\'t be applied. Try another approach.'
                exit(1)
        return pair_dist_list

    @staticmethod
    def numbering_is_equal(location, res_list):
        return range(location.start, location.end+1) == sorted(list(set(res_list)))

    def find_residues_to_graft(self, contact_threshold):
        res_list = [r1 for r1, _, dist in self.pair_dist_list
                    if dist <= contact_threshold]
        self.residues_to_graft = sorted(list(set(res_list)))


class ManualGrafting(Grafting):
    """Humanize sequence by grafting of manually specified non-human residues
    onto the human template sequence.
    """

    def __init__(self, sequence, template, given_residues):
        super(ManualGrafting, self).__init__(sequence, template)
        self.find_residues_to_graft(given_residues)
        self.graft()

    def find_residues_to_graft(self, given_residues):
        self.residues_to_graft = [int(res) for res in given_residues.split(',')]


def command_line():
    parser = argparse.ArgumentParser(description='Initial humanization of an antibody by CDR or SDR Grafting')
    parser.add_argument('--pdb_id', '-p',
                        required=False,
                        help='PDB ID for using structure to determine antigen contacts of the murine antibody domain '
                            +'(default: use CDRs defined by Kabat)')
    parser.add_argument('--chains', '-c',
                        required=False,
                        help='Chains of antiboy domain to be humanized (first) and of its cognate antigen (second), '
                            +'e.g. DB (D: antibody chain, B: antigen chain)')
    parser.add_argument('--contact_threshold', '-x',
                        required=False,
                        help='Distance threshold used to define protein contacts '
                            +'(default: 6)')
    parser.add_argument('--manual', '-m',
                        required=False,
                        help='The non-human residues to be grafted onto the human template sequence, '
                            +'specified by a comma-separated list of residue numbers, e.g. 14,15,16,79,80 '
                            +'(residue numbers starting from 1)')
    parser.add_argument('--sequence', '-s',
                        required=True,
                        help='The sequence to be humanized (in fasta format). '
                            +'If you intend to apply SDR Grafting, please provide the location of the domain '
                            +'on the sequence in the header of the fasta entry, separated from the id '
                            +'by a forward slash (e.g. >some_id/20-140)')
    parser.add_argument('--template', '-t',
                        required=True,
                        help='The human sequence to be used as a template (in fasta format)')
    parser.add_argument('--output', '-o',
                        required=True,
                        help='Output File')
    args = parser.parse_args()

    required_together = [getattr(args, x) for x in ('pdb_id', 'chains')]
    if not all(required_together) and not required_together == 2*[None]:
        print >> stderr, 'If structural information is used, both pdb id and chains must be specified.'
        sys.exit(-1)

    # read in sequences
    with open(args.sequence, 'rU') as seq_file, open(args.template, 'rU') as tem_file:
        seq_record, tem_record = SeqIO.read(seq_file, 'fasta'), SeqIO.read(tem_file, 'fasta')

    if args.manual is not None:
        grafting = ManualGrafting(seq_record, tem_record, args.manual)
    elif args.pdb_id is not None and args.contact_threshold is not None:
        grafting = SDRGrafting(seq_record, tem_record, args.pdb_id, args.chains, contact_threshold)
    elif args.pdb_id is not None and args.contact_threshold is None:
        grafting = SDRGrafting(seq_record, tem_record, args.pdb_id, args.chains)
    else:
        grafting = CDRGrafting(seq_record, tem_record)

    grafting.to_file(args.output)

if __name__ == '__main__':
    command_line()
