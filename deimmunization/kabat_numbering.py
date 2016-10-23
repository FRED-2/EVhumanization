#!/usr/bin/env python
# encoding: utf-8

import sys
from HTMLParser import HTMLParser
from urllib import urlopen
import argparse
from Bio import SeqIO
from collections import defaultdict

html_data = []

class SimpleTextParser(HTMLParser):
    def handle_data(self, data):
        html_data.append(data)

class KabatNumbering():
    """Call antibody Kabat numbering web server.

    Attributes:
        sequence -- SeqRecord of sequence being numbered
        kabat_list -- list of tuples (kabat identifier, residue)
        kabat_dict -- dictionary: key = kabat identifier as (kabat number, kabat letter),
                                  value = residue number (starting from 1)
    """

    def __init__(self, sequence):
        self.sequence = sequence
        abnum_url = 'http://www.bioinf.org.uk/cgi-bin/abnum/abnum.pl?plain=1&aaseq=%s&scheme=-k' % sequence.seq
        kabat_parser = SimpleTextParser()
        try:
            kabat_parser.feed(urlopen(abnum_url).read())
            kabat_parser.close()
        except IOError:
            print >> stderr, 'connection to %s refused' % abnum_url
            sys.exit(-1)
        self._rawdata = ''.join(html_data)
        self.kabat_list = [(x[:-2],x[-1:]) for x in self._rawdata.split('\n') if not x == '']
        self.kabat_dict = defaultdict(lambda: -1)
        res_num = 1
        for kabat, res in self.kabat_list:
            if not res == '-':
                kabat_key = (int(kabat[1:-1]),kabat[-1]) if kabat[-1].isalpha() else (int(kabat[1:]),'')
                self.kabat_dict[kabat_key] = res_num
                res_num += 1

    def __str__(self):
        return 'KabatNumbering [sequence=%s, kabat_list=%s, kabat_dict=%s]' % (self.sequence, self.kabat_list, self.kabat_dict)

    def to_file(self, output):
        with open(output, 'w') as out_file:
            out_file.write('# Kabat numbering\n')
            out_file.write('# id: %s\n# name: %s\n# description: %s\n'
                % (self.sequence.id,self.sequence.name,self.sequence.description))
            out_file.write(''.join(self._rawdata))
        print 'Kabat file written to %s' % output

def command_line():
    parser = argparse.ArgumentParser(description='Wrapper for Kabat numbering of variable antibody sequence (http://www.bioinf.org.uk/abs/abnum/)')
    parser.add_argument('--sequence', '-s',
        required=True,
        help='The sequence to be numbered (in fasta format)')
    parser.add_argument('--output', '-o',
        required=True,
        help='Output file')
    args = parser.parse_args()

    with open(args.sequence, 'rU') as seq_file:
        seq_record = SeqIO.read(seq_file, 'fasta')
    kabat_numb = KabatNumbering(seq_record)
    kabat_numb.to_file(args.output)

if __name__ == '__main__':
    command_line()
