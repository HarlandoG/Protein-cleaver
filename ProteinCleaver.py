#!/usr/bin/python3

import argparse

parser = argparse.ArgumentParser(description = "Processes FASTA file and places peptide sequences after cleavage by enzyme in output file")
parser.add_argument("-i", "--input", help="input fasta filename")
parser.add_argument("-o", "--output", help="New or existing filename to store peptides")
parser.add_argument("-t", "--trypsin", help="cut fasta sequences using Trypsin and write to output file", action="store_true")
parser.add_argument("-l", "--endolys", help="cut fasta sequences using Endoproteinase Lys-C and write to output file", action="store_true")
parser.add_argument("-a", "--endoarg", help = "cut fasta sequences using Endoproteinase Arg-C and write to output file", action="store_true")
parser.add_argument("-v", "--v8prot", help = "cut fasta sequences using V8 proteinase and write to output file", action="store_true")

args = parser.parse_args()
f = open(args.output, "a+")

def ReadFastaFile(f):
    f = open(args.input, 'r')
    sequences = []
    sequencenames = []
    seq = ""
    for line in f:
        if line.startswith(">"):
            sequencenames.append(line)
            if seq:
                sequences.append(seq)
            seq = ""
        else:
            seq += line.rstrip()
    if seq:
        sequences.append(seq)
    f.close()
    return (sequences, sequencenames)


def enzcut(bases):
    cuts = []
    for sequence in bases:
        seqi = []
        for i, base in enumerate(sequence):
            next = i+1
            if next >= len(sequence) - 1:
                next = i
            if args.trypsin:
                enz = enznames[0]
                if base == 'K'and not sequence[next] == 'P':
                    seqi.append(i)
                if base == 'R' and not sequence[next] == 'P':
                    seqi.append(i)
            if args.endolys:
                enz = enznames[1]
                if base == 'K'and not sequence[next] == 'P':
                    seqi.append(i)
            if args.endoarg:
                enz = enznames[2]
                if base == 'R'and not sequence[next] == 'P':
                    seqi.append(i)
            elif args.v8prot:
                enz = enznames[3]
                if base == 'E'and not sequence[next] == 'P':
                    seqi.append(i)
        cuti = []
        for i, cut in enumerate(seqi):
            if i == 0:
                cuti.append(sequence[:cut+1])
            if i == (len(seqi)-1):
                cuti.append(sequence[cut+1:])
            else:
                cuti.append(sequence[cut+1:seqi[i+1]+1])
        cuts.append(cuti)
    return cuts

sequences, seqnames = ReadFastaFile(args.input)
cuts = enzcut(sequences)


for i, name in enumerate(seqnames):
    for j, cut in enumerate(cuts[i]):
        print("{0} Peptide {1:1d}\n{2}\n".format(''.join(name.split()), enz, j, cut), file=f)
