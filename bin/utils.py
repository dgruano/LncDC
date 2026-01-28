import re
import sys
import gzip

def load_fasta(filename):
    '''
    load_fasta
    Inputs:
        filename - name of FASTA file to load
    Returns:
        a list of sequences
    '''
    sequences = []
    seq = ''
    
    if filename.endswith((".gz", ".Z", ".z")):
        fd = gzip.open(filename, 'rt')
    else:
        fd = open(filename, 'r')

    for line in fd:
        if line.startswith('>'):
            if seq != '':
                sequences.append(seq.replace('a','A').replace('t','T').replace('g','G').replace('c','C'))
            seq = ''
        else:
            seq += line.strip()
    if seq != '':
        sequences.append(seq.replace('a','A').replace('t','T').replace('g','G').replace('c','C'))
    
    fd.close()
    return sequences


def load_precomputed_ss(ss_file):
    seqs_plus_ss = load_fasta(ss_file)
    pattern = r"([acgtnACGTN]+)([\.()]+)(-?\d+(\.\d+)?)$"
    sequences = []
    secondary_structures = []
    for entry in seqs_plus_ss:
        match = re.match(pattern, entry)
        if match:
            sequences.append(match.group(1).replace('a','A').replace('t','T').replace('g','G').replace('c','C'))
            secondary_structures.append(match.group(2))
        else:
            sys.stderr.write("ERROR: The format of the secondary structure file is incorrect! \n")
            sys.exit(1)
    return sequences, secondary_structures