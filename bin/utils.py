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
    pattern = r"([acgtunACGTUN]+)([\.()]+)(-?\d+(\.\d+)?)$"
    sequences = []
    secondary_structures = []
    invalid_count = 0
    for idx, entry in enumerate(seqs_plus_ss):
        match = re.match(pattern, entry)
        if match:
            seq = match.group(1).upper().replace('U', 'T')
            ss = match.group(2)
            if len(seq) != len(ss):
                print(f"[load_precomputed_ss] WARNING: Sequence/structure length mismatch at entry {idx}: seq_len={len(seq)}, ss_len={len(ss)}\nSequence: {seq[:50]}...\nStructure: {ss[:50]}...")
                invalid_count += 1
                continue
            sequences.append(seq)
            secondary_structures.append(ss)
        else:
            print(f"[load_precomputed_ss] ERROR: Format error at entry {idx}: {entry}")
            invalid_count += 1
            continue
    print(f"[load_precomputed_ss] Loaded {len(sequences)} valid entries, {invalid_count} invalid entries skipped.")
    if len(sequences) == 0:
        sys.stderr.write("ERROR: No valid sequence/structure pairs loaded from file!\n")
        sys.exit(1)
    return sequences, secondary_structures