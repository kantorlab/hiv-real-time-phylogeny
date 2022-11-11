import sys
from Bio import SeqIO

fafile, seq_start, outfile = sys.argv[1:]
seq_start = int(seq_start)

sequences = [seq for seq in SeqIO.parse(fafile, "fasta") if not seq.id.startswith("O")]
print(len(sequences), "sequences")
print("starting at new sequence ID", seq_start)

def length(s):
    return sum(c != "-" for c in s)

def distance(s1, s2):
    return sum(c1 != c2 for (c1, c2) in zip(s1, s2) if c1 != "-" and c2 != "-")

with open(outfile, "w") as f:
    print("SequenceID1,SequenceID2,min_length,distance", file=f)
    for i, seq1 in enumerate(sequences):
        for seq2 in sequences[i:]:
            seq_id = int(seq2.id.partition("_")[2])
            if seq_id >= seq_start and seq_id < 10000:
                print(",".join([seq1.id,
                                seq2.id,
                                str(min(length(seq1.seq), length(seq2.seq))),
                                str(distance(str(seq1.seq).upper(), str(seq2.seq).upper()))]),
                      file=f)

# vim: syntax=python expandtab sw=4 ts=4
