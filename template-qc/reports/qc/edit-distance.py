import matplotlib.colors
import matplotlib.pyplot
import pandas as pd
import sys
from Bio import SeqIO
from htmlmin import minify

distfile, fafile, sierrafile, css_file = sys.argv[1:5]
outfiles = sys.argv[5:]

nts = frozenset(("A", "C", "T", "G"))
names = ["Same Subtype / Low Distance", "Same Subtype / High Distance", "Different Subtype / Low Distance"]
assert len(names) == len(outfiles)

dist = pd.read_csv(distfile)

dist["distance"] /= dist.min_length.astype(float)
dist["PatientID1"] = dist.SequenceID1.str.partition("_")[0]
dist["PatientID2"] = dist.SequenceID2.str.partition("_")[0]

align = SeqIO.index(fafile, "fasta")

sierra = pd.read_csv(sierrafile, usecols=["SequenceID", "Subtype"])

css = open(css_file).read()

dist = dist.merge(sierra.rename(columns={"SequenceID": "SequenceID1", "Subtype": "Subtype1"}), how="left", on="SequenceID1")\
           .merge(sierra.rename(columns={"SequenceID": "SequenceID2", "Subtype": "Subtype2"}), how="left", on="SequenceID2")

inter = dist[(dist.PatientID1 != dist.PatientID2)]\
          .sort_values(["distance", "SequenceID1", "SequenceID2"])

same_subtype = inter.Subtype1 == inter.Subtype2

reports = [
    inter.loc[same_subtype,:].sort_values("distance").head(100),
    inter.loc[same_subtype,:].sort_values("distance", ascending=False).head(100),
    inter.loc[~same_subtype,:].sort_values("distance").head(100)
]

def diff_seq(seq1, seq2):
    index = [i for i in range(len(seq1)) if seq1[i] != seq2[i] and seq1[i] != "-" and seq2[i] != "-"]
    return "".join(seq1[i] for i in index), "".join(seq2[i] for i in index)

for name, report, outfile in zip(names, reports, outfiles):
    html = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        "<title>Real-time Phylogeny QC Report: {}</title>".format(name),
        "<style>{} th {{ text-align: left; }}</style>".format(css),
        "</head>",
        "<body>",
        "<h1>Sequences with {}</h1>".format(name),
        '<table class="table">',
        "<thead><th>Priority</th><th>Difference</th><th>ID</th><th>Subtype</th><th>Length</th><th>Gaps</th><th>Ambiguous</th><th>Differing Sites</th></thead>"]

    # HTML data rows
    for i, row in enumerate(report.fillna("-").itertuples()):
        seq1, seq2 = diff_seq(align[row.SequenceID1], align[row.SequenceID2])
        len1 = sum(c != "-" for c in align[row.SequenceID1])
        len2 = sum(c != "-" for c in align[row.SequenceID2])
        gap1 = sum(c == "-" for c in align[row.SequenceID1])
        gap2 = sum(c == "-" for c in align[row.SequenceID2])
        amb1 = sum(c != "-" and c.upper() not in nts for c in align[row.SequenceID1])
        amb2 = sum(c != "-" and c.upper() not in nts for c in align[row.SequenceID2])
        html.append('<tr><th rowspan=2 valign="top">{}</th><th rowspan=2 valign="top">{:0.2f}%</th><td>{}</td><td>{}</td><td>{:,d}</td><td>{:,d}</td><td>{:,d}</td><td><tt>{}</tt></td></tr>'.format(i+1, 100*row.distance, row.SequenceID1, row.Subtype1, len1, gap1, amb1, seq1))
        html.append('<tr><td>{}</td><td>{}</td><td>{:,d}</td><td>{:,d}</td><td>{:,d}</td><td><tt>{}</tt></td></tr>'.format(row.SequenceID2, row.Subtype2, len2, gap2, amb2, seq2))

    html += ["</table>", "</body>", "</html>"]

    with open(outfile, "w") as f:
        f.write(minify("\n".join(html)))

# vim: syntax=python expandtab sw=4 ts=4
