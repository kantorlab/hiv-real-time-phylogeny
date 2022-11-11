import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys
from matplotlib.ticker import FuncFormatter

distfile, sierrafile, outfile = sys.argv[1:]

dist = pd.read_csv(distfile)

dist["distance"] /= dist.min_length
dist["PatientID1"] = dist.SequenceID1.str.partition("_")[0]
dist["PatientID2"] = dist.SequenceID2.str.partition("_")[0]

sierra = pd.read_csv(sierrafile, usecols=["SequenceID", "Subtype"])

dist = dist.merge(sierra.rename(columns={"SequenceID": "SequenceID1", "Subtype": "Subtype1"}), how="left", on="SequenceID1")\
           .merge(sierra.rename(columns={"SequenceID": "SequenceID2", "Subtype": "Subtype2"}), how="left", on="SequenceID2")

intra = dist[(dist.PatientID1 == dist.PatientID2) & (dist.Subtype1 == dist.Subtype2)]
subtypes = intra.Subtype1.value_counts()
subtypes = subtypes[subtypes > 50]
intra = intra[intra.Subtype1.isin(subtypes.index)]

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,6))

ax1.set_title("Intra-patient distance by subtype", loc="left")
g = sns.FacetGrid(intra, hue="Subtype1")
g.map(sns.distplot, "distance", hist=False, ax=ax1)
ax1.set_xlabel("% difference")
ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, _: "{:.0%}".format(x)))
ax1.legend(title="Subtype (N>50)")

inter = dist[dist.PatientID1 != dist.PatientID2].drop_duplicates(["PatientID1", "PatientID2"], keep="last")
inter["Same subtype"] = inter.Subtype1 == inter.Subtype2

ax2.set_title("Inter-patient distance by subtype", loc="left")
g = sns.FacetGrid(inter, hue="Same subtype")
g.map(sns.distplot, "distance", hist=False, ax=ax2)
ax2.set_xlabel("% difference")
ax2.xaxis.set_major_formatter(FuncFormatter(lambda x, _: "{:.0%}".format(x)))
ax2.legend(title="Same subtype")

f.tight_layout()
f.savefig(outfile)

# vim: syntax=python expandtab sw=4 ts=4
