import pandas as pd
import sys
from numpy import nan

in_files = sys.argv[1:-1]
out_file = sys.argv[-1]

phylogenies = []

tokens = in_files[0].split(".")
assert tokens[-1] == "tsv" and tokens[-2] == "clusters"
phylogeny = tokens[-3]
phylogenies.append(phylogeny)
clusters = pd.read_csv(in_files[0], sep="\t").rename(columns={"ClusterID": phylogeny})

for in_file in in_files[1:]:
    tokens = in_file.split(".")
    assert tokens[-1] == "tsv" and tokens[-2] == "clusters"
    phylogeny = tokens[-3]
    phylogenies.append(phylogeny)
    clusters = clusters.merge(pd.read_csv(in_file, sep="\t").rename(columns={"ClusterID": phylogeny}),
                              how="outer",
                              on=["StudyID", "SequenceID"])

# Merge clusters
assert clusters["StudyID"].is_unique
clusters = clusters.set_index("StudyID")
clusters["ClusterID"] = nan
tips = set(clusters.index)
cluster_id = 1
while tips:
    index = tips.pop()
    clusters.loc[index, "ClusterID"] = cluster_id
    for phylogeny in phylogenies:
        index_cluster = clusters.loc[index, phylogeny]
        if index_cluster != nan:
            for member in clusters[clusters[phylogeny] == index_cluster].index:
                clusters.loc[member, "ClusterID"] = cluster_id
                if member in tips:
                    tips.remove(member)
    cluster_id += 1

clusters.to_csv(out_file, sep="\t", float_format="%g")

# vim: syntax=python expandtab sw=4 ts=4
