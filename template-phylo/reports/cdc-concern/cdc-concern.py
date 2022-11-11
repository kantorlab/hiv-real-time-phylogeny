import pandas as pd
import sys

cluster_file, patient_file, out_file = sys.argv[1:]

clusters = pd.read_csv(cluster_file, sep="\t", index_col="StudyID")
patients = pd.read_csv(patient_file, index_col="StudyID")

assert clusters.index.is_unique
assert patients.index.is_unique

clusters = clusters.join(patients, how="left")

print(clusters.HIVDxYear.isnull().sum(), "of", len(clusters), "have missing HIV Dx date")

recent = clusters[clusters.HIVDx12mo == 1].ClusterID.value_counts()
print(recent)

concern = recent[recent >= 3].index.unique()

clusters[clusters.ClusterID.isin(concern)].to_csv(out_file)

# vim: syntax=python expandtab sw=4 ts=4
