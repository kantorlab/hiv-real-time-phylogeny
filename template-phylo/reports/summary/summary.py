import json
import pandas as pd
import numpy as np
import sys

cluster_file, cdc_file, patient_file, seq_file, sierra_file, out_file = sys.argv[1:]

clusters = pd.read_csv(cluster_file, sep="\t", index_col="StudyID")
cdc      = pd.read_csv(cdc_file, index_col="StudyID")
patients = pd.read_csv(patient_file, index_col="StudyID")
seqs     = pd.read_csv(seq_file, usecols=["StudyID", "SequenceID", "AgeAtSeq", "Year", "TreatmentStatus"], index_col="StudyID")
sierra   = pd.read_csv(sierra_file)

clusters = clusters[["ClusterID"]]
assert clusters.index.is_unique

cdc["CDC"] = 1
cdc = cdc[["CDC"]]
assert cdc.index.is_unique

assert patients.index.is_unique
assert seqs.index.is_unique

sierra["StudyID"] = sierra.SequenceID.str.partition("_")[0]
sierra = sierra[["StudyID", "Subtype", "PI", "NRTI", "NNRTI", "SDRM_PR", "SDRM_RT"]].set_index("StudyID")
print(sierra.Subtype.value_counts())
assert sierra.index.is_unique

patients = seqs.join(patients).join(clusters).join(cdc).join(sierra)

columns = [patients,
           patients[patients.ClusterID.notnull()],
           patients[patients.ClusterID.isnull()],
           patients[patients.CDC.notnull()]]

with open(out_file, "w") as f:
    print("Variable", "All", "Clustered", "Non-Clustered", "CDC Concern", sep=",", file=f)
    for v, X in [
        ["N", [len(c) for c in columns]],
        ["Clusters", [c.ClusterID.nunique() for c in columns]],
        ["Male", [(c.Gender == 1).sum() for c in columns]],
        ["Hispanic", [(c.Ethnicity == 1).sum() for c in columns]],
        ["White", [(c.Race == 1).sum() for c in columns]],
        ["African-American", [(c.Race == 2).sum() for c in columns]],
        ["Asian", [(c.Race == 3).sum() for c in columns]],
        ["MSM", [(c.MSM == 1).sum() for c in columns]],
        ["Prior Incarceration", [(c.EverAtACI == 1).sum() for c in columns]],
        ["Substance Use", [(c.EverSubstanceUse == 1).sum() for c in columns]],
        ["Psychiatric Dx", [(c.EverPsychotic == 1).sum() for c in columns]],
        ["Injected Drugs", [(c.EverIDU == 1).sum() for c in columns]],
        ["Subtype B", [(c.Subtype == "B").sum() for c in columns]],
        ["Treatment Naive", [(c.TreatmentStatus == "N").sum() for c in columns]],
        ["DRMs in Treated", [((c.TreatmentStatus == "T") & (c.PI.notnull() | c.NRTI.notnull() | c.NNRTI.notnull())).sum() for c in columns]],
        ["SDRMs in Naive", [((c.TreatmentStatus == "N") & (c.SDRM_PR.notnull() | c.SDRM_RT.notnull())).sum() for c in columns]],
        ["Median Dx Year", [c.HIVDxYear.median() for c in columns]]]:
        print(v, *X, sep=",", file=f)

# vim: syntax=python expandtab sw=4 ts=4
