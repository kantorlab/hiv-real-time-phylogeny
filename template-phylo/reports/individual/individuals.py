import pandas as pd
import numpy as np
import sys

cluster_file, cdc_file, hivtrace_file, patient_file, seq_file, sierra_file, recency_file, pvl_file, out_file = sys.argv[1:]

clusters = pd.read_csv(cluster_file, sep="\t", index_col="StudyID")
cdc      = pd.read_csv(cdc_file, index_col="StudyID")
hivtrace = pd.read_csv(hivtrace_file, sep="\t", index_col="StudyID")
patients = pd.read_csv(patient_file, index_col="StudyID")
seqs     = pd.read_csv(seq_file, usecols=["StudyID", "SequenceID", "AgeAtSeq", "Year", "TreatmentStatus"], index_col="StudyID")
sierra   = pd.read_csv(sierra_file)
recency  = pd.read_csv(recency_file)
pvl      = pd.read_csv(pvl_file)

clusters = clusters[["ClusterID"]]
assert clusters["ClusterID"].notnull().all()
assert clusters.index.is_unique

cdc["CDCClusterID"] = cdc["ClusterID"]
cdc = cdc[["CDCClusterID"]]
assert cdc.index.is_unique

hivtrace["TRACEClusterID"] = hivtrace["ClusterID"]
hivtrace = hivtrace[["TRACEClusterID"]]
assert hivtrace["TRACEClusterID"].notnull().all()
assert hivtrace.index.is_unique

assert patients.index.is_unique
assert seqs.index.is_unique

sierra["StudyID"] = recency.SequenceID.str.partition("_")[0]
sierra = sierra[["StudyID", "Subtype", "PI", "NRTI", "NNRTI", "SDRM_PR", "SDRM_RT"]].set_index("StudyID")
assert sierra.index.is_unique

recency["Recency"] = recency.Ambiguous / recency.Length
recency["StudyID"] = recency.SequenceID.str.partition("_")[0]
recency = recency[["StudyID", "Recency"]].set_index("StudyID")
assert recency.index.is_unique

pvl = pvl.drop_duplicates("StudyID", keep="last").rename(columns={"Year": "PVLYear"}).set_index("StudyID")

patients.join(seqs).join(clusters).join(cdc).join(hivtrace).join(sierra).join(recency).join(pvl).to_csv(out_file)

# vim: syntax=python expandtab sw=4 ts=4
