import dendropy
import numpy as np
import os
import pandas as pd
import subprocess
import sys
from datetime import datetime

case_file, tree_file, r_plot_file, new_ids_file, treatment_file, dataset, out_file = sys.argv[1:]

cases = pd.read_csv(case_file, index_col="StudyID")
out_prefix = os.path.splitext(out_file)[0]
new_ids = pd.read_csv(new_ids_file)

NA = "-"
gender = {1: "M", 2: "F", 3: "T", "-": NA, NA: NA}
race = {1: "White", 2: "Black", 3: "Asian", 4: "Other", "White": "White", "Multiple Races": "Mult.", "Hispanic": "Hisp.", "-": NA, NA: NA}
hispanic = {1: r"\CIRCLE", 2: r"\Circle", "Hispanic": r"\CIRCLE", "Non-Hispanic": r"\Circle", "Unknown": NA, "-": NA, NA: NA}
yesno = {1: r"\CIRCLE", 0: r"\Circle", 9: "-", "Y": r"\CIRCLE", "N": r"\Circle", "-": NA, NA: NA}
countries = {
  "United States Of America": "USA",
  "Puerto Rico": "PR",
  "Dominican Republic": "DR",
  "Bahamas, The": "Bahamas"
}
country = lambda x : countries.get(x, x)

def year(x):
    try:
        return int(x)
    except:
        if isinstance(x, float) and x != np.nan:
            return "{:g}".format(x)
        else:
            return NA

def approx(x):
    if isinstance(x, float) and x != np.nan:
        return r"$\approx$ {:,d}".format(int(x))
    else:
        return NA

cases = cases[cases.SequenceID.notnull()]
cases["seq_id"] = cases.SequenceID.str.partition("_")[2].astype(int)
cases["CDCClustered"] = cases.CDCClusterID.notnull().astype(int).apply(yesno.get)
cases["Clustered"] = cases.ClusterID.notnull().astype(int).apply(yesno.get)
cases["TRACEClustered"] = cases.TRACEClusterID.notnull().astype(int).apply(yesno.get)
cases["IndexCase"] = cases.HIVDx6mo.apply(yesno.get)
cases["Age"] = (cases.AgeAtSeq + (datetime.now().year - cases.Year)).fillna(NA).apply(year)
cases["Gender"] = cases.Gender.apply(gender.get)
cases["Race"] = cases.Race.fillna(NA).str.split(",", 1).str[0].apply(race.get)
cases["Hispanic"] = cases.Ethnicity.apply(hispanic.get)
cases["CountryOfBirth"] = cases.CountryOfBirth.fillna(NA).apply(country)
cases["DxYear"] = cases.HIVDxYear.fillna(NA).apply(year)
cases["LastNegativeYear"] = cases.YearOfLastNegativeTest.fillna(NA).apply(year)
cases["SequenceYear"] = cases.Year.fillna(NA).apply(year)
cases["MSM"] = cases.MSM.apply(yesno.get)
cases["Subtype"] = cases.Subtype.fillna(NA).str.replace("_", r"\_")
cases["TreatmentNaive"] = cases.TreatmentStatus.apply(yesno.get)
cases["DRM"] = (cases.PI.notnull() | cases.NRTI.notnull() | cases.NNRTI.notnull()).astype(int).apply(yesno.get)
cases["SDRM"] = (cases.SDRM_PR.notnull() | cases.SDRM_RT.notnull()).astype(int).apply(yesno.get)
cases["ViralLoad"] = cases.PVL.fillna(NA).apply(approx)
cases["ViralLoadYear"] = cases.PVLYear.fillna(NA).apply(year)
cases["EverAtACI"] = cases.EverAtACI.apply(yesno.get)
cases["EverSubstanceUse"] = cases.EverSubstanceUse.apply(yesno.get)
cases["EverPsychotic"] = cases.EverPsychotic.apply(yesno.get)
cases["EverIDU"] = cases.EverIDU.apply(yesno.get)

def percent(cluster, field):
    return "{:.01f}%".format(100.0 * cluster[field] / cluster["SequenceID"])

# Table helper function
def format_table(header, rows, header_offset=0, highlight=None):
    labels = list(header.values())[header_offset:]
    columns = list(header.keys())[header_offset:]
    tex = [r"\begin{tabular}{l%s}" % ("l" * len(labels)),
           " &",
           " & ".join(r"\rotatebox{90}{\bf %s}" % col for col in labels),
           r"\\",
           r"\toprule"]
    for row in rows:
        if highlight and row.Index not in highlight:
            tex.append("%s &" % row.Index)
        else:
            tex.append(r"\bf %s &" % row.Index)
        tex += [" & ".join(str(getattr(row, col)) for col in columns),
                r"\\"]
    tex += [r"\bottomrule",
            r"\end{tabular}"]
    return tex

def format_table_no_index(header, rows, header_offset=0):
    labels = ["Cluster Member"] + list(header.values())[header_offset:]
    columns = ["Cluster"] + list(header.keys())[header_offset:]
    tex = [r"\begin{tabular}{%s}" % ("l" * len(labels)),
           " & ".join(r"\rotatebox{90}{\bf %s}" % col for col in labels),
           r"\\",
           r"\toprule"]
    for row in rows:
        tex += [" & ".join(str(getattr(row, col)) for col in columns),
                r"\\"]
    tex += [r"\bottomrule",
            r"\end{tabular}"]
    return tex

# Summary table
cases["ThisMonth"] = False
cases.loc[new_ids.StudyID, "ThisMonth"] = True
summary = cases.loc[new_ids.StudyID].fillna(NA)
summary = summary.sort_values(["CDCClustered", "Clustered", "IndexCase"])
header = {"CDCClustered": "CDC Cluster of Concern",
          "Clustered": "Phylogenetic Cluster",
          "TRACEClustered": "HIV-TRACE Cluster",
          "IndexCase": "Index Case",
          "Age": "Age",
          "Gender": "Gender",
          "Race": "Race",
          "Hispanic": "Hispanic",
          "CountryOfBirth": "Country of Birth",
          "DxYear": "HIV Diagnosis Year",
          "LastNegativeYear": "Last Negative Test Year",
          "SequenceYear": "Earliest Sequence Year",
          "MSM": "MSM",
          "Subtype": "Subtype",
          "TreatmentNaive": "Treatment Naive",
          "DRM": "DRMs",
          "SDRM": "SDRMs",
          "ViralLoad": "Viral Load",
          "ViralLoadYear": "Viral Load Year",
          "EverAtACI": "Prior Incarceration",
          "EverSubstanceUse": "Substance Use",
          "EverPsychotic": "Mental Health Diagnosis",
          "EverIDU": "Injection Drug Use"
}
tex = [r"\subsection*{{{} New Sequences}}".format(len(summary))]
tex += format_table(header, summary.itertuples())
tex.append(r"\newpage")

# CDC cluster of concern summaries
clusters = set()
for case in summary.itertuples():
    if case.CDCClusterID != NA and case.CDCClusterID not in clusters:
        clusters.add(case.CDCClusterID)
        cluster = cases[cases.CDCClusterID == case.CDCClusterID].fillna(NA)
        tex += [r"\subsection*{{CDC Cluster of Concern {}}}".format(len(clusters)),
                r"\vspace{-0.25in}",
                r"\hspace{2.35in}"]
        tex += format_table(header, cluster.itertuples(), header_offset=3, highlight=summary.index.tolist())

# Phylogenetic cluster summaries

def extract_cluster(tree_file, cluster):
    """
    """
    tree = dendropy.Tree.get(file=open(tree_file), schema="newick")
    cluster = [x.partition("_")[0] for x in cluster]
    taxa = []
    for i in cluster:
        for taxon in tree.taxon_namespace:
            if taxon.label.startswith(i):
                taxa.append(taxon)
    root = tree.mrca(taxa=taxa)
    assert root is not None
    root = root.parent_node
    return [n.taxon.label.partition(" ")[0] for n in root.leaf_nodes()]

clusters = set()
for case in summary.itertuples():
    if case.ClusterID != NA and case.ClusterID not in clusters:
        clusters.add(case.ClusterID)
        cluster_file = "{}_{}.csv".format(out_prefix, case.Index)
        cluster_ids = cases.loc[cases.ClusterID == case.ClusterID, "SequenceID"]
        new_cluster_ids = extract_cluster(tree_file, cluster_ids.tolist())
        cluster = cases[cases.index.isin(new_cluster_ids)].fillna(NA)
        cluster["Cluster"] = (cluster["ClusterID"] == case.ClusterID)
        cluster.to_csv(cluster_file, columns=["SequenceID", "ThisMonth"])
        tips_file = "{}_{}.tips.csv".format(out_prefix, case.Index)
        subprocess.run(["Rscript", r_plot_file, case.Index, tree_file, cluster_file, treatment_file, dataset, "{}_{}.pdf".format(out_prefix, case.Index), tips_file], check=True)
        tips = pd.read_csv(tips_file, names=["SequenceID"], skiprows=1).reset_index()
        cluster = cluster.reset_index().merge(tips, on="SequenceID").set_index("StudyID")
        cluster["Cluster"] = cluster["Cluster"].apply({True: r"\CIRCLE", False: r"\Circle"}.get)
        tex += [r"\vspace{0.25in}",
                r"\subsection*{{Phylogenetic Cluster {}}}".format(len(clusters)),
                r"\vspace{-0.5in}",
                r"\begin{minipage}{\textwidth}",
                r"\hspace{-0.3in}",
                r"\begin{minipage}{4.75in}",
                r"\vspace{1.83in}",
                r"\includegraphics[width=6in]{{../individual/{}_{}}}".format(os.path.basename(out_prefix), case.Index),
                r"\end{minipage}",
                r"\begin{minipage}{9.19in}"]
        tex += format_table_no_index(header, cluster.sort_values("index", ascending=False).itertuples(), header_offset=3)
        tex += [r"\end{minipage}",
                r"\end{minipage}"]

with open(out_file, "w") as f:
    f.write("\n".join(tex))

# vim: syntax=python expandtab sw=4 ts=4
