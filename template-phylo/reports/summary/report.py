import json
import os
import sys
from datetime import datetime

individuals_file, segmental_file, dataset, out_file = sys.argv[1:]

def diff(x):
    if x > 0:
        return r"\textcolor{{ForestGreen}}{{+{:,d}}}".format(x)
    elif x < 0:
        return r"\textcolor{{Red}}{{{:,d}}}".format(x)
    else: 
        return r"\textcolor{gray}{+0}"

def percdiff(x):
    if x >= 0.1:
        return r"\textcolor{{ForestGreen}}{{+{:.1f}\%}}".format(x)
    elif x <= -0.1:
        return r"\textcolor{{Red}}{{{:.1f}\%}}".format(x)
    else: 
        return r"\textcolor{gray}{+0.0\%}"

stats = json.load(open(stats_file))
clusters = ["{:,d}".format(int(stats["clusters"][0])), diff(int(stats["clusters"][1]))]
cdc_clusters = ["{:,d}".format(int(stats["cdc_clusters"][0])), diff(int(stats["cdc_clusters"][1]))]

dataset_ver, dataset, _ = dataset.split("_")
dataset = "{}-{}-{}".format(dataset[:4], dataset[4:6], dataset[6:])

tex = [r"\documentclass{article}[11pt]",
       r"\usepackage[papersize={14in, 14in},margin=0.5in]{geometry}",
       r"\usepackage{helvet}",
       r"\usepackage{graphicx}",
       r"\usepackage[dvipsnames]{xcolor}",
       r"\usepackage{wasysym}",
       r"\usepackage{booktabs}",
       r"\renewcommand*\familydefault{\sfdefault}",
       r"\renewcommand{\arraystretch}{1.1}",
       r"\begin{document}",
       r"\section*{{Dataset {} ({})}}".format(dataset, dataset_ver),
       r"\begin{tabular}{ll}",
       r"\bf Report generated on: & {}\\".format(datetime.today().strftime("%Y-%m-%d")),
       r"\bf Method for phylogenetic clusters: & Clustered by any of the following:\\",
       r" & RAxML with $\geq 80\%$ bootstrap support and $\leq 4.5 \%$ genetic distance\\",
       r" & IQ-TREE2 (Ultrafast) with $\geq 95\%$ approximate bootstrap support and $\leq 3.0 \%$ genetic distance\\",
       r" & FastTree with $\geq 80\%$ bootstrap support and $\leq 4.5 \%$ genetic distance\\",
       r" & FastTree (aLRT) with $\geq 90\%$ aLRT support and $\leq 3.0 \%$ genetic distance\\",
       r" & MEGA with $\geq 80\%$ bootstrap support and $\leq 4.5 \%$ genetic distance\\",
       r"\bf Method for HIV-TRACE clusters: & HIV-TRACE with $\leq$ 1.5\% genetic distance\\",
       r"\bf Method for CDC clusters of concern: & HIV-TRACE with $\leq$ 0.5\% genetic distance and 3 newly diagnosed cases in 12 months\\",
       r"\bf Definition of index case: & Diagnosed in prior 6 months\\",
       r"\end{tabular}",
       r"\subsection*{{{} Phylogenetic Clusters ({})}}".format(*clusters),
       r"\subsection*{{{} CDC Clusters of Concern ({})}}".format(*cdc_clusters)]

tex.append(open(individuals_file).read())
tex.append(r"\includegraphics[width=13in]{{{}}}".format(os.path.basename(segmental_file)))

# Summary table
rows = []
for row in stats["rows"]:
    if row[0] == "N":
        rows.append(r"\bf {{\em N}} Patients & {:,d} & +{:,d} & {}".format(int(row[1][0]), int(row[1][1]), percdiff(row[1][2])))
        rows += [r" & {:,d} & \em {:.1f}\% & {}".format(int(x[0]), x[1], percdiff(x[2])) for x in row[2:]]
    elif row[0] == "Median Dx Year":
        rows.append(r"\bf {}".format(row[0]))
        for x in row[1:]:
            try:
                cell1 = "{:d}".format(int(x[0]))
            except ValueError:
                cell1 = "-"
            try:
                cell2 = "{:d}".format(int(x[2]))
            except ValueError:
                cell2 = "-"
            rows.append(r" & {} & & \em ({})".format(cell1, cell2))
    else:
        rows.append(r"\bf {}".format(row[0]))
        rows += [" & {:,d} & \em {:.1f}\% & {}".format(int(x[0]), x[1], percdiff(x[2])) for x in row[1:]]
    rows.append(r" \\")

tex.append(r"\end{document}")

with open(out_file, "w") as f:
    f.write("\n".join(tex))

# vim: syntax=python expandtab sw=4 ts=4
