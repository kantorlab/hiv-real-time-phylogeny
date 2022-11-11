import pandas as pd
import sys

log_file, out_file = sys.argv[1:]

with open(out_file, "w") as f:

    print("StudyID", "SequenceID", "ClusterID", sep="\t", file=f)

    log = open(log_file)

    line = next(log)
    while not line.startswith("Found"):
        line = next(log)

    df = pd.read_csv(log, sep="\t", comment="-")
    print(df.head())

    for _, row in df.iterrows():
        tips = row["TipNames"][1:-1].split(", ")
        print(tips)
        for tip in tips:
            study_id = tip.partition("_")[0]
            print(study_id, tip, row["ClusterNumber"], sep="\t", file=f)

# vim: syntax=python expandtab sw=4 ts=4
