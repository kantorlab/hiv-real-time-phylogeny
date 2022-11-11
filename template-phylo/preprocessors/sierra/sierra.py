import json
import pandas as pd
import sys

seq_file, sierra_file, out_file = sys.argv[1:]

sequences = pd.read_csv(seq_file, usecols=["StudyID", "SequenceID"])
ids = frozenset(sequences.SequenceID)
sierra = json.load(open(sierra_file))

columns = ["SequenceID", "Subtype", "Apobec", "StopCodon", "Unusual",
           "PI", "NRTI", "NNRTI", "SDRM_PR", "SDRM_RT"]

drugs = ["ATV/r", "DRV/r", "FPV/r", "IDV/r", "LPV/r", "NFV", "SQV/r", "TPV/r", # PI
         "3TC", "ABC", "AZT", "D4T", "DDI", "FTC", "TDF", "ZDV", # NRTI
         "DOR", "EFV", "ETR", "NVP", "RPV"] # NNRTI

columns += drugs

data = {x: [] for x in columns}

for record in sierra:

    if record["inputSequence"]["header"] not in ids: continue

    data["SequenceID"].append(record["inputSequence"]["header"])

    # Subtype
    data["Subtype"].append(record["subtypeText"].partition(" ")[0])

    # Apobec
    apobec = 0
    for gene in record["alignedGeneSequences"]:
        for mutation in gene["mutations"]:
            apobec += mutation["isApobecDRM"]
    data["Apobec"].append(apobec)

    # Validation
    stopcodon = 0
    unusual = 0
    if "validationResults" in record:
        for validation in record["validationResults"]:
            message = validation["message"]
            if message.startswith("There is 1 stop"):
                stopcodon += 1
            elif message.startswith("There is 1 unusual"):
                unusual += 1
            elif message.startswith("There are "):
                tokens = message.split(" ")
                if tokens[3] == "stop":
                    stopcodon += int(tokens[2])
                elif tokens[3] == "unusual":
                    unusual += int(tokens[2])
    data["StopCodon"].append(stopcodon)
    data["Unusual"].append(unusual)

    # Drug resistance
    mutations = {"PI": set(), "NRTI": set(), "NNRTI": set()}
    sdrms = {"PR": "", "RT": ""}
    scores = {drug: None for drug in drugs}
    for gene in record["drugResistance"]:
        for drug in gene["drugScores"]:
            if drug["drugClass"]["name"] == "INSTI": continue
            scores[drug["drug"]["displayAbbr"]] = drug["score"]
            for partial in drug["partialScores"]:
                for mutation in partial["mutations"]:
                    mutations[drug["drugClass"]["name"]].add(mutation["text"])
    data["PI"].append(",".join(sorted(mutations["PI"])))
    data["NRTI"].append(",".join(sorted(mutations["NRTI"])))
    data["NNRTI"].append(",".join(sorted(mutations["NNRTI"])))
    for gene in record["alignedGeneSequences"]:
        if len(gene["SDRMs"]) > 0:
            sdrms[gene["gene"]["name"]] = ",".join(sdrm["text"] for sdrm in gene["SDRMs"])
    data["SDRM_PR"].append(sdrms["PR"])
    data["SDRM_RT"].append(sdrms["RT"])
    for drug in scores:
        data[drug].append(scores[drug])

df = sequences.merge(pd.DataFrame(data=data, columns=columns), how="inner", on="SequenceID")
assert len(df) == len(sequences)

print(df.Subtype.value_counts())
print(df.Apobec.value_counts())
print(df.StopCodon.value_counts())
print(df.Unusual.value_counts())

df.to_csv(out_file, index=False)

# vim: syntax=python expandtab sw=4 ts=4
