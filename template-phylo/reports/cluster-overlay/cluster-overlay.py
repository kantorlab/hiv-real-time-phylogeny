import csv
import sys

cluster_file1, cluster_file2, out_file = sys.argv[1:]

color1 = "color=red"
color2 = "color=blue"
color_both = "color=purple"
weight_both = "penwidth=5"

nodes = {}
edges = {}

with open(cluster_file1) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        n = row["SequenceID"]
        nodes[n] = {"color": color1}
        for nn in row["Edges"].split(","):
            if nn:
                e = (min(n, nn), max(n, nn))
                edges[e] = {"color": color1, "style": "style=dashed"}

clusters = {}

with open(cluster_file2) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        n = row["SequenceID"]
        if n in nodes:
            nodes[n]["color"] = color_both
            nodes[n]["weight"] = weight_both
        else:
            nodes[n] = {"color": color2}
        cluster = row["ClusterID"]
        if cluster in clusters:
            for nn in clusters[cluster]:
                e = (min(n, nn), max(n, nn))
                if e in edges and edges[e]["color"] == color1:
                    edges[e]["color"] = color_both
                    edges[e]["weight"] = weight_both
                    edges[e]["style"] = "style=solid"
                else:
                    edges[e] = {"color": color2, "style": "style=dashed"}
            #clusters[cluster].append(n)
        else:
            clusters[cluster] = [n]

with open(out_file, "w") as f:
    print("graph g {", file=f)
    for n in nodes:
        nodes[n]["all"] = 'label="",shape=circle'
        print("{} [{}]".format(n, ",".join(sorted(nodes[n].values()))), file=f)
    for e in edges:
        print("{} -- {} [{}]".format(e[0], e[1], ",".join(sorted(edges[e].values()))), file=f)
    print("}", file=f)

# vim: syntax=python expandtab sw=4 ts=4
