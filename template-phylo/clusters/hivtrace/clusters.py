import json
import sys
from collections import defaultdict

json_file, out_file = sys.argv[1:]

with open(json_file) as f:
    trace = json.load(f)

# Build edge lookup
edges = defaultdict(set)
for edge in trace["trace_results"]["Edges"]:
    nodes = edge["sequences"]
    assert len(nodes) == 2
    edges[nodes[0]].add(nodes[1])
    edges[nodes[1]].add(nodes[0])

with open(out_file, "w") as f:
    print("StudyID", "SequenceID", "ClusterID", "Edges", sep="\t", file=f)
    for n in trace["trace_results"]["Nodes"]:
        print(n["id"].partition("_")[0], n["id"], "HIVTRACE_{}".format(n["cluster"]), ",".join(sorted(edges[n["id"]])), sep="\t", file=f)

# vim: syntax=python expandtab sw=4 ts=4
