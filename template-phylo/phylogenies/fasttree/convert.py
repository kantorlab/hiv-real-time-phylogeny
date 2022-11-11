import dendropy
import sys
tree = dendropy.Tree.get(path=sys.argv[1], schema="newick")
tree.write(
    path=sys.argv[2],
    schema="newick",
    suppress_rooting=True,
    edge_label_compose_fn=(lambda e: "{:f}".format(e.length)))
