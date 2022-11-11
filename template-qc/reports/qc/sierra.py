import matplotlib.colors
import matplotlib.pyplot
import pandas as pd
import sys
from htmlmin import minify

sierrafile = sys.argv[1]
cssfile    = sys.argv[2]
seq_start  = int(sys.argv[3])
outfiles   = sys.argv[4:]

columns = ["Apobec", "StopCodon", "Unusual"]
assert len(outfiles) == len(columns)

sierra = pd.read_csv(sierrafile,
                     usecols=["SequenceID", "Apobec", "StopCodon", "Unusual", "Subtype",
                              "NRTI", "NNRTI", "PI", "SDRM_PR", "SDRM_RT"])

css = open(cssfile).read()

# Logic for coloring values
cmap = matplotlib.pyplot.get_cmap("autumn")
color_columns = {
    "Apobec": 1.0 / sierra.Apobec.max(),
    "StopCodon": 1.0 / sierra.StopCodon.max(),
    "Unusual": 1.0 / sierra.Unusual.max()
}
def color_value(x, column):
    if column in color_columns:
        scaled = int((cmap.N - 1) * (1.0 - x * color_columns[column]))
        color = matplotlib.colors.rgb2hex(cmap(scaled)[:3])
        return '<td style="text-align: right"><span style="font-weight: bold; color: white; background-color: {}; padding: 5px">{}</span></td>'.format(color, x)
    return "<td>{}</td>".format(x)

# Generate report ordered by each column of interest
for outfile, column in zip(outfiles, columns):

    report = sierra.sort_values([column, "SequenceID"], ascending=False)

    html = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        "<title>Real-time Phylogeny QC Report: Sierra {}</title>".format(column),
        "<style>{} th {{ text-align: left; }}</style>".format(css),
        "</head>",
        "<body>",
        "<h1>Sequences by descending {}</h1>".format(column),
        '<table class="table">',
        "<thead><th>Priority</th>{}</thead>".format("".join("<th>{}</th>".format(c) for c in sierra.columns))]

    # HTML data rows
    i = 1
    for row in report.fillna("-").itertuples():
        if int(row.SequenceID.partition("_")[2]) >= seq_start:
            html.append("<tr><td>{}</td>{}</tr>".format(i, "".join(color_value(getattr(row, c), c) for c in sierra.columns)))
            i += 1

    html += ["</table>", "</body>", "</html>"]

    with open(outfile, "w") as f:
        f.write(minify("\n".join(html)))

# vim: syntax=python expandtab sw=4 ts=4
