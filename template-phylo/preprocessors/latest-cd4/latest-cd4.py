import pandas as pd
import sys

in_file, out_file = sys.argv[1:]

pd.read_csv(in_file).drop_duplicates("StudyID", keep="last").to_csv(out_file, index=False)

