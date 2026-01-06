#!/usr/bin/env python3

import pandas as pd
import sys

df_list = []

for fn in sys.argv[1:]:
    df = pd.read_csv(fn)
    df_list.append(df)

df_concat = pd.concat(df_list, ignore_index=True)

df_concat.to_csv(f"rAMELY_summary_statistics.csv", index=False)
