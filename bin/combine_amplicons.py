#!/usr/bin/env python3

import glob
import pandas as pd
from pathlib import Path

tsv_files = glob.glob('./*.amplicon_coverage.bed')
df_list = []
for f in tsv_files:
    f = Path(f)
    name = f.name.split('.amplicon_coverage.bed')[0]
    df = pd.read_csv(f, sep='\t')
    df = df[['amplicon_id', 'fraction_covered']]
    df.rename(columns={'fraction_covered': name}, inplace=True)
    df_list.append(df)

df = pd.DataFrame().join([d.set_index('amplicon_id') for d in df_list], how='outer').dropna().reset_index()
df.to_csv('merged_amplicons.csv', index=False)
