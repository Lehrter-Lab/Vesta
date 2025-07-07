# CMikolaitis @ USA/DISL

import pandas as pd
import re

## Input params
infile = 'particle.pth'

## Read data

with open(infile) as f:
    lines, data, t = f.readlines(), [], None
    for line in lines:
        if m := re.match(r"\s*(\d+\.\d+)\s+\d+", line): t = float(m[1])
        elif m := re.match(r"\s*(\d+)\s+([\dE+-.]+)\s+([\dE+-.]+)\s+([\dE+-.]+)", line):
            data.append([t, int(m[1]), float(m[2]), float(m[3]), float(m[4])])

df = pd.DataFrame(data, columns=["time", "index", "x", "y", "z"])
print(df.head())