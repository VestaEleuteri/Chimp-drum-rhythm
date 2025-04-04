import thebeat
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np
from theme_drumming import theme_drumming

mpl.rcParams.update(theme_drumming())

df = pd.read_csv(os.path.join("dataframes", "drumming_long.csv"))

western = [757, 935, 980, 990, 605]
eastern = [23, 513, 515, 686]


for bout in western:
    seq = thebeat.Sequence(df[df.Drumming_bout == bout].IBI_ms.values)
    with plt.style.context(theme_drumming()):
        fig, ax = plt.subplots(figsize=(6, 2), tight_layout=True)
        ax.set_prop_cycle(color=["#7E1900"])
    seq.plot_sequence(ax=ax, linewidth=[5] * len(seq.onsets))
    ax.set_xlabel("Time (ms)", fontsize=12)
    ax.set_xlim(0, 2200)
    fig.savefig(os.path.join("plots", "example event plots", "western", f"bout_{bout}.svg"))

for bout in eastern:
    seq = thebeat.Sequence(df[df.Drumming_bout == bout].IBI_ms.values)
    with plt.style.context(theme_drumming()):
        fig, ax = plt.subplots(figsize=(6, 2), tight_layout=True)
        ax.set_prop_cycle(color=["#54B2CF"])
    seq.plot_sequence(ax=ax, linewidth=[5] * len(seq.onsets))
    ax.set_xlabel("Time (ms)", fontsize=12)
    ax.set_xlim(0, 2200)
    fig.savefig(os.path.join("plots", "example event plots", "eastern", f"bout_{bout}.svg"))
