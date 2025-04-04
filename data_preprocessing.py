"""
Script for preprocessing raw drumming data (which come from annotations in ELAN).
The script reads in the raw data, cleans it, and calculates various measures (CV, nPVI, entropy, time since recording).
"""

import thebeat
import numpy as np
import pandas as pd
import os
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import datetime


# Melt raw data
df_orig = pd.read_csv(os.path.join("dataframes", "drumming_raw.csv"))
id_vars = [col for col in df_orig.columns if not "Interbeat" in col]
df_molten = pd.melt(df_orig, value_vars=[f"Interbeat {i}" for i in range(1, 27)], id_vars=id_vars)
df_molten.rename({'variable': 'IBI_i', 'value': 'IBI'}, axis=1, inplace=True)
df_molten.IBI_i = df_molten.IBI_i.str.replace("Interbeat ", "").astype(int)
df_molten.Drumming_bout = df_molten.Drumming_bout.astype(int)
#df_molten = df_molten.sort_values(by=['Drumming_bout'])
df_molten = df_molten[~df_molten.IBI.isna()]
df_molten['IBI_ms'] = df_molten.IBI * 1000
df_molten = df_molten.reset_index(drop=True)

# Rename columns and names
df_molten.Subspecies = df_molten.Subspecies.str.replace("EAC", "Eastern")
df_molten.Subspecies = df_molten.Subspecies.str.replace("WAC", "Western")

df_molten = df_molten.rename(columns={"N_beats": "N_hits"})

# Add columns with combined names
df_molten['Pop_com_code'] = df_molten.Population + " (" + df_molten.Community + ")"
df_molten['Com_indiv_code'] = df_molten.Community + " (" + df_molten.Individual + ")"

# Save molten DataFrame
df_molten.to_csv(os.path.join("dataframes", "drumming_long_before_cleaning.csv"), index=False)


# Data cleaning
df = df_molten.copy()
# Use only resting and traveling and remove Community V
df = df[(df.N_hits > 1) & (df.Context.isin(["Traveling", "Resting"])) & (df.Community != "Community V")]
# Remove bouts with fewer than three beats
df = df[df.N_hits > 2]

# Remove individuals with fewer than 9 IBIs (for being able to calculate entropy)
fewer_than_nine_ibis = df.groupby("Com_indiv_code").IBI.nunique() < 9
fewer_than_nine_ibis = fewer_than_nine_ibis[fewer_than_nine_ibis].index
df = df[~df.Com_indiv_code.isin(fewer_than_nine_ibis)]

# Remove individuals with fewer than 3 bouts
fewer_than_three_bouts = df.groupby("Com_indiv_code").Drumming_bout.nunique() < 3
fewer_than_three_bouts = fewer_than_three_bouts[fewer_than_three_bouts].index
df = df[~df.Com_indiv_code.isin(fewer_than_three_bouts)]

# Calculate coefficient of variation and nPVI
for bout, bout_df in df.groupby('Drumming_bout'):
    iois = bout_df['IBI_ms'].values
    seq = thebeat.Sequence(iois)
    npvi = thebeat.stats.get_npvi(seq)
    cov = thebeat.stats.get_cov(seq)
    df.loc[df.Drumming_bout == bout, 'Bout_cv'] = cov
    df.loc[df.Drumming_bout == bout, 'Bout_npvi'] = npvi

# Do k-means clustering and calculate entropy
for indiv, indiv_df in df.groupby("Com_indiv_code"):
    # Get indiidivdual's IBIs, reshape for clustering
    all_ibis = indiv_df.IBI.values.reshape(-1, 1)
    # Keep track of individual's silhouette scores (higher is better)
    sil_scores = {}
    # try k=2 and k=3
    for k in (2, 3):
        clustering = KMeans(n_clusters=k, n_init=100)
        clustering.fit(all_ibis)
        sil_scores[k] = silhouette_score(all_ibis, clustering.labels_)
    # choose the best fit
    best_k = max(sil_scores, key=sil_scores.get)
    clustering = KMeans(n_clusters=best_k, n_init=100).fit(all_ibis)
    if best_k == 2:
        labels = pd.Series(clustering.labels_.astype(str)).replace({"0": "short", "1": "long"})
    else:
        labels = pd.Series(clustering.labels_.astype(str)).replace({"0": "short", "1": "medium", "2": "long"})

    df.loc[df.Com_indiv_code == indiv, "IBI_cluster"] = labels.values

for bout, bout_df in df.groupby("Drumming_bout"):
    clustered_lengths = bout_df.IBI_cluster.to_list()
    pk = [clustered_lengths.count(length) / len(clustered_lengths) for length in clustered_lengths]
    entropy = np.abs(-np.sum(pk * np.log2(pk)))
    df.loc[df.Drumming_bout == bout, "Bout_entropy"] = entropy


# Calculate time since recording
for bout, bout_df in df.groupby("Drumming_bout"):
    date = bout_df.Date.unique()[0]

    if "." in date and len(date) == 10:
        format = "%d.%m.%Y"
    elif "." in date and len(date) == 8:
        format = "%d.%m.%y"
    elif "-" in date:
        format = "%Y-%m-%d"

    if not date == "NV":
        date = datetime.datetime.strptime(date, format)
        Days_before_2024 = (datetime.datetime.strptime('01-01-24', "%d-%m-%y") - date).days
    else:
        Days_before_2024 = np.nan

    df.loc[df.Drumming_bout == bout, "Days_before_2024"] = Days_before_2024
    df.loc[df.Drumming_bout == bout, "Date_year"] = date.year
    df.loc[df.Drumming_bout == bout, "Date_month"] = str(date.month).zfill(2)
    df.loc[df.Drumming_bout == bout, "Date_day"] = str(date.day).zfill(2)


# Export

# Sort DataFrame
df = df.sort_values(by=["Drumming_bout", "IBI_i"]).reset_index(drop=True)

# Save dataframes
df.to_csv(os.path.join("dataframes", "drumming_long.csv"), index=False)

# Export wide-format DataFrame for pDFA
df_unmolten = df.copy()
df_unmolten = df_unmolten.drop(
    columns=[
        "IBI_i",
        "IBI",
        "IBI_ms",
        "IBI_cluster",
    ]
).drop_duplicates()

new_column_order = ['Drumming_bout', 'Subspecies', 'Population', 'Community', 'Individual', 'Pop_com_code', 'Com_indiv_code',
       'Date', 'Date_year', 'Date_month', 'Date_day', 'Days_before_2024', 'With_Ph', 'Start_Ph',
       'End_Ph', 'Start_Ph_N', 'End_Ph_N', 'Context_Before',
       'Context_After', 'Context', 'Total_bout_duration', 'N_hits',
       'Bout_cv', 'Bout_npvi', 'Bout_entropy']

df_unmolten = df_unmolten[new_column_order]

n_interbeat_columns = df.N_hits.max()
for i in range(1, n_interbeat_columns + 1):
    df_unmolten[f"Interbeat {i}"] = np.nan

for bout, bout_df in df.groupby("Drumming_bout"):
    for i in range(1, len(bout_df) + 1):
        df_unmolten.loc[df_unmolten.Drumming_bout == bout, f"Interbeat {i}"] = bout_df.IBI.values[
            i - 1
        ]

df_unmolten.to_csv(os.path.join("dataframes", "drumming_wide.csv"), index=False)


# Make ratios DataFrame
ratios_df = pd.DataFrame(columns=['Population', 'Community', 'Individual', 'Pop_com_code', 'Com_indiv_code', 'Drumming_bout', 'Ratio_i', 'Ratio'])

for bout, bout_df in df.groupby('Drumming_bout'):
    ratios = thebeat.Sequence(bout_df['IBI_ms'].values).interval_ratios_from_dyads
    if len(ratios) < 1:
        continue
    bout_dict = {
        'Subspecies': bout_df.Subspecies.unique()[0],
        'Population': bout_df.Population.unique()[0],
        'Community': bout_df.Community.unique()[0],
        'Individual': bout_df.Individual.unique()[0],
        'Pop_com_code': bout_df.Pop_com_code.unique()[0],
        'Com_indiv_code': bout_df.Com_indiv_code.unique()[0],
        'Drumming_bout': bout,
        'Days_before_2024': bout_df.Days_before_2024.unique()[0],
        'Ratio_i': list(range(1, len(ratios) + 1)),
        'Ratio': ratios
    }

    ratios_df = pd.concat([ratios_df, pd.DataFrame(bout_dict)])

ratios_df.to_csv(os.path.join("dataframes", "drumming_ratios.csv"), index=False)
