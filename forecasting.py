"""
Usage:
  forecasting.py <infile> <outfolder> [--from_meta|--from_lineage] [--n_days_for_forecast=<n>]

Options:
    --from_meta     Read from metadata input. Will be inferred to be true if input is contains "metadata" but not "lineage"
    --from_lineage  Read from metadata_lineage file. Will be inferred to be true is input contains "lineage"
    --n_days_for_forecast=<n>  Number of days to forecast [default: 90]
"""

import pandas as pd
import numpy as np
from docopt import docopt
import var_classification_helper as varclass
import var_ranking_helper as helper
import utils
import parse_gisaid as gisaid
import os

today = utils.today

def read_input(in_file, from_meta, from_lineage, filter_last_n_days):
    if from_lineage:
        print(f"Reading gisaid lineage summary: {in_file}")
        df = read_lineage_table(in_file, filter_last_n_days=filter_last_n_days)
        print(f"Reading {len(df)} rows from lineage file")
        return df
    elif from_meta:
        print(f"Reading gisaid metadata summary: {in_file}")
        return gisaid.read_gisaid_assummary(fname=in_file, filter_last_n_days=filter_last_n_days)
    else:
        if filter_last_n_days is not None:
            raise ValueError("Cannot filter by granular date if reading from summary file")
        print("Reading directly from summary table")
        return pd.read_csv(in_file, sep="\t")

def format_lineage_table(lineage_table):
    df_tmp = lineage_table[
        ["AA_Substitution", "country", "pango_lineage", "GISAID_clade", "date"]
    ].copy()

    df_tmp["monthdate"] = (
        df_tmp["date"].str.split("-").str[:2].str.join("-")
    )

    df_tmp = df_tmp.rename(
        columns={"AA_Substitution": "haplotype", "country": "location"}
    )

    df_tmp["haplotype"] = (
        df_tmp["haplotype"]
        .str.replace(r"(", "", regex=False)
        .str.replace(r")", "", regex=False)
    )

    collected_counts = df_tmp.groupby(["location", "monthdate"]).size()
    haplotype_counts = (
        df_tmp.groupby(
            ["haplotype", "location", "monthdate", "pango_lineage", "GISAID_clade"]
        )
        .size()
        .rename("haplotype_counts")
    )
    final = (
        haplotype_counts.reset_index()
        .set_index(["location", "monthdate"])
        .join(collected_counts.rename("collected_counts"))
        .reset_index()
    )
    final = final[final["haplotype"].str.len() > 0]
    return final

def read_lineage_table(path, filter_last_n_days=None):
    lineage_table = pd.read_table(path)

    lineage_table = gisaid.filter_by_date(lineage_table, filter_last_n_days)

    return format_lineage_table(lineage_table)


def df2score(df, months, keepall=False):
    """
    Return mutation scores in the specified window
    """
    if months is not None:
        df_mo = df[df["monthdate"].isin(months)]
    else:
        df_mo = df
    feature_df = varclass.calculate_features(df_mo)
    if keepall:
        return feature_df
    return feature_df["EpiScore"]


def df2topscores(df, months):
    return (
        df2score(df, months)
        .loc[lambda x: x > x.quantile(0.95)]
        .sort_values(ascending=False)
    )


def df2pred(df, months):
    """
    Turn mutation scores into a set of predicted mutations (ordered by score)
    """
    return df2topscores(df, months).index





def retrieve_hap_w_vars(mutations, df):
    these_haplos = pd.Series(np.unique(df["haplotype"]))
    matches, haps = zip(
        *sorted(
            zip(
                [
                    len(set(xx.split(", ")).intersection(mutations)) / len(mutations)
                    for xx in these_haplos
                ],
                these_haplos,
            ),
            reverse=True,
        )
    )
    haps = [hh for mm, hh in zip(matches, haps) if mm == np.max(matches)]
    _, best_hap = sorted(
        zip([len(set(hh.split(", ")).difference(mutations)) for hh in haps], haps)
    )[0]
    return best_hap


def write_summaries(in_file, out_folder, from_meta=False, from_lineage=False, n_days_for_forecast=None):
    df_updatepred = read_input(
        in_file, from_meta, from_lineage, filter_last_n_days=n_days_for_forecast)

    haps_bymonth = (
        df_updatepred.groupby("monthdate")["haplotype_counts"].sum().iloc[-8:]
    )

    months_updatepred = haps_bymonth.iloc[-4:].index

    scores_updatepred = df2score(df_updatepred, months_updatepred, keepall=True)

    # Write out predicted mutations with scores
    print("Writing out scores...")
    scores_updatepred.to_csv(f"{out_folder}/scores_{today()}.csv")



def var2site_df(mutations, input_df):
    sites = helper.var2site(mutations)
    new_df = input_df.reindex(sites).copy()
    new_df.index = mutations
    return new_df


if __name__ == "__main__":
    print("Entered main")
    arguments = docopt(__doc__)

    if "lineage" in arguments["<infile>"]:
        arguments["--from_lineage"] = True
    elif "metadata" in arguments["<infile>"]:
        arguments["--from_meta"] = True

    print("Arguments are:")
    print(arguments)
    
    if not os.path.exists(arguments["<outfolder>"]):
        os.mkdir(arguments["<outfolder>"])
    
    write_summaries(
        arguments["<infile>"],
        arguments["<outfolder>"],
        from_meta=arguments["--from_meta"],
        from_lineage=arguments["--from_lineage"],
        n_days_for_forecast=int(arguments["--n_days_for_forecast"]),
    )
