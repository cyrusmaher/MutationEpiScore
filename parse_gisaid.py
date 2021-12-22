import pandas as pd
from utils import athome
from collections import defaultdict, Counter
from tqdm import tqdm
import var_ranking_helper as helper

def get_hap_and_muts(row):
    hap = row["AA Substitutions"]
    hap = hap.replace("(", "").replace(")", "")

    muts = helper.split_mutstring(
        hap
    )  # [xx.strip() for xx in hap.split(",") if xx.strip() and "X" not in xx]

    # Make sure you don't double count
    # the same haplotype in a different order
    hap = ",".join(sorted(muts))

    return hap, muts


def count_occurrences(occur_dict):
    return {kk: len(vv) for kk, vv in occur_dict.items()}


def muts_perlineage(x, min_rel_freq=0.1):
    return list((x / x.max()).loc[lambda x: x > min_rel_freq].index)


def _build_episcore_matrix(
    df, mut_haplo_count, mut_countries_count, mut_counts, all_haplos
):
    mut_scores = pd.concat(
        [
            pd.Series(mut_haplo_count, name="FracHaplos") / len(all_haplos),
            pd.Series(mut_countries_count, name="NCountries"),
            pd.Series(mut_counts, name="Prev") / len(df),
        ],
        axis=1,
    )
    mut_scores["EpiScore"] = (10 ** mut_scores.rank(pct=True)).mean(axis=1)

    mut_scores["Gene"] = mut_scores.index.str.split("_").str[0]
    mut_scores["Mut"] = mut_scores.index.str.split("_").str[1]
    return mut_scores


def summarize_mutations(df, months):
    if months is not None:
        assert pd.Series(months).isin(df["year-month"]).all()
        df = df[df["year-month"].isin(months)]
    all_haplos = {}
    mut_counts = Counter()

    mut_regions = defaultdict(dict)
    mut_countries = defaultdict(dict)
    mut_haplos = defaultdict(dict)
    pango_counts = defaultdict(Counter)
    mut_country_counts = defaultdict(int)
    country_counts = defaultdict(int)

    for _, row in tqdm(df.iterrows(), total=len(df)):

        hap, muts = get_hap_and_muts(row)
        # save all unique haplotypes
        all_haplos[hap] = ""
        region, country = row["Location"].split(" / ")[:2]

        pango = row["Pango lineage"]
        country_counts[country] += 1

        for mm in muts:
            mut_country_counts[(mm, country)] += 1
            mut_counts[mm] += 1
            mut_haplos[mm][hap] = ""
            mut_regions[mm][region] = ""
            mut_countries[mm][country] = ""
            pango_counts[pango][mm] += 1
            # mut_states[mm][state] = ""

    mut_haplo_count = count_occurrences(mut_haplos)
    mut_countries_count = count_occurrences(mut_countries)
    mut_county_df = pd.Series(mut_country_counts).unstack(level=1).fillna(0)

    return (
        mut_county_df,
        pd.Series(country_counts),
        _build_episcore_matrix(
            df, mut_haplo_count, mut_countries_count, mut_counts, all_haplos
        ),
    )


def read_gisaid_metadata(fname=athome("Data/SARS2/metadata_oct2021.tsv")):
    print(fname)
    df = pd.read_table(fname).dropna(subset=["AA Substitutions", "Location"])

    # A small number of sequences are short (<5000bp)
    df = df[
        (df["Sequence length"] > 28_000)
        & df["Collection date"].str.contains("-")  # Some samples only have the year
    ]

    df["year-month"] = df["Collection date"].str.split("-").str[:2].str.join("-")
    assert (df["Type"].dropna() == "betacoronavirus").all()
    return df


def gisaid2haplosummary(df, states=False):
    df_tmp = df[
        ["AA Substitutions", "Location", "Pango lineage", "Clade", "year-month"]
    ].copy()

    if states:
        df_tmp = df_tmp[df_tmp["Location"].str.contains("USA")]
        df_tmp["Location"] = (
            df_tmp["Location"].str.split("/").str[2].str.strip().str.title()
        )
        states = df_tmp["Location"].value_counts().index[:51]
        df_tmp = df_tmp[df_tmp["Location"].isin(states)]
    else:
        df_tmp["Location"] = (
            df_tmp["Location"].str.split("/").str[1].str.strip().str.title()
        )
    df_tmp["AA Substitutions"] = (
        df_tmp["AA Substitutions"]
        .str.replace(r"(", "", regex=False)
        .str.replace(r")", "", regex=False)
    )
    df_tmp = df_tmp.rename(
        columns={
            "AA Substitutions": "haplotype",
            "Clade": "GISAID_clade",
            "Pango lineage": "pango_lineage",
            "year-month": "monthdate",
            "Location": "location",
        }
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


def filter_by_date(raw_table, filter_last_n_days):
    if filter_last_n_days is None:
        return raw_table
    else:
        date_col = "Submission date" if "Submission date" in raw_table.columns else "date"
        # Some samples only have the year. Omit these samples
        raw_table = raw_table[raw_table[date_col].str.contains("-")]
        date = pd.to_datetime(raw_table[date_col])
        return raw_table[((date.max() - date).dt.days <= filter_last_n_days)]

def read_gisaid_assummary(
    fname=athome("Data/SARS2/metadata_nov2021.tsv"), states=False,
    filter_last_n_days=None
):
    df = filter_by_date(
        read_gisaid_metadata(fname),
        filter_last_n_days
    )

    return gisaid2haplosummary(df, states=states)

def read_lineage_table(path):
    lineage_table = pd.read_table(path)
    df_tmp = lineage_table[["AA_Substitution", "country", "pango_lineage", "GISAID_clade", "year-month"]].copy()
    df_tmp = df_tmp.rename(
        columns={
            "AA_Substitution": "haplotype",
            "year-month": "monthdate",
            "country": "location",
        }
    )

    df_tmp["haplotype"] = (
        df_tmp["haplotype"].str.replace(r"(", "", regex=False).str.replace(r")", "", regex=False)
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
