import pandas as pd
import numpy as np
from collections import Counter
from utils import athome
from collections import defaultdict
import re
from tqdm import tqdm
from copy import deepcopy


var_pat = re.compile("[A-z]([0-9]+)")



def filter2spike(mutations, drop_gene=False):
    if drop_gene:
        return [
            xx.split("_")[1]
            for xx in mutations
            if (xx is not None) and ("Spike_" in xx)
        ]

    else:
        return [xx for xx in mutations if (xx is not None) and ("Spike_" in xx)]


def filter2spike_df(df, drop_gene=False):
    mask = [(xx is not None) and ("Spike_" in xx) for xx in df.index]
    sm_df = df[mask].copy()
    if drop_gene:
        sm_df.index = [xx.replace("Spike_", "") for xx in sm_df.index]

        assert sm_df.index.str.contains("_").sum() == 0
    return sm_df


def get_all_muts(df, min_count=1):
    muts_dict = defaultdict(int)

    for kk, row in tqdm(df.iterrows(), total=len(df)):
        muts = row["haplotype"].split(",")
        for mm in muts:
            muts_dict[mm] += 1
    muts = sorted([kk for kk, vv in muts_dict.items() if vv > min_count])
    print("Found", len(muts), "mutations")
    valid_muts = []
    print("Validating...")
    for mm in tqdm(muts):
        try:
            site_grouper(mm)
            valid_muts.append(mm)
        except Exception as e1:
            pass
    print("Kept", len(valid_muts), "valid mutations")
    return valid_muts


def site_grouper(x):
    if "_" in x:
        try:
            prot, mut = x.split("_")
            match = var_pat.search(mut).group(1)
            assert match
            return prot + "_" + match
        except:
            print(f"'{x}'")
            raise Exception
    else:
        return int(var_pat.search(x).group(1))


def _var2site(vv, pat):
    return int(pat.search(vv).group(1))


def var2site(vec, pat=var_pat, drop_gene=False):
    labels = [site_grouper(vv) for vv in vec]
    if drop_gene:
        return [int(xx.split("_")[1]) for xx in labels if (xx is not None)]
    else:
        return labels


def read_haplotype_summary(fname):
    return pd.read_csv(fname, sep="\t").query(
        "monthdate > '2019'"
    )  # Filter rare invalid dates


def site2tuple(ll):
    return [(int(xx[1:-1]), xx[0], xx[-1]) for xx in ll]


def tuple2site(llist):
    return [ll[1] + str(ll[0]) + ll[2] for ll in llist]


def uniquify(llist):
    new_ll = []

    for ll in llist:
        if ll not in new_ll:
            new_ll.append(ll)
    return new_ll


def nuniq_clades(df):
    return df["GISAID_clade"].dropna().nunique()


def haploswvars(df, variants):
    return df[df["haplotype"].apply(lambda x: np.any([vv in x for vv in variants]))]


def split_mutstring(string):
    return [xx.strip() for xx in string.split(",") if xx.strip() and "X" not in xx]


def count_variants_per_haplotype(df, pangolin_lineage):
    # Subset to the lineage of interest
    df_hap = df[df["pango_lineage"] == pangolin_lineage]

    # Count all the variants in the lineage
    var_counts = defaultdict(int)
    for kk, row in df_hap.iterrows():
        for hh in split_mutstring(row["haplotype"]):
            var_counts[hh] += row["haplotype_counts"]
    return pd.Series(var_counts)


def extend_VOCs(df: pd.DataFrame, VOCs: dict):
    """
    Pull in additional mutations that co-occur with the CDC variants of concern
    Args:
        df: the haplotype dataframe
        VOCs: a dictionary mapping CDC VOC pangolin lineages to a list of variants
    """

    def sort_vars(x):
        try:
            return (x.split("_")[0], site_grouper(x.split("_")[1]))
        except:
            return ("_", "_")

    VOCs = deepcopy(VOCs)
    for vv in tqdm(VOCs):
        var_counts = count_variants_per_haplotype(df, vv)

        # Must occur at least 80% as often as the most common variant
        vars_tokeep = (
            pd.Series(var_counts).pipe(lambda x: x / x.max()).loc[lambda x: x > 0.8]
        )
        toadd = vars_tokeep.index.difference(VOCs[vv])
        VOCs[vv].extend(list(toadd))

        # Sort in positional order
        VOCs[vv] = sorted(  # Sort by gene name, then by integer location
            VOCs[vv], key=sort_vars
        )

    # Remove redundant lineages (e.g. two california lineages w the same mutations)
    all_variants = []
    keys = list(VOCs.keys())
    for vv in keys:
        variants = tuple(VOCs[vv])
        if variants in all_variants:
            del VOCs[vv]
        else:
            all_variants.append(variants)
    return VOCs



def summarize_haplotypes_and_variants(df):
    """
    Print information about the number of distinct variants and haplotypes
    """
    n_haplos = df["haplotype"].nunique()
    print(f"There are {n_haplos} distinct haplotypes in GISAID")

    all_vars = []
    for hh in df["haplotype"]:
        all_vars.extend([hh.strip() for hh in hh.split(",")])
    all_vars = np.unique(all_vars)
    print(f"There are {len(all_vars)} distinct variants")
    return all_vars


def get_variants_persite(df_before):
    site_vars = defaultdict(dict)
    for hh in df_before["haplotype"]:
        vv = hh.split(", ")
        for variant in vv:
            site_vars[site_grouper(variant)][variant] = ""
    return site_vars


def nvariants_persite(df_before):
    site_vars = get_variants_persite(df_before)
    variants_persite = {}

    for kk, vv in site_vars.items():
        for v2, _ in vv.items():
            variants_persite[v2] = len(vv)
    return pd.Series(variants_persite)


def calculate_n_haplotypes_wherepresent(df):
    var_hap_obs = defaultdict(dict)
    country_counter = Counter()
    var_n_obs = Counter()
    collected_counts = {}
    n_haplos = df["haplotype"].nunique()
    for ii in df.index:
        hh = df.loc[ii, "haplotype"]
        cc = df.loc[ii, "location"]
        dd = df.loc[ii, "monthdate"]
        collected_counts[cc, dd] = df.loc[ii, "collected_counts"]
        for variant in hh.split(","):
            variant = variant.strip()

            if (
                variant[-1] == "_"
            ):  # Skip the case where they didn't put in the mutation after the gene name
                continue

            var_n_obs[variant] += df.loc[ii, "haplotype_counts"]
            var_hap_obs[variant][hh] = ""
            country_counter[variant, cc] += df.loc[
                ii, "haplotype_counts"
            ]  # used to be 1

    total_collected = pd.Series(collected_counts).sum()
    percent_haplos_present = (
        pd.Series(
            {kk: len(var_hap_obs[kk]) for kk in var_hap_obs.keys()},
            name="Frac_HaplosWherePresent",
        )
        / n_haplos
    )

    vars_persite = percent_haplos_present.groupby(site_grouper).size().to_dict()
    vars_persite_ser = pd.Series(
        [vars_persite[site_grouper(ii)] for ii in percent_haplos_present.index],
        index=percent_haplos_present.index,
    )
    return (
        percent_haplos_present,
        (pd.Series(country_counter).unstack(level=1) > 1)
        .sum(axis=1)
        .rename("N_Countries"),  # calculate the number of countries where observed,
        (pd.Series(var_n_obs) / total_collected).rename("Frac_Vars"),
        pd.Series(var_n_obs).pipe(lambda x: x / x.sum()).rename("RelFrac_Vars"),
        pd.Series(vars_persite_ser).rename("VarsPerSite"),
        pd.Series(var_n_obs).rename("NCounts"),
    )


def calculate_epi_features(df):
    percent_haplos_present, countries, var_counts = calculate_n_haplotypes_wherepresent(
        df
    )
    return pd.concat(
        [
            percent_haplos_present,
            countries.rename("N_countries"),
            var_counts.rename("Frac_Counts"),
        ],
        axis=1,
    )
