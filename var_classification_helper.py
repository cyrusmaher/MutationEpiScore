import var_ranking_helper as helper
import pandas as pd
import numpy as np


def _validate_traintest_months(df, train_months, test_months):
    # Make sure date format is as expected
    for mm in train_months + test_months:
        assert mm[:4] in ["2019", "2020", "2021"], mm

    # All training months are before the test months, with no overlap
    for tt in train_months:
        for te in test_months:
            assert tt < te

    # Make sure all the requested dates are there'
    obs_months = [str(xx).split("T")[0]
                  for xx in df["monthdate"].dropna().unique()]
    dates_notfound = set(train_months + test_months).difference(obs_months)
    assert (
        len(dates_notfound) == 0
    ), f"Some dates not found in inputdata: {dates_notfound}"

    print(
        "Train and test periods were valid. Train period is prior to test period with no overlap..."
    )


def split_traintest(df, train_months, test_months):
    """
    Split an input dataframe into train and test by month
    Args:
        df: the DataFrame
        train_months: a list of YYYY-mm-dd dates to train on
        test_months: a list of YYYY-mm-dd dates to train on

    train months must occur before test months
    """
    _validate_traintest_months(df, train_months, test_months)
    df_train = df[df["monthdate"].isin(train_months)].copy()
    assert len(df_train) > 0
    df_test = df[df["monthdate"].isin(test_months)].copy()
    assert len(df_test) > 0
    return df_train, df_test

def get_all_obs_varnsites(df):
    all_obs_vars = calculate_features(df, classify=False).query("NCounts > 1").index
    all_obs_sites = sorted(np.unique(helper.var2site(all_obs_vars)).astype(int))
    return all_obs_vars, all_obs_sites


def read_haplotype_summary(fname):
    return pd.read_csv(helper.haplo_file, sep="\t").query(
        "monthdate > '2019'"
    )  # Filter rare invalid dates



def variant_summary_bymonth_and_country(df):
    """
    Calculate summaries per country per month
    """
    all_tmp_c = []
    for (kk, ll), dd in df.groupby(["monthdate", "location"]):
        tmp = pd.concat(helper.calculate_n_haplotypes_wherepresent(dd), axis=1)
        tmp["monthdate"] = kk
        tmp["location"] = ll
        all_tmp_c.append(tmp)
    return pd.concat(all_tmp_c)


def get_total_delta(data):
    """
    Calculate the change from the beginning to end of the period
    """
    months = sorted(data.columns)
    assert str(months[0])[:3] == "202", "Not expected date format"
    return data[months[1]] - data[months[0]]


def get_total_fc(data):
    """
    Calculate the change from the beginning to end of the period
    """
    months = sorted(data.columns)
    assert str(months[0])[:3] == "202", "Not expected date format"
    x1 = data[months[1]]
    x2 = data[months[0]]
    summ = x1 + x2
    return pd.Series(np.where(summ == 0, 0, x1 / summ), index=summ.index)


def get_topn(dd, higher_better=True, topn=3):
    """
    Surface metrics from the topn countries
    topn is set to 3 because most variants are present in 3 or fewer countries
    """
    # Take the top N countries.
    dd = dd.sort_values(ascending=not higher_better).head(topn)
    dd.index = [f"Top{xx+1}" for xx in range(min(len(dd), topn))]
    return dd


def calculate_change_features(df_before, topn_fc=3, topn_delta=2, higher_better=True):
    """
    Extract rate of change featurizations
    """
    month_summary = variant_summary_bymonth_and_country(df_before)

    feature_df_change = []

    for ff in "Frac_HaplosWherePresent", "Frac_Vars":
        # pick a variable and move months to the columns
        data_bycountry = (
            month_summary.set_index(["monthdate", "location"], append=True)[ff]
            .unstack(level="monthdate")
            .fillna(0)
        )

        df_change = (
            get_total_fc(data_bycountry)
            .groupby(level=0)
            .apply(get_topn, topn=topn_fc, higher_better=higher_better)
            .unstack(level=1)
            .fillna(0)
            .add_prefix("FC_")
        )
        df_change2 = (
            get_total_delta(data_bycountry)
            .groupby(level=0)
            .apply(get_topn, topn=topn_delta, higher_better=higher_better)
            .unstack(level=1)
            .fillna(0)
            .add_prefix("Delta_")
        )
        feature_df_change.append((df_change.join(df_change2).add_prefix(ff + "_")))

    feature_df_change = pd.concat(feature_df_change, axis=1)

    return feature_df_change


def calculate_features(df_before, change_features=False, classify=True, **kws):
    """
    Calculate and join cross-sectional and rate-of-change features
    """
    feature_df_cross = pd.concat(
        helper.calculate_n_haplotypes_wherepresent(df_before), axis=1
    )

    if classify:
        epi_cols = ["Frac_HaplosWherePresent", "N_Countries", "Frac_Vars"]
        feature_df_cross = feature_df_cross[epi_cols]
        feature_df_cross["EpiScore"] = (10 ** feature_df_cross.rank(pct=True)).mean(
            axis=1
        )
        feature_df_cross["EpiZScore"] = (
            feature_df_cross[epi_cols]
            .apply(lambda x: (x - x.mean()) / x.std())
            .mean(axis=1)
        )

    if change_features:
        feature_df_change = calculate_change_features(df_before, **kws)
        return feature_df_cross.join(feature_df_change)
    else:
        return feature_df_cross
