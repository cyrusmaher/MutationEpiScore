import var_classification_helper as varclass
import pandas as pd

def read_test_data():
    df = pd.read_csv("./test_data_haplos.csv")

    train_months = ["2020-10-01", "2020-09-01"]
    test_months = ["2020-11-01", "2020-12-01"]

    df_train, df_test = varclass.split_traintest(df, train_months, test_months)
    return df, df_train, df_test


def _test_haplo(df_train, features_train):
    all_haplos = pd.Series(df_train["haplotype"].unique())
    all_variants = df_train["haplotype"].str.split(", ").explode().unique()
    expected_haplo = pd.Series(
        {vv: all_haplos.str.contains(vv).mean() for vv in all_variants}
    ).loc[features_train.index]

    assert (features_train["Frac_HaplosWherePresent"] == expected_haplo).all()


def _test_countries(df_train, features_train):
    all_variants = df_train["haplotype"].str.split(", ").explode().unique()

    expected_countries = pd.Series(
        {
            vv: len(
                df_train[df_train["haplotype"].str.contains(vv)]["location"].unique()
            )
            for vv in all_variants
        }
    ).loc[features_train.index]

    assert (features_train["N_Countries"] == expected_countries).all()


def _test_counts(df_train, features_train):
    expected_counts = pd.Series(
        {
            ii: df_train[df_train["haplotype"].str.contains(ii)][
                "haplotype_counts"
            ].sum()
            for ii in features_train.index
        }
    ).loc[features_train.index]


    frac_denom = df_train.groupby(["location", "monthdate"])["collected_counts"].first()
    expected_frac = (expected_counts / frac_denom.sum()).loc[features_train.index]

    assert (expected_frac == features_train["Frac_Vars"]).sum()

    expected_relfrac = expected_counts.pipe(lambda x: x / x.sum())




def test_features():
    df, df_train, df_test = read_test_data()
    features_train = varclass.calculate_features(df_train)
    _test_countries(df_train, features_train)
    _test_haplo(df_train, features_train)
    _test_counts(df_train, features_train)


def count_variant(df, variant, countries=["United_Kingdom", "USA"]):
    var_count = (
        df[df["haplotype"].str.contains(variant) & df["location"].isin(countries)]
        .groupby("location")["haplotype_counts"]
        .sum()
    )
    return var_count

if __name__ == "__main__":
    test_features()
    print("All tests pass!")