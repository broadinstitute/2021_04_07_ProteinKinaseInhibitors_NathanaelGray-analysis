import os
import glob
import pandas as pd
import numpy as np
import random
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import matplotlib.pyplot as plt


def load_data(exp, plate, filetype):
    """load all data from a single experiment into a single dataframe"""
    path = os.path.join('../profiles',
                        f'{exp}',
                        f'{plate}',
                        f'*_{filetype}')
    files = glob.glob(path)
    df = pd.concat(pd.read_csv(_, low_memory=False) for _ in files)
    return df


def get_metacols(df):
    """return a list of metadata columns"""
    return [c for c in df.columns if c.startswith("Metadata_")]


def get_featurecols(df):
    """returna  list of featuredata columns"""
    return [c for c in df.columns if not c.startswith("Metadata")]


def get_metadata(df):
    """return dataframe of just metadata columns"""
    return df[get_metacols(df)]


def get_featuredata(df):
    """return dataframe of just featuredata columns"""
    return df[get_featurecols(df)]


def remove_negcon_empty_wells(df):
    """return dataframe of non-negative control wells"""
    df = (
        df.query('Metadata_control_type!="negcon"')
        .dropna(subset=['Metadata_broad_sample'])
        .reset_index(drop=True)
    )
    return df


def remove_all_control_empty_wells(df):
    """return dataframe of treatment wells"""
    df = (
        df.query('Metadata_pert_type=="trt"')
            .reset_index(drop=True)
    )
    return df


def remove_empty_wells(df):
    """return dataframe of non-empty wells"""
    df = (
        df.dropna(subset=['Metadata_broad_sample'])
        .reset_index(drop=True)
    )
    return df

def concat_profiles(df1, df2):
    """Concatenate dataframes"""
    if df1.shape[0] == 0:
        df1 = df2.copy()
    else:
        frames = [df1, df2]
        df1 = pd.concat(frames, ignore_index=True, join="inner")

    return df1


def concat_profiles_index(df1, df2):
    """Concatenate dataframes"""
    if df1.shape[0] == 0:
        df1 = df2.copy()
    else:
        frames = [df1, df2]
        df1 = pd.concat(frames, join="inner")

    return df1


def consensus(profiles_df, group_by_feature):
    metadata_df = (
        get_metadata(profiles_df)
            .drop_duplicates(subset=[group_by_feature])
    )

    feature_cols = [group_by_feature] + get_featurecols(profiles_df)
    profiles_df = profiles_df[feature_cols].groupby([group_by_feature]).median().reset_index()

    profiles_df = (
        metadata_df.merge(profiles_df, on='Metadata_broad_sample')
            .drop(columns=['Metadata_Well'])
    )

    return profiles_df


def percent_score(null_dist, corr_dist, how):
    """
    Calculates the Percent strong or percent recall scores
    :param null_dist: Null distribution
    :param corr_dist: Correlation distribution
    :param how: "left", "right" or "both" for using the 5th percentile, 95th percentile or both thresholds
    :return: proportion of correlation distribution beyond the threshold
    """
    if how == 'right':
        perc_95 = np.nanpercentile(null_dist, 95)
        above_threshold = corr_dist > perc_95
        return np.mean(above_threshold.astype(float))*100, perc_95
    if how == 'left':
        perc_5 = np.nanpercentile(null_dist, 5)
        below_threshold = corr_dist < perc_5
        return np.mean(below_threshold.astype(float))*100, perc_5
    if how == 'both':
        perc_95 = np.nanpercentile(null_dist, 95)
        above_threshold = corr_dist > perc_95
        perc_5 = np.nanpercentile(null_dist, 5)
        below_threshold = corr_dist < perc_5
        return (np.mean(above_threshold.astype(float)) + np.mean(below_threshold.astype(float)))*100, perc_95, perc_5


def corr_between_replicates(df, group_by_feature):
    """
        Correlation between replicates
        Parameters:
        -----------
        df: pd.DataFrame
        group_by_feature: Feature name to group the data frame by
        Returns:
        --------
        list-like of correlation values
     """
    replicate_corr = []
    replicate_grouped = df.groupby(group_by_feature)
    for name, group in replicate_grouped:
        group_features = get_featuredata(group)
        corr = np.corrcoef(group_features)
        if len(group_features) == 1:  # If there is only one replicate on a plate
            replicate_corr.append(np.nan)
        else:
            np.fill_diagonal(corr, np.nan)
            replicate_corr.append(np.nanmedian(corr))  # median replicate correlation
    return replicate_corr


def corr_between_replicates_df(df, group_by_feature):
    """
    Correlation between replicates
    :param df: pd.DataFrame
    :param group_by_feature: Feature name to group the data frame by
    :return: Dataframe of correlation values and groups
    """
    replicate_corr_df = pd.DataFrame()
    replicate_grouped = df.groupby(group_by_feature)
    for name, group in replicate_grouped:
        group_features = get_featuredata(group)
        corr = np.corrcoef(group_features)
        if len(group_features) == 1:  # If there is only one replicate on a plate
            replicate_corr = np.nan
        else:
            np.fill_diagonal(corr, np.nan)
            replicate_corr = np.nanmedian(corr)  # median replicate correlation
        replicate_corr_df = replicate_corr_df.append({group_by_feature: name,
                                                      'replicate_correlation': replicate_corr},
                                                      ignore_index=True)
    return replicate_corr_df


def corr_between_non_replicates(df, n_samples, n_replicates, metadata_compound_name):
    """
        Null distribution between random "replicates".
        Parameters:
        ------------
        df: pandas.DataFrame
        n_samples: int
        n_replicates: int
        metadata_compound_name: Compound name feature
        Returns:
        --------
        list-like of correlation values, with a  length of `n_samples`
    """
    df.reset_index(drop=True, inplace=True)
    null_corr = []
    while len(null_corr) < n_samples:
        compounds = random.choices([_ for _ in range(len(df))], k=n_replicates)
        sample = df.loc[compounds].copy()
        if len(sample[metadata_compound_name].unique()) == n_replicates:
            sample_features = get_featuredata(sample)
            corr = np.corrcoef(sample_features)
            np.fill_diagonal(corr, np.nan)
            null_corr.append(np.nanmedian(corr))  # median replicate correlation
    return null_corr


def distribution_plot(df, output_file, metric):

    if metric == 'Percent Replicating':
        metric_col = 'Percent_Replicating'
        null = 'Null_Replicating'
        null_label = 'non-replicates'
        signal = 'Replicating'
        signal_label = 'replicates'
        x_label = 'Replicate correlation'
    elif metric == 'Percent Matching':
        metric_col = 'Percent_Matching'
        null = 'Null_Matching'
        null_label = 'non-matching perturbations'
        signal = 'Matching'
        signal_label = 'matching perturbations'
        x_label = 'Correlation between perturbations targeting the same gene'

    n_experiments = len(df)

    plt.rcParams['figure.facecolor'] = 'white'  # Enabling this makes the figure axes and labels visible in PyCharm Dracula theme
    plt.figure(figsize=[12, n_experiments * 6])

    for i in range(n_experiments):
        plt.subplot(n_experiments, 1, i + 1)
        plt.hist(df.loc[i, f'{null}'], label=f'{null_label}', density=True, bins=20, alpha=0.5)
        plt.hist(df.loc[i, f'{signal}'], label=f'{signal_label}', density=True, bins=20, alpha=0.5)
        plt.axvline(df.loc[i, 'Value_95'], label='95% threshold')
        plt.legend(fontsize=20)
        plt.title(
            f"{df.loc[i, 'Description']}\n" +
            f"{metric} = {df.loc[i, f'{metric_col}']}",
            fontsize=25
        )
        plt.ylabel("density", fontsize=25)
        plt.xlabel(f"{x_label}", fontsize=25)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        sns.despine()
    plt.tight_layout()
    plt.savefig(f'figures/{output_file}')


def draw_plates(df):
    df = (
        df.assign(row=lambda x: x.well_position.str[0:1])
        .assign(col=lambda x: x.well_position.str[1:])
    )
    wells = df[['row', 'col', 'control_type', 'pert_type']].copy()
    wells['col'] = wells.col.astype('int')
    wells['control_type'] = wells.control_type.fillna(wells['pert_type'])
    wells['cat'], uniques = pd.factorize(wells.control_type, sort=True)
    wells_pivot = wells.pivot('row', 'col', 'cat')

    sns.set(rc={'figure.figsize': (24, 16)})
    sns.set(font_scale=3)

    n_cat = wells.control_type.nunique()
    cat = list(uniques)
    color = 'tab10'
    colors = sns.color_palette(color, n_cat+1)[::-1]
    colors.pop(0)

    if n_cat != len(list(wells.control_type.drop_duplicates())):
        n_cat += 1
        cat.insert(0, 'empty')
        colors = sns.color_palette(color, n_cat)[::-1]

    cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

    ax = sns.heatmap(wells_pivot, cmap=cmap, linewidths=.5, linecolor='lightgray', square=True, cbar_kws={'shrink': 1/(2*n_cat), 'aspect': n_cat})

    colorbar = ax.collections[0].colorbar

    r = colorbar.vmax - colorbar.vmin
    colorbar.set_ticks([colorbar.vmin + 0.5 * r / (n_cat) + r * i / (n_cat) for i in range(n_cat)])
    colorbar.set_ticklabels(cat)

    ax.xaxis.tick_top()
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_xlabel('')


def draw_any_plate(df, feature):
    wells = df[['row', 'col', feature]].copy()
    wells['col'] = wells.col.astype('int')
    wells['cat'], uniques = pd.factorize(wells[feature], sort=True)
    wells_pivot = wells.pivot('row', 'col', 'cat')

    sns.set(rc={'figure.figsize': (24, 16)})
    sns.set(font_scale=3)

    n_cat = wells[feature].nunique()
    cat = list(uniques)
    color = 'tab10'
    colors = sns.color_palette(color, n_cat+1)[::-1]
    colors.pop(0)

    if n_cat != len(list(wells[feature].drop_duplicates())):
        n_cat += 1
        cat.insert(0, 'empty')
        colors = sns.color_palette(color, n_cat)[::-1]

    cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

    ax = sns.heatmap(wells_pivot, cmap=cmap, linewidths=.5, linecolor='lightgray', square=True, cbar_kws={'shrink': 1/(2*n_cat), 'aspect': n_cat})

    colorbar = ax.collections[0].colorbar

    r = colorbar.vmax - colorbar.vmin
    colorbar.set_ticks([colorbar.vmin + 0.5 * r / (n_cat) + r * i / (n_cat) for i in range(n_cat)])
    colorbar.set_ticklabels(cat)

    ax.xaxis.tick_top()
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_xlabel('')