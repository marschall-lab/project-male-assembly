#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import collections as col

import pandas as pd
import numpy as np
import sklearn.preprocessing as prep
import sklearn.linear_model as model


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        '--data-table',
        '-d',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='data_table',
    )
    parser.add_argument(
        '--regression-target',
        '-t',
        type=str,
        choices=['ctgnum', 'assmlen', 'ctgng50'],
        dest='target'
    )
    parser.add_argument(
        '--out-data',
        '-od',
        type=lambda x: pl.Path(x).resolve(),
        dest='out_data',
    )
    parser.add_argument(
        '--out-model',
        '-om',
        type=lambda x: pl.Path(x).resolve(),
        dest='out_model'
    )
    parser.add_argument(
        '--num-cpu',
        '-n',
        type=int,
        default=1,
        dest='num_cpu'
    )
    parser.add_argument(
        '--cv-folds',
        '-cv',
        type=int,
        default=5,
        dest='cv_folds'
    )
    args = parser.parse_args()
    return args


def load_data_table(table, regress_target):

    df = pd.read_csv(table, sep='\t', header=0)
    # drop QV, because two samples are missing
    # a QV estimate (no short reads)
    df.drop('qv', axis=1, inplace=True)
    df.set_index('sample', drop=True, inplace=True)

    # drop haplogroup and project info; e.g., haplogroup A
    # is represented with two samples, both of which are
    # assembled end-to-end, which would inflate the importance
    # of haplogroup A for a high-quality assembly
    drop_one_hot_variables = [c for c in df.columns if c.startswith('is_')]
    df.drop(drop_one_hot_variables, axis=1, inplace=True)

    possible_targets = [
        'assembly_length_bp',
        'largest_contig_bp',
        'contig_NG50',
        'contigs_num'
    ]

    drop_targets = [t for t in possible_targets if t != regress_target]
    df.drop(drop_targets, axis=1, inplace=True)

    scaler = prep.StandardScaler()
    rescaled = pd.DataFrame(
        scaler.fit_transform(df),
        index=df.index,
        columns=df.columns
    )

    targets = rescaled.loc[:, regress_target].copy()
    targets.index.rename('sample', inplace=True)
    features = rescaled[[c for c in df.columns if c not in possible_targets]].copy()
    features.index.rename('sample', inplace=True)

    return df, features, targets


def main():

    target_map = {
        'assmlen': 'assembly_length_bp',
        'ctgng50': 'contig_NG50',
        'ctgnum': 'contigs_num'
    }

    args = parse_command_line()

    data_table, features, targets = load_data_table(args.data_table, target_map[args.target])
    assert pd.notnull(targets).all()
    assert pd.notnull(features).all(axis=1).all()

    l1_ratios = [0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.9999]
    l1_indices = dict((v, i) for i, v in enumerate(l1_ratios))

    en_init = model.ElasticNetCV(
        l1_ratio=l1_ratios,
        fit_intercept=False,
        max_iter=1000,
        cv=args.cv_folds,
        copy_X=True,
        n_jobs=args.num_cpu
    )

    en_fit = en_init.fit(features, targets)

    model_summary = col.OrderedDict({
        'model_type': 'elastic_net',
        'regression_target': target_map[args.target]
    })
    # value for alpha chosen by CV
    alpha = en_fit.alpha_
    model_summary['model_alpha'] = alpha
    # value for L1-ratio chosen by CV
    l1_ratio = en_fit.l1_ratio_
    model_summary['model_L1_ratio'] = l1_ratio

    for feat_name, feat_coef in zip(features.columns, en_fit.coef_):
        model_summary[f'model_coef_{feat_name}'] = feat_coef

    # this indexing is needed to extract the correct
    # values for the MSE along the regression path
    l1_idx = l1_indices[l1_ratio]
    # select all alphas for this L1 ratio
    alphas = en_fit.alphas_[l1_idx, :]
    alpha_idx = np.argwhere(np.isclose(alphas, alpha)).flatten()
    assert alpha_idx.size == 1
    alpha_idx = alpha_idx[0]

    mse_vals = en_fit.mse_path_[l1_idx, alpha_idx, :].flatten()
    assert mse_vals.size == args.cv_folds

    for fold, mse in zip(range(1, args.cv_folds + 1), mse_vals):
        model_summary[f'test_mse_fold{fold}'] = mse

    model_summary['test_mse_avg'] = np.mean(mse_vals)
    model_summary = pd.DataFrame.from_records([model_summary])
    model_summary.to_csv(args.out_model, header=True, index=False, sep='\t')

    data_table.to_csv(args.out_data, header=True, index=True, sep='\t')

    return 0


if __name__ == '__main__':
    main()
