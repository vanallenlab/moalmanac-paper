import argparse
import numpy as np
import pandas as pd
import pickle

import models as models
from metrics import Metrics
from plots import AveragePrecision, AveragePrecisionK

case = 'case'
comparison = 'comparison'
n_sensitive_union = 'n_sensitive_union'
n_shared = 'n_shared'


def format_gdsc_pairs(dataframe):
    return (dataframe
                .sort_values([case, comparison], ascending=[True, False])
                .set_index([case, comparison])
                .rename(columns={n_sensitive_union: n_shared})
                .loc[:, n_shared])


def format_inputs_sample_names(inputs, names_to_use, samples):
    for data_type in ['variants', 'copy_number_alterations']:
        tmp = inputs[data_type]
        if names_to_use != 'ccle_name':
            sample_map = summary.loc[:, ['ccle_name', names_to_use]].dropna().set_index('ccle_name')[names_to_use].to_dict()
            tmp['model_id'] = tmp['model_id'].replace(sample_map)
        inputs[data_type] = tmp[tmp['model_id'].isin(samples)].reset_index(drop=True)
    for data_type in ['fusions', 'fusions_gene1', 'fusions_gene2', 'gdsc']:
        tmp = inputs[data_type]
        if names_to_use != 'sanger':
            sample_map = summary.loc[:, [names_to_use, 'sanger']].dropna().set_index('sanger')[names_to_use].to_dict()
            tmp['model_id'] = tmp['model_id'].replace(sample_map)
        inputs[data_type] = tmp[tmp['model_id'].isin(samples_to_use)].reset_index(drop=True)
    for data_type in ['gdsc_pairwise']:
        tmp = inputs[data_type]
        if names_to_use != 'sanger':
            sample_map = summary.loc[:, [names_to_use, 'sanger']].dropna().set_index('sanger')[names_to_use].to_dict()
            tmp['case'] = tmp['case'].replace(sample_map)
            tmp['comparison'] = tmp['comparison'].replace(sample_map)
        tmp = tmp[tmp['case'].isin(samples) & tmp['comparison'].isin(samples)].reset_index(drop=True)
        inputs[data_type] = tmp
    return inputs


def write_pickle(handle, output):
    file = open(handle, 'wb')
    pickle.dump(output, file)
    file.close()


def main(inputs, samples):
    np.random.seed(seed=42)

    models_list = [
        models.AlmanacGenes,
        models.AlmanacFeatureTypes,
        models.AlmanacFeatures,
        models.CGC,
        models.CGCFeatureTypes,
        models.Compatibility,
        models.NonsynVariantCount,
        models.PCAonAlmanac,
        models.PCAonCGC,
        models.RankedSortAlmanacEvidenceCGC,
        models.SNFbyEvidenceCGC,
        models.SNFTypesCGC,
        models.SNFTypesCGCwithEvidence,
        models.SNFTypesAlmanac,
        models.Tree
    ]

    calculated = [model.calculate(inputs, samples) for model in models_list]
    model_names = [model.label for model in models_list]
    model_descriptions = {}
    for model in models_list:
        model_descriptions[model.label] = model.description

    distances = pd.concat(calculated, axis=1)
    labels = format_gdsc_pairs(inputs['gdsc_pairwise'])

    distances.loc[distances.index, n_shared] = labels.loc[distances.index]
    labeled = distances
    labeled.to_csv('tables/distances/models.labeled.txt', sep='\t')

    evaluated_models_dictionary = Metrics.evaluate_models(samples, labeled, model_names, model_descriptions)
    write_pickle('tables/models.evaluated.pkl', evaluated_models_dictionary)
    AveragePrecision.plot(evaluated_models_dictionary, model_names)
    AveragePrecisionK.plot(evaluated_models_dictionary, model_names)


if __name__ == "__main__":
    moalmanac_method_datasources = '/Users/brendan/Github/moalmanac/moalmanac/datasources'
    
    arg_parser = argparse.ArgumentParser(prog='Evaluate matchmaking models for Molecular Oncology Almanac',
                                         description='Format and annotate CCLE for use with Molecular Oncology Almanac')
    arg_parser.add_argument('--variants',
                            help='File handle to annotated somatic variants',
                            default='annotated/ccle.variants.evaluated.txt')
    arg_parser.add_argument('--copy_number_alterations',
                            help='File handle to annotated copy number alterations',
                            default='annotated/ccle.copy-numbers.evaluated.txt')
    arg_parser.add_argument('--fusions',
                            help='File handle to annotated fusions',
                            default='annotated/sanger.fusions.evaluated.txt')
    arg_parser.add_argument('--fusions_gene1',
                            help='File handle to annotated fusions, just gene1',
                            default='annotated/sanger.fusions.gene1.evaluated.txt')
    arg_parser.add_argument('--fusions_gene2',
                            help='File handle to annotated fusions, just gene2',
                            default='annotated/sanger.fusions.gene2.evaluated.txt')
    arg_parser.add_argument('--summary',
                            help='File handle to model information and which samples to use',
                            default='formatted/cell-lines.summary.txt')
    arg_parser.add_argument('--gdsc',
                            help='File handle GDSC therapeutic screen data',
                            default='formatted/sanger.gdsc.txt')
    arg_parser.add_argument('--gdsc_pairwise',
                            help='File handle GDSC pairwise sensitivities',
                            default='formatted/sanger.gdsc.pairwise-sensitive.txt')
    arg_parser.add_argument('--almanac',
                            help='File handle to molecular oncology almanac',
                            default=f'{moalmanac_method_datasources}/moalmanac/moalmanac.json')
    arg_parser.add_argument('--cgc',
                            help='File handle to cancer gene census', 
                            default=f'{moalmanac_method_datasources}/cancergenecensus/cancer_gene_census_v85.tsv')
    arg_parser.add_argument('--therapy_mappings',
                            help='File handle to mappings of Molecular Oncology Almanac therapy to GDSC.',
                            default='almanac-gdsc-mappings.json')
    arg_parser.add_argument('--naming',
                            help='Naming convention for sample names.',
                            default='broad',
                            choices=['broad', 'sanger', 'ccle_name'])
    args = arg_parser.parse_args()

    input_handles = {
        'variants': args.variants,
        'copy_number_alterations': args.copy_number_alterations,
        'fusions': args.fusions,
        'fusions_gene1': args.fusions_gene1,
        'fusions_gene2': args.fusions_gene2,
        'summary': args.summary,
        'gdsc': args.gdsc,
        'gdsc_pairwise': args.gdsc_pairwise,
        'almanac': args.almanac,
        'cgc': args.cgc,
        'therapy_mappings': args.therapy_mappings,
    }

    inputs_dictionary = {}
    dtypes_txt = ['variants', 'copy_number_alterations', 'fusions', 'fusions_gene1', 'fusions_gene2',
                  'summary', 'gdsc', 'gdsc_pairwise', 'cgc']
    for dtype in dtypes_txt:
        inputs_dictionary[dtype] = pd.read_csv(input_handles[dtype], sep='\t')
    inputs_dictionary['almanac'] = input_handles['almanac']

    sample_col_to_use = args.naming
    summary = inputs_dictionary['summary']
    samples_to_use = summary[(summary['use_evaluate'].eq(1))][sample_col_to_use].sort_values().tolist()
    summary['model_id'] = summary[sample_col_to_use]
    inputs_dictionary['summary'] = summary[summary[sample_col_to_use].isin(samples_to_use)]
    inputs_dictionary = format_inputs_sample_names(inputs_dictionary, sample_col_to_use, samples_to_use)

    main(inputs_dictionary, samples_to_use)
