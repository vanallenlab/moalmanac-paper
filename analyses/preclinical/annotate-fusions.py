# RUN THIS FROM moalmanac/moalmanac, instead of moalmanac/moalmanac/datasources/gdsc
# Comment line 687 to not add addito

import time
import argparse
import pandas as pd

import annotator
import datasources
import features
import evaluator
import illustrator
import investigator
import matchmaker
import ontologymapper
import reporter
import writer

from config import COLNAMES
from config import CONFIG

snv_handle = 'snv_handle'
indel_handle = 'indel_handle'
bases_covered_handle = 'bases_covered_handle'
cnv_handle = 'cnv_handle'
fusion_handle = 'fusion_handle'
germline_handle = 'germline_handle'
validation_handle = 'validation_handle'

snv_input = 'snv_input'
indel_input = 'indel_input'
seg_input = 'seg_input'
fusion_input = 'fusion_input'
germline_input = 'germline_input'

patient_section = 'patient'
patient_id = COLNAMES[patient_section]['patient_id']
tumor_type = COLNAMES[patient_section]['tumor_type']
stage = COLNAMES[patient_section]['stage']
description = COLNAMES[patient_section]['description']
purity = COLNAMES[patient_section]['purity']
ploidy = COLNAMES[patient_section]['ploidy']
wgd = COLNAMES[patient_section]['wgd']
ms_status = COLNAMES[patient_section]['ms_status']

oncotree_section = 'oncotree'
ontology = COLNAMES[oncotree_section]['ontology']
code = COLNAMES[oncotree_section]['code']

feature_type_section = 'feature_types'
feature_type_mut = CONFIG[feature_type_section]['mut']
feature_type_germline = CONFIG[feature_type_section]['germline']
feature_type_cna = CONFIG[feature_type_section]['cna']
feature_type_fusion = CONFIG[feature_type_section]['fusion']
feature_type_burden = CONFIG[feature_type_section]['burden']
feature_type_signature = CONFIG[feature_type_section]['signature']
feature_type_microsatellite = CONFIG[feature_type_section]['microsatellite']
feature_type_aneuploidy = CONFIG[feature_type_section]['aneuploidy']
feature_types = {
    'mutation': feature_type_mut,
    'germline': feature_type_germline,
    'copynumber': feature_type_cna,
    'fusion': feature_type_fusion,
    'burden': feature_type_burden,
    'signature': feature_type_signature,
    'microsatellite': feature_type_microsatellite,
    'aneuploidy': feature_type_aneuploidy
}

parser = argparse.ArgumentParser(description='Process output directory')
parser.add_argument('--directory', type=str, help='output directory')
args = parser.parse_args()
root = args.directory

patient = {patient_id: 'GDSC', tumor_type: 'NA', stage: 'NA', description: 'NA', purity: 'NA', ploidy: 'NA',
           ms_status: 'mss', wgd: False, code: 'NA', ontology: 'NA'}

dbs = datasources.Datasources.generate_db_dict(CONFIG)
almanac_json = datasources.Almanac.import_ds(dbs)
almanac_genes = almanac_json.table('genes').all()[0]['genes']

table = 'Rearrangement'
db = pd.DataFrame(almanac_json.table(table).all())

evidence_map = {
    'FDA-Approved': 5.0, 'Guideline': 4.0, 'Clinical trial': 3.0,
    'Clinical evidence': 2.0, 'Preclinical': 1.0, 'Inferential': 0.0}
inv_evidence_map = {float(v): k for k, v in evidence_map.items()}

db_columns = ['feature_display', 'gene1', 'gene2', 'rearrangement_type', 'predictive_implication']
db = db.loc[:, db_columns].drop_duplicates()
db.rename(columns={
                   'rearrangement_type': 'alteration_type',
                   'predictive_implication': 'evidence'},
          inplace=True)
db['evidence_map'] = db['evidence'].replace(evidence_map)
db.sort_values(['evidence_map', 'feature_display'], ascending=[False, True], inplace=True)
db['merged'] = 1

df = pd.read_csv('{}/formatted/sanger.fusions.txt'.format(root), sep='\t')
df['feature_type'] = 'Rearrangement'
df['alteration_type'] = 'Fusion'
df['feature_match'] = 0
df['feature_match_1'] = 0
df['feature_match_2'] = 0
df['feature_match_3'] = 0
df['feature_match_4'] = 0
df['evidence'] = pd.NA

df_gene1 = df.copy(deep=True)
df_gene2 = df.copy(deep=True)

df_gene1['alteration'] = df['feature'] + '--' + df['partner']
df_gene2['alteration'] = df['partner'] + '--' + df['feature']
df_gene2.rename(columns={'partner': 'feature', 'feature': 'partner'}, inplace=True)

db_gene1 = db.copy(deep=True).rename(columns={'gene1': 'feature', 'gene2': 'partner'}).fillna('')
db_gene2 = db.copy(deep=True).rename(columns={'gene1': 'partner', 'gene2': 'feature'}).fillna('')


def annotate_data(df, db, genes, consider_partner=False):
    # Biologically Relevant, no evidence
    idx_match = df['feature'].isin(genes)
    df.loc[idx_match, 'feature_match_1'] = 1

    # Gene and Dtype match; feature_match = 2, will have evidence
    match_columns = ['feature']
    df = (df
          .merge(db
                 .loc[:, match_columns + ['evidence_map', 'merged']]
                 .sort_values('evidence_map', ascending=False)
                 .drop_duplicates(match_columns, keep='first'),
                 on=match_columns,
                 how='left'
                 )
          )
    index_match = df['merged'].eq(1)
    df.loc[index_match, 'feature_match_2'] = 1
    df.loc[index_match, 'evidence'] = df.loc[index_match, 'evidence_map']
    df.drop(['merged', 'evidence_map'], axis=1, inplace=True)

    # Gene, Dtype, AType match; feature_match = 3, will have evidence
    match_columns = ['feature', 'alteration_type']
    df = (df
          .merge(db[~db['alteration_type'].eq('')]
                 .loc[:, match_columns + ['evidence_map', 'merged']]
                 .sort_values('evidence_map', ascending=False)
                 .drop_duplicates(match_columns, keep='first'),
                 on=match_columns,
                 how='left'
                 )
          )
    index_match = df['merged'].eq(1)
    df.loc[index_match, 'feature_match_3'] = 1
    df.loc[index_match, 'evidence'] = df.loc[index_match, 'evidence_map']
    df.drop(['merged', 'evidence_map'], axis=1, inplace=True)

    if consider_partner:
        # Gene, Dtype, AType match, and alteration match; feature_match = 4, will have evidence
        match_columns = ['feature', 'alteration_type', 'partner']
        df = (df
          .merge(db[~db['alteration_type'].eq('') & ~db['partner'].eq('')]
                 .loc[:, match_columns + ['evidence_map', 'merged']]
                 .sort_values('evidence_map', ascending=False)
                 .drop_duplicates(match_columns, keep='first'),
                 on=match_columns,
                 how='left'
                 )
          )
        index_match = df['merged'].eq(1)
        df.loc[index_match, 'feature_match_4'] = 1
        df.loc[index_match, 'evidence'] = df.loc[index_match, 'evidence_map']
        df.drop(['merged', 'evidence_map'], axis=1, inplace=True)

    feature_match_columns = ['feature_match_1', 'feature_match_2', 'feature_match_3', 'feature_match_4']
    df['feature_match'] = df.loc[:, feature_match_columns].sum(axis=1)

    inv_map = {float(v): k for k, v in evidence_map.items()}
    df['evidence_map'] = df['evidence'].copy(deep=True)
    df['evidence'] = df['evidence'].fillna('').replace(inv_map)

    return df


# Gene1 and Gene 2
group1 = annotate_data(df_gene1, db_gene1, almanac_genes, True)
group2 = annotate_data(df_gene1, db_gene2, almanac_genes, True)
group3 = annotate_data(df_gene2, db_gene1, almanac_genes, True)
group4 = annotate_data(df_gene2, db_gene2, almanac_genes, True)

values = pd.concat([
    group1.rename(columns={'evidence_map': 'group1'})['group1'],
    group2.rename(columns={'evidence_map': 'group2'})['group2'],
    group3.rename(columns={'evidence_map': 'group3'})['group3'],
    group4.rename(columns={'evidence_map': 'group4'})['group4'],
], axis=1)
values = values.fillna(-1).idxmax(axis=1)

idx_group1 = values[values.eq('group1')].index
idx_group2 = values[values.eq('group2')].index
idx_group3 = values[values.eq('group3')].index
idx_group4 = values[values.eq('group4')].index

columns = ['evidence', 'feature_match_1', 'feature_match_2', 'feature_match_3', 'feature_match_4', 'feature_match']
pairs = [(group1, idx_group1), (group2, idx_group2), (group3, idx_group3), (group4, idx_group4)]

df['which_match'] = ''
for index in [idx_group1, idx_group2]:
    df.loc[index, 'which_match'] = df.loc[index, 'feature']
for index in [idx_group3, idx_group4]:
    df.loc[index, 'which_match'] = df.loc[index, 'partner']

for group, index in pairs:
    df.loc[index, columns] = group.loc[index, columns]

columns = ['feature', 'partner']
for index in df.index:
    df.loc[index, columns] = sorted(df.loc[index, columns].astype(str).tolist())

df = (df
      .sort_values(['evidence', 'feature_match'], ascending=False)
      .drop_duplicates(['feature', 'partner', 'model_id'], keep='first')
      )

df['alteration'] = ''
df = annotator.CancerHotspots.annotate(df, dbs)
df = annotator.CancerHotspots3D.annotate(df, dbs)
df = annotator.CancerGeneCensus.annotate(df, dbs)
df = annotator.Cosmic.annotate(df, dbs)
df = annotator.GSEACancerPathways.annotate(df, dbs)
df = annotator.GSEACancerModules.annotate(df, dbs)

df.to_csv('{}/annotated/sanger.fusions.evaluated.txt'.format(root), sep='\t', index=False)


# Gene 1
group1 = annotate_data(df_gene1, db_gene1, almanac_genes, False)
group2 = annotate_data(df_gene1, db_gene2, almanac_genes, False)

values = pd.concat([
    group1.rename(columns={'evidence_map': 'group1'})['group1'],
    group2.rename(columns={'evidence_map': 'group2'})['group2'],
], axis=1)
values = values.fillna(-1).idxmax(axis=1)

idx_group1 = values[values.eq('group1')].index
idx_group2 = values[values.eq('group2')].index

columns = ['evidence', 'feature_match']
group1.loc[idx_group2, columns] = group2.loc[idx_group2, columns]

group1 = (group1
          .sort_values(['evidence', 'feature_match'], ascending=False)
          .drop_duplicates(['feature', 'model_id'], keep='first')
          )

group1['alteration'] = ''
group1 = annotator.CancerHotspots.annotate(group1, dbs)
group1 = annotator.CancerHotspots3D.annotate(group1, dbs)
group1 = annotator.CancerGeneCensus.annotate(group1, dbs)
group1 = annotator.Cosmic.annotate(group1, dbs)
group1 = annotator.GSEACancerPathways.annotate(group1, dbs)
group1 = annotator.GSEACancerModules.annotate(group1, dbs)

group1.to_csv('{}/annotated/sanger.fusions.gene1.evaluated.txt'.format(root), sep='\t', index=False)

# Gene 2
group3 = annotate_data(df_gene2, db_gene1, almanac_genes, False)
group4 = annotate_data(df_gene2, db_gene2, almanac_genes, False)

values = pd.concat([
    group3.rename(columns={'evidence_map': 'group3'})['group3'],
    group4.rename(columns={'evidence_map': 'group4'})['group4'],
], axis=1)
values = values.fillna(-1).idxmax(axis=1)

idx_group3 = values[values.eq('group3')].index
idx_group4 = values[values.eq('group4')].index

columns = ['evidence', 'feature_match']
group3.loc[idx_group4, columns] = group4.loc[idx_group4, columns]

group3 = (group3
          .sort_values(['evidence', 'feature_match'], ascending=False)
          .drop_duplicates(['feature', 'model_id'], keep='first')
          )

group3['alteration'] = ''
group3 = annotator.CancerHotspots.annotate(group3, dbs)
group3 = annotator.CancerHotspots3D.annotate(group3, dbs)
group3 = annotator.CancerGeneCensus.annotate(group3, dbs)
group3 = annotator.Cosmic.annotate(group3, dbs)
group3 = annotator.GSEACancerPathways.annotate(group3, dbs)
group3 = annotator.GSEACancerModules.annotate(group3, dbs)

group3.to_csv('{}/annotated/sanger.fusions.gene2.evaluated.txt'.format(root), sep='\t', index=False)
