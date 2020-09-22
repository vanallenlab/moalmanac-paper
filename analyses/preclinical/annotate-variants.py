# RUN THIS FROM moalmanac/moalmanac, instead of moalmanac/moalmanac/datasources/gdsc

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

table = 'Somatic Variant'
db = pd.DataFrame(almanac_json.table(table).all())

evidence_map = {
    'FDA-Approved': 5.0, 'Guideline': 4.0, 'Clinical trial': 3.0,
    'Clinical evidence': 2.0, 'Preclinical': 1.0, 'Inferential': 0.0}
inv_evidence_map = {float(v): k for k, v in evidence_map.items()}

db_columns = ['feature_display', 'gene', 'variant_annotation', 'protein_change', 'predictive_implication']
db['variant_annotation'].replace({'Oncogenic Mutations': '', 'Activating mutation': ''}, inplace=True)
db = db.loc[:, db_columns].drop_duplicates()
db.rename(columns={'gene': 'feature',
                   'variant_annotation': 'alteration_type',
                   'protein_change': 'alteration',
                   'predictive_implication': 'evidence'},
          inplace=True)
db['evidence_map'] = db['evidence'].replace(evidence_map)
db.sort_values(['evidence_map', 'feature_display'], ascending=[False, True], inplace=True)
db['merged'] = 1

df = pd.read_csv('{}/formatted/ccle.variants.txt'.format(root), sep='\t')

df['feature_match_1'] = 0
df['feature_match_2'] = 0
df['feature_match_3'] = 0
df['feature_match_4'] = 0
df['evidence'] = pd.NA

# Biologically Relevant, no evidence
idx_match = df['feature'].isin(almanac_genes)
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

# Gene, Dtype, AType match, and alteration match; feature_match = 4, will have evidence
match_columns = ['feature', 'alteration_type', 'alteration']
df = (df
      .merge(db[~db['alteration_type'].eq('') & ~db['alteration'].eq('')]
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

inv_map = {float(v): k for k, v in evidence_map.items()}
df['evidence_map'] = df['evidence'].copy(deep=True)
df['evidence'] = df['evidence'].fillna('').replace(inv_map)

df = annotator.CancerHotspots.annotate(df, dbs)
df = annotator.CancerHotspots3D.annotate(df, dbs)
df = annotator.CancerGeneCensus.annotate(df, dbs)
df = annotator.Cosmic.annotate(df, dbs)
df = annotator.GSEACancerPathways.annotate(df, dbs)
df = annotator.GSEACancerModules.annotate(df, dbs)

df.to_csv('{}/annotated/ccle.variants.evaluated.txt'.format(root), sep='\t', index=False)
