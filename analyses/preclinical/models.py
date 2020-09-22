import pandas as pd
from sklearn import neighbors
import sklearn
import snf
import tinydb


class Models:
    variants = 'variants'
    cnas = 'copy_number_alterations'
    fusions = 'fusions'
    fusions_gene1 = 'fusions_gene1'
    fusions_gene2 = 'fusions_gene2'
    summary = 'summary'

    feature = 'feature'
    feature_type = 'feature_type'
    model_id = 'model_id'

    case = 'case'
    comparison = 'comparison'

    variant = 'Somatic Variant'
    copy_number = 'Copy Number'
    rearrangement = 'Rearrangement'

    jaccard = neighbors.DistanceMetric.get_metric('jaccard')

    @staticmethod
    def aggregate_columns_by_row(dataframe, delimiter):
        return dataframe.agg(delimiter.join, axis=1)

    @staticmethod
    def calculate_distance(dataframe, metric):
        distance = metric(dataframe.values)
        return pd.DataFrame(distance, columns=dataframe.index.tolist(), index=dataframe.index.tolist())

    @classmethod
    def format_label_for_output(cls, string):
        return string.lower().replace(': ', '-').replace(' ', '_')

    @staticmethod
    def reset_multi_indexed_dataframe(dataframe, remapped_labels):
        return dataframe.reset_index().rename(columns=remapped_labels)

    @classmethod
    def series_concat(cls, list_of_series, drop_duplicates=True):
        series = pd.concat(list_of_series)
        return cls.series_to_list(series, drop_duplicates)

    @staticmethod
    def series_to_list(series, drop_duplicates=True):
        if drop_duplicates:
            return series.drop_duplicates().dropna().sort_values().tolist()
        else:
            return series.dropna().sort_values().tolist()

    @classmethod
    def stack_distances(cls, dataframe, label):
        # Commented out sorting and dropping of duplicates to leave all N x N comparisons
        stacked = dataframe.stack().reset_index()
        # stacked.drop_duplicates(subset=[cls.case, cls.comparison], keep='first', inplace=True)
        stacked.rename(columns={'level_0': cls.case, 'level_1': cls.comparison, 0: label}, inplace=True)

        stacked.to_csv('tables/distances/{}.stacked.txt'.format(label), sep='\t')
        stacked.set_index([cls.case, cls.comparison], inplace=True)
        return stacked.loc[stacked.index, label]


class Almanac(Models):
    almanac = 'almanac'
    genes = 'genes'
    gene = 'gene'
    gene1 = 'gene1'
    gene2 = 'gene2'
    variant_annotation = 'variant_annotation'
    protein_change = 'protein_change'
    direction = 'direction'

    alteration_type = 'alteration_type'
    alteration = 'alteration'
    partner = 'partner'

    missense = 'Missense'
    truncating_types = ['Nonsense', 'Nonstop', 'Frameshift', 'Splice Site']
    truncating = 'Truncating'
    fusion = 'Fusion'

    variant_missense = 'variant_missense'
    variant_truncating = 'variant_truncating'

    feature_match = 'feature_match'
    feature_match_1 = 'feature_match_2'
    feature_match_2 = 'feature_match_2'
    feature_match_3 = 'feature_match_3'
    feature_match_4 = 'feature_match_4'

    string = 'feature_string'

    @classmethod
    def import_dbs(cls, inputs):
        db = tinydb.TinyDB(inputs[cls.almanac])
        tables = {}
        for table in [cls.variant, cls.copy_number, cls.rearrangement]:
            tmp = pd.DataFrame(db.table(table).all())
            tables[table] = tmp
        tables[cls.genes] = pd.Series(db.table(cls.genes).all()[0][cls.genes])
        return tables

    @classmethod
    def generate_features(cls, dbs):
        missense = cls.generate_features_missense(dbs[cls.variant], cls.gene, cls.variant_annotation,
                                                  cls.protein_change)
        # cls.generate_features_missense_aggregated(dbs[cls.variant], cls.gene, cls.variant_annotation, cls.protein_change)
        truncating = cls.generate_features_truncating_aggregated(dbs[cls.variant], cls.gene, cls.variant_annotation)
        # cls.generate_features_truncating(dbs[cls.variant], cls.gene, cls.variant_annotation)
        nonspecific_variant = cls.generate_features_nonspecific_variant(dbs[cls.variant])
        copy_number = cls.generate_features_copy_number(dbs[cls.copy_number], cls.gene, cls.direction)
        fusions = cls.generate_features_fusions(dbs[cls.rearrangement])
        fusion_partners = cls.generate_features_fusions_partners(dbs[cls.rearrangement])
        series = (pd
                  .Index(missense)
                  .union(pd.Index(truncating))
                  .union(pd.Index(nonspecific_variant))
                  .union(pd.Index(copy_number))
                  .union(pd.Index(fusions))
                  .union(pd.Index(fusion_partners)))
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_copy_number(cls, db, column_gene, column_type):
        series = cls.aggregate_columns_by_row(db.loc[:, [column_gene, column_type]], ' ')
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_nonspecific_variant(cls, db):
        idx = (db[cls.variant_annotation]
               .replace({'Oncogenic Mutations': '', 'Activating mutation': ''})
               .fillna('')
               .eq('')
               )
        series = db.loc[idx, cls.gene] + f' {cls.variant}'
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_missense(cls, db, column_gene, column_type, column_protein):
        idx = db[column_type].eq(cls.missense) & ~db[column_protein].fillna('').eq('')
        series = cls.aggregate_columns_by_row(db.loc[idx, [column_gene, column_protein]], ' ')
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_missense_aggregated(cls, db, column_gene, column_type, column_protein):
        idx = db[column_type].eq(cls.missense) & ~db[column_protein].fillna('').eq('')
        series = db.loc[idx, column_gene] + f' {cls.missense}'
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_truncating(cls, db, column_gene, column_type):
        idx = db[column_type].isin(cls.truncating_types)
        series = cls.aggregate_columns_by_row(db.loc[idx, [column_gene, column_type]], ' ')
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_truncating_aggregated(cls, db, column_gene, column_type):
        idx = db[column_type].isin(cls.truncating_types)
        series = db.loc[idx, column_gene] + f' {cls.truncating}'
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_fusions(cls, db):
        columns = [cls.gene1, cls.gene2]
        paired = cls.sort_fusions(db.loc[:, columns].dropna(), cls.gene1, cls.gene2)
        series = cls.aggregate_columns_by_row(paired.loc[:, columns], '--')
        return cls.series_to_list(series, drop_duplicates=True)

    @classmethod
    def generate_features_fusions_partners(cls, db):
        series1 = db.loc[db[cls.gene2].fillna('').eq(''), cls.gene1] + f' {cls.fusion}'
        series2 = db.loc[db[cls.gene2].fillna('').eq(''), cls.gene2] + f' {cls.fusion}'
        return cls.series_to_list(pd.concat([series1, series2]), drop_duplicates=True)

    @classmethod
    def generate_gene_features(cls, db):
        return cls.series_to_list(db[cls.genes])

    @classmethod
    def generate_gene_features_dtype(cls, dbs):
        variants = dbs[cls.variant][cls.gene].drop_duplicates() + f' {cls.variant}'
        copy_number = dbs[cls.copy_number][cls.gene].drop_duplicates() + f' {cls.copy_number}'
        rearrangements = (pd
                          .concat([dbs[cls.rearrangement][cls.gene1], dbs[cls.rearrangement][cls.gene2]])
                          .drop_duplicates()
                          .dropna()
                          .sort_values()
                          ) + f' {cls.rearrangement}'
        return cls.series_concat([variants, copy_number, rearrangements])

    @classmethod
    def sort_fusions(cls, df, gene1_col, gene2_col):
        for index in df.index:
            df.loc[index, [gene1_col, gene2_col]] = sorted(df.loc[index, [gene1_col, gene2_col]].tolist())
        return df


class AlmanacFeatures(Almanac):
    label = 'jaccard-almanac-features'
    label_output = Models.format_label_for_output(label)
    description = 'We sort by agreement based measure (jaccard) by considering all somatic variant, copy number, and ' \
                  'rearrangement molecular features catalogued in the Molecular Oncology Almanac.'

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        almanac = cls.import_dbs(input_dtypes)
        boolean_dataframe = cls.create_boolean_table(input_dtypes, samples_list, almanac)
        distance_dataframe = cls.calculate_distance(boolean_dataframe, cls.jaccard.pairwise)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        boolean_dataframe.to_csv('tables/features/{}.boolean.txt'.format(cls.label_output), sep='\t')
        return stacked_dataframe

    @classmethod
    def create_bool(cls, dataframe, dtype_column):
        df = dataframe.loc[:, [cls.string, cls.model_id]].drop_duplicates()
        df[dtype_column] = 1
        return df.set_index([cls.string, cls.model_id])

    @classmethod
    def create_boolean_table(cls, inputs, samples, almanac):
        features = cls.generate_features(almanac)
        index = pd.MultiIndex.from_product([features, samples])
        columns = [cls.variant_missense, cls.variant_truncating, 'variant_nonspecific',
                   cls.copy_number,
                   cls.fusions, cls.fusions_gene1, cls.fusions_gene2
                   ]
        df = pd.DataFrame(index=index, columns=columns)

        variants = inputs[cls.variants]
        missense = cls.subset_missense(variants)
        truncating = cls.subset_truncating(variants)
        nonspecific = cls.subset_nonspecific_variant(variants, almanac)
        copy_number_alterations = cls.subset_copy_number(inputs[cls.cnas])
        fusions = cls.subset_fusions(inputs[cls.fusions])
        fusions_gene1 = cls.subset_fusion_members(inputs[cls.fusions_gene1])
        fusions_gene2 = cls.subset_fusion_members(inputs[cls.fusions_gene2])

        missense[cls.string] = cls.aggregate_columns_by_row(missense.loc[:, [cls.feature, cls.alteration]], ' ')
        # missense[cls.feature] + ' Missense'
        truncating[cls.string] = truncating[cls.feature] + f' {cls.truncating}'
        # cls.aggregate_columns_by_row(truncating.loc[:, [cls.feature, cls.alteration_type]], ' ')
        nonspecific[cls.string] = nonspecific[cls.feature] + f' {cls.variant}'
        copy_number_alterations[cls.string] = cls.aggregate_columns_by_row(
            copy_number_alterations.loc[:, [cls.feature, cls.alteration_type]], ' ')
        fusions[cls.string] = cls.aggregate_columns_by_row(fusions.loc[:, [cls.feature, cls.partner]], '--')
        fusions_gene1[cls.string] = fusions_gene1[cls.feature] + f' {cls.fusion}'
        fusions_gene2[cls.string] = fusions_gene2[cls.feature] + f' {cls.fusion}'

        df[cls.variant_missense] = cls.create_bool(missense, cls.variant_missense)
        df[cls.variant_truncating] = cls.create_bool(truncating, cls.variant_truncating)
        df['variant_nonspecific'] = cls.create_bool(nonspecific, 'variant_nonspecific')
        df[cls.copy_number] = cls.create_bool(copy_number_alterations, cls.copy_number)
        df[cls.fusions] = cls.create_bool(fusions, cls.fusions)
        df[cls.fusions_gene1] = cls.create_bool(fusions_gene1, cls.fusions_gene1)
        df[cls.fusions_gene2] = cls.create_bool(fusions_gene2, cls.fusions_gene2)

        df = cls.reset_multi_indexed_dataframe(df.fillna(0), {'level_0': cls.string, 'level_1': cls.model_id})
        df[cls.label] = df.loc[:, columns].sum(axis=1).astype(bool).astype(int)
        return df.pivot_table(index=cls.model_id, columns=cls.string, values=cls.label)

    @classmethod
    def subset_copy_number(cls, df):
        return df[df[cls.feature_match_3].eq(1)].reset_index(drop=True)

    @classmethod
    def subset_fusions(cls, df):
        dataframe = df[df[cls.feature_match_4].eq(1)]
        return cls.sort_fusions(dataframe, cls.feature, cls.partner).reset_index(drop=True)

    @classmethod
    def subset_fusion_members(cls, df):
        return df[df[cls.feature_match_2].eq(1)].reset_index(drop=True)

    @classmethod
    def subset_nonspecific_variant(cls, df, db):
        table = db[cls.variant]
        idx1 = table['variant_annotation'].replace({'Oncogenic Mutations': '', 'Activating mutation': ''}).fillna(
            '').eq('')
        idx2 = table[cls.protein_change].fillna('').eq('')
        relevant_genes = table.loc[(idx1 & idx2), cls.gene]
        return df[df[cls.feature_match_2].eq(1) & df[cls.feature].isin(relevant_genes)].reset_index(drop=True)

    @classmethod
    def subset_missense(cls, df):
        return df[(
                df[cls.feature_match_4].eq(1)
                & df[cls.alteration_type].eq(cls.missense)
                & ~df[cls.alteration_type].fillna('').eq('')
        )].reset_index(drop=True)

    @classmethod
    def subset_truncating(cls, df):
        return df[df[cls.feature_match_3].eq(1) & df[cls.alteration_type].isin(cls.truncating_types)].reset_index(
            drop=True)


class AlmanacFeatureTypes(Almanac):
    label = 'jaccard-almanac-feature-types'
    label_output = Models.format_label_for_output(label)
    description = 'We sort by agreement based measure (jaccard) by considering both gene and data type for all ' \
                  'somatic variants, copy number alterations, and rearrangements catalogued in the ' \
                  'Molecular Oncology Almanac (e.g. CDKN2A copy number alterations match but not a ' \
                  'CDKN2A deletion and CDKN2A nonsense somatic variant).'

    @classmethod
    def create_bool(cls, dataframe, value_column, dtype_column, string):
        df = dataframe[dataframe[value_column].astype(int).eq(1)].loc[:, [cls.feature, cls.model_id]].drop_duplicates()
        df[cls.feature] = df[cls.feature] + string
        df[dtype_column] = 1
        return df.set_index([cls.feature, cls.model_id])

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        almanac = cls.import_dbs(input_dtypes)
        boolean_dataframe = cls.create_boolean_table(input_dtypes, samples_list, almanac)
        distance_dataframe = cls.calculate_distance(boolean_dataframe, cls.jaccard.pairwise)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        boolean_dataframe.to_csv('tables/features/{}.boolean.txt'.format(cls.label_output), sep='\t')
        return stacked_dataframe

    @classmethod
    def create_boolean_table(cls, inputs, samples, almanac):
        features = cls.generate_gene_features_dtype(almanac)
        index = pd.MultiIndex.from_product([features, samples])

        variants = inputs[cls.variants]
        copy_number_alterations = inputs[cls.cnas]
        fusions_gene1 = inputs[cls.fusions_gene1]
        fusions_gene2 = inputs[cls.fusions_gene2]

        df = pd.DataFrame(columns=[cls.variants, cls.cnas, cls.fusions_gene1, cls.fusions_gene2], index=index)
        df[cls.variants] = cls.create_bool(variants, cls.feature_match_2, cls.variants, ' {}'.format(cls.variant))
        df[cls.cnas] = cls.create_bool(copy_number_alterations, cls.feature_match_2, cls.cnas,
                                       ' {}'.format(cls.copy_number))
        df[cls.fusions_gene1] = cls.create_bool(fusions_gene1, cls.feature_match_2, cls.fusions_gene1,
                                                ' {}'.format(cls.rearrangement))
        df[cls.fusions_gene2] = cls.create_bool(fusions_gene2, cls.feature_match_2, cls.fusions_gene2,
                                                ' {}'.format(cls.rearrangement))
        df = cls.reset_multi_indexed_dataframe(df.fillna(0), {'level_0': cls.feature, 'level_1': cls.model_id})

        sum_columns = [cls.variants, cls.cnas, cls.fusions_gene1, cls.fusions_gene2]
        df[cls.label] = df.loc[:, sum_columns].sum(axis=1).astype(bool).astype(int)
        return df.pivot_table(index=cls.model_id, columns=cls.feature, values=cls.label)


class AlmanacGenes(Almanac):
    label = 'jaccard-almanac-genes'
    label_output = Models.format_label_for_output(label)
    description = 'We sort by agreement based measure (jaccard) by considering any somatic variant, ' \
                  'copy number alteration, and rearrangement in any gene catalogued in Molecular Oncology Almanac.'

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        almanac = cls.import_dbs(input_dtypes)
        boolean_dataframe = cls.create_boolean_table(input_dtypes, samples_list, almanac)
        distance_dataframe = cls.calculate_distance(boolean_dataframe, cls.jaccard.pairwise)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        boolean_dataframe.to_csv('tables/features/{}.boolean.txt'.format(cls.label_output), sep='\t')
        return stacked_dataframe

    @classmethod
    def create_bool(cls, dataframe, value_column, dtype_column):
        df = dataframe[dataframe[value_column].astype(int).eq(1)].loc[:, [cls.feature, cls.model_id]].drop_duplicates()
        df[dtype_column] = 1
        return df.set_index([cls.feature, cls.model_id])

    @classmethod
    def create_boolean_table(cls, inputs, samples, almanac):
        features = cls.generate_gene_features(almanac)
        index = pd.MultiIndex.from_product([features, samples])

        variants = inputs[cls.variants]
        copy_number_alterations = inputs[cls.cnas]
        fusions_gene1 = inputs[cls.fusions_gene1]
        fusions_gene2 = inputs[cls.fusions_gene2]

        df = pd.DataFrame(columns=[cls.variants, cls.cnas, cls.fusions_gene1, cls.fusions_gene2], index=index)
        df[cls.variants] = cls.create_bool(variants, cls.feature_match_1, cls.variants)
        df[cls.cnas] = cls.create_bool(copy_number_alterations, cls.feature_match_1, cls.cnas)
        df[cls.fusions_gene1] = cls.create_bool(fusions_gene1, cls.feature_match_1, cls.fusions_gene1)
        df[cls.fusions_gene2] = cls.create_bool(fusions_gene2, cls.feature_match_1, cls.fusions_gene2)
        df = cls.reset_multi_indexed_dataframe(df.fillna(0), {'level_0': cls.feature, 'level_1': cls.model_id})

        sum_columns = [cls.variants, cls.cnas, cls.fusions_gene1, cls.fusions_gene2]
        df[cls.label] = df.loc[:, sum_columns].sum(axis=1).astype(bool).astype(int)
        return df.pivot_table(index=cls.model_id, columns=cls.feature, values=cls.label)


class AlmanacEvidence(Almanac):
    label = 'almanac-by-evidence'
    predictive_implication = 'predictive_implication'

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        almanac = cls.import_dbs(input_dtypes)
        almanac_subset = cls.subset_almanac_by_evidence(almanac, ['FDA-Approved'])
        boolean_dataframe = AlmanacFeatures.create_boolean_table(input_dtypes, samples_list, almanac_subset)
        distance_dataframe = AlmanacFeatures.calculate_distance(boolean_dataframe, cls.jaccard.pairwise)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)
        return stacked_dataframe

    @classmethod
    def subset_almanac_by_evidence(cls, db, evidence_list):
        for table in [cls.variant, cls.copy_number, cls.rearrangement]:
            db_table = db[table]
            db[table] = db_table[db_table[cls.predictive_implication].isin(evidence_list)]
        return db


class CGC(Models):
    label = 'jaccard-cgc-genes'
    label_output = Models.format_label_for_output(label)
    input_label = 'cgc'
    description = 'We sort by agreement based measure (jaccard) by considering any variant in a ' \
                  'Cancer Gene Census gene.'

    cgc_bin = 'cgc_bin'
    gene = 'Gene Symbol'

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        boolean_dataframe = cls.create_boolean_table(input_dtypes, samples_list)
        distance_dataframe = cls.calculate_distance(boolean_dataframe, cls.jaccard.pairwise)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        boolean_dataframe.to_csv('tables/features/{}.boolean.txt'.format(cls.label_output), sep='\t')
        return stacked_dataframe

    @classmethod
    def create_bool(cls, dataframe, value_column, dtype_column):
        df = dataframe[dataframe[value_column].astype(int).eq(1)].loc[:, [cls.feature, cls.model_id]].drop_duplicates()
        df[dtype_column] = 1
        return df.set_index([cls.feature, cls.model_id])

    @classmethod
    def create_boolean_table(cls, inputs, samples):
        cgc = inputs[CGC.input_label]
        features = cls.create_features_list(cgc)

        variants = inputs[cls.variants]
        copy_number_alterations = inputs[cls.cnas]
        fusions_gene1 = inputs[cls.fusions_gene1]
        fusions_gene2 = inputs[cls.fusions_gene2]

        index = pd.MultiIndex.from_product([features, samples])
        df = pd.DataFrame(columns=[cls.variants, cls.cnas, cls.fusions_gene1, cls.fusions_gene2], index=index)
        df[cls.variants] = cls.create_bool(variants, cls.cgc_bin, cls.variants)
        df[cls.cnas] = cls.create_bool(copy_number_alterations, cls.cgc_bin, cls.cnas)
        df[cls.fusions_gene1] = cls.create_bool(fusions_gene1, cls.cgc_bin, cls.fusions_gene1)
        df[cls.fusions_gene2] = cls.create_bool(fusions_gene2, cls.cgc_bin, cls.fusions_gene2)
        df = cls.reset_multi_indexed_dataframe(df.fillna(0), {'level_0': cls.feature, 'level_1': cls.model_id})

        sum_columns = [cls.variants, cls.cnas, cls.fusions_gene1, cls.fusions_gene2]
        df[cls.label] = df.loc[:, sum_columns].sum(axis=1).astype(bool).astype(int)
        return df.pivot_table(index=cls.model_id, columns=cls.feature, values=cls.label)

    @classmethod
    def create_features_list(cls, db):
        return cls.series_to_list(db[cls.gene], drop_duplicates=True)


class CGCFeatureTypes(Models):
    label = 'jaccard-cgc-feature-types'
    label_output = Models.format_label_for_output(label)
    description = 'We sort by agreement based measure (jaccard) by considering variants in a ' \
                  'Cancer Gene Census gene and feature type (e.g. CDKN2A copy number alterations match but not a ' \
                  'CDKN2A deletion and CDKN2A nonsense somatic variant).'

    cgc_bin = 'cgc_bin'
    gene = 'Gene Symbol'

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        boolean_dataframe = cls.create_boolean_table(input_dtypes, samples_list)
        distance_dataframe = cls.calculate_distance(boolean_dataframe, cls.jaccard.pairwise)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        boolean_dataframe.to_csv('tables/features/{}.boolean.txt'.format(cls.label_output), sep='\t')
        return stacked_dataframe

    @classmethod
    def create_bool(cls, dataframe, value_column, dtype_column, string):
        df = dataframe[dataframe[value_column].astype(int).eq(1)].loc[:, [cls.feature, cls.model_id]].drop_duplicates()
        df[cls.feature] = df[cls.feature] + string
        df[dtype_column] = 1
        return df.set_index([cls.feature, cls.model_id])

    @classmethod
    def create_boolean_table(cls, inputs, samples):
        cgc = inputs[CGC.input_label]
        features = cls.create_features_list(cgc)

        variants = inputs[cls.variants]
        copy_number_alterations = inputs[cls.cnas]
        fusions_gene1 = inputs[cls.fusions_gene1]
        fusions_gene2 = inputs[cls.fusions_gene2]

        index = pd.MultiIndex.from_product([features, samples])
        df = pd.DataFrame(columns=[cls.variants, cls.cnas, cls.fusions_gene1, cls.fusions_gene2], index=index)
        df[cls.variants] = cls.create_bool(variants, cls.cgc_bin, cls.variants, f' {cls.variant}')
        df[cls.cnas] = cls.create_bool(copy_number_alterations, cls.cgc_bin, cls.cnas, f' {cls.copy_number}')
        df[cls.fusions_gene1] = cls.create_bool(fusions_gene1, cls.cgc_bin, cls.fusions_gene1, f' {cls.rearrangement}')
        df[cls.fusions_gene2] = cls.create_bool(fusions_gene2, cls.cgc_bin, cls.fusions_gene2, f' {cls.rearrangement}')
        df = cls.reset_multi_indexed_dataframe(df.fillna(0), {'level_0': cls.feature, 'level_1': cls.model_id})

        sum_columns = [cls.variants, cls.cnas, cls.fusions_gene1, cls.fusions_gene2]
        df[cls.label] = df.loc[:, sum_columns].sum(axis=1).astype(bool).astype(int)
        return df.pivot_table(index=cls.model_id, columns=cls.feature, values=cls.label)

    @classmethod
    def create_features_list(cls, db):
        genes = db[cls.gene].drop_duplicates()
        series = pd.concat([
            genes + f' {cls.variant}',
            genes + f' {cls.copy_number}',
            genes + f' {cls.rearrangement}',
        ], ignore_index=True)
        return cls.series_to_list(series, drop_duplicates=True)


class Compatibility(Almanac):
    label = 'compatibility'
    label_output = Models.format_label_for_output(label)
    description = 'Inspired by dating algorithms, we weight each molecular feature (or question) based on strength ' \
                  'of the match (e.g. a BRAF deletion only matches BRAF p.V600E by gene). With these relative ' \
                  'weights, we calculate a max score for each sample and compare compare against other cell lines.'

    mapped_alt = 'mapped_alteration'
    score = 'score_contribution'

    base_weight = 10.0
    type_relative = 2.5
    specific_relative = 3.0

    gene_weight = base_weight  # 10
    types_weight = base_weight * type_relative  # 25
    alt_weight = base_weight * type_relative * specific_relative  # 75

    @classmethod
    def append_alt_weights(cls, db, obs_genes, obs_alts, df, dtype):
        for index in db.index:
            tmp_gene = db.loc[index, cls.gene]
            tmp_alt = db.loc[index, cls.mapped_alt]
            if (tmp_gene in obs_genes) & (tmp_alt in obs_alts):
                df.loc[cls.create_slice(genes=tmp_gene, dtype=dtype, alt=tmp_alt), cls.score] = cls.alt_weight
        return df

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        distance_dataframe = cls.calculate_compatibility(input_dtypes)
        distance_dataframe = distance_dataframe.loc[samples_list, samples_list]
        distance_dataframe.index = distance_dataframe.index.tolist()
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        distance_dataframe.to_csv(f'tables/distances/{cls.label_output}.txt', sep='\t')
        return stacked_dataframe

    @classmethod
    def calculate_compatibility(cls, inputs):
        contributions = cls.create_contributions_dataframe(inputs)
        contributions_reset = (contributions
                               .reset_index()
                               .rename(columns={'level_0': cls.model_id,
                                                'level_1': cls.gene,
                                                'level_2': cls.alteration_type,
                                                'level_3': cls.alteration}))
        contributions_reset[cls.string] = contributions_reset[cls.gene] + ' ' + contributions_reset[cls.mapped_alt]
        pivoted = (contributions_reset
                   .pivot_table(values=cls.score, index=cls.model_id, columns=[cls.string], fill_value=0))
        max_values = pivoted.sum(axis=1)
        distance = cls.calculate_compatibility_distance(pivoted, max_values)
        return distance.multiply(distance.T)

    @classmethod
    def calculate_compatibility_distance(cls, df, max_scores):
        list_of_series = []
        for model in df.index:
            features = df.loc[model, :][~df.loc[model, :].eq(0)].index
            series = df.loc[:, features].sum(axis=1).divide(max_scores.loc[model]).sort_values(ascending=False)
            series.name = model
            list_of_series.append(series)
        return (1 - pd.concat(list_of_series, axis=1).T)

    @classmethod
    def create_contributions_dataframe(cls, inputs):
        almanac_dict = cls.prepare_almanac(inputs)
        almanac_genes = almanac_dict[cls.genes]
        almanac_copy_number = almanac_dict[cls.copy_number]
        almanac_fusion = almanac_dict[cls.fusion]
        almanac_variant = almanac_dict[cls.variant]
        almanac_copy_number_genes = cls.extract_genes(almanac_copy_number)
        almanac_fusion_genes = cls.extract_genes(almanac_fusion)
        almanac_variant_genes = cls.extract_genes(almanac_variant)

        df = cls.prepare_observed_alterations(inputs)

        df.loc[(cls.create_slice(genes=almanac_genes)), cls.score] = cls.gene_weight
        df.loc[(cls.create_slice(genes=almanac_copy_number_genes, dtype=cls.copy_number)), cls.score] = cls.types_weight
        df.loc[(cls.create_slice(genes=almanac_fusion_genes, dtype=cls.fusion)), cls.score] = cls.types_weight
        df.loc[(cls.create_slice(genes=almanac_variant_genes, dtype=cls.variant)), cls.score] = cls.types_weight

        obs_genes = df.index.get_level_values(level=1).tolist()
        obs_alts = df.index.get_level_values(level=3).tolist()

        df = cls.append_alt_weights(almanac_copy_number, obs_genes, obs_alts, df, cls.copy_number)
        df = cls.append_alt_weights(almanac_fusion, obs_genes, obs_alts, df, cls.fusion)
        df = cls.append_alt_weights(almanac_variant, obs_genes, obs_alts, df, cls.variant)
        return df

    @classmethod
    def create_slice(cls, models=None, genes=None, dtype=None, alt=None):
        slice_models = slice(None) if models is None else models
        slice_genes = slice(None) if genes is None else genes
        slice_dtype = slice(None) if dtype is None else dtype
        slice_alt = slice(None) if alt is None else alt
        return slice_models, slice_genes, slice_dtype, slice_alt

    @classmethod
    def extract_genes(cls, df):
        return df[cls.gene].drop_duplicates().sort_values().tolist()

    @classmethod
    def prepare_almanac(cls, inputs):
        almanac = cls.import_dbs(inputs)
        db = cls.prepare_almanac_alterations(almanac)
        db[cls.genes] = almanac[cls.genes].tolist()
        return db

    @classmethod
    def prepare_almanac_alterations(cls, db):
        copy_number_alterations = cls.prepare_almanac_copy_number_alterations(db[cls.copy_number])
        fusion_alterations = cls.prepare_almanac_fusions(db[cls.rearrangement])
        variant_alterations = cls.prepare_almanac_variants(db[cls.variant])
        return {
            cls.copy_number: copy_number_alterations,
            cls.fusion: fusion_alterations,
            cls.variant: variant_alterations
        }

    @classmethod
    def prepare_almanac_gene_dtypes(cls, genes):
        return pd.MultiIndex.from_product([genes, [cls.variant, cls.copy_number, cls.rearrangement]])

    @classmethod
    def prepare_almanac_copy_number_alterations(cls, db):
        db[cls.feature_type] = cls.copy_number
        return (db
                .rename(columns={cls.direction: cls.mapped_alt})
                .loc[:, [cls.gene, cls.mapped_alt, cls.feature_type]]
                .dropna()
                .drop_duplicates()
                .reset_index(drop=True))

    @classmethod
    def prepare_almanac_fusions(cls, db):
        db[cls.feature_type] = cls.fusion
        idx_paired = db[cls.gene2].notnull()
        db.loc[idx_paired, cls.mapped_alt] = (db
                                              .loc[idx_paired, [cls.gene1, cls.gene2]]
                                              .apply(lambda x: '--'.join(reversed(sorted(x))), axis=1)
                                              )
        return (pd
                .concat([db.loc[idx_paired, [cls.gene1, cls.mapped_alt, cls.feature_type]]
                        .rename(columns={cls.gene1: cls.gene}),
                         db.loc[idx_paired, [cls.gene2, cls.mapped_alt, cls.feature_type]]
                        .rename(columns={cls.gene2: cls.gene})],
                        ignore_index=True)
                .loc[:, [cls.gene, cls.mapped_alt, cls.feature_type]]
                .sort_values([cls.gene, cls.mapped_alt])
                .drop_duplicates()
                .dropna()
                .reset_index(drop=True))

    @classmethod
    def prepare_almanac_variants(cls, db):
        db[cls.feature_type] = cls.variant
        idx_missense = cls.return_almanac_variants_missense_index(db)
        idx_truncating = cls.return_almanac_variants_truncating_index(db)
        db.loc[idx_missense, cls.mapped_alt] = db.loc[idx_missense, cls.protein_change]
        db.loc[idx_truncating, cls.mapped_alt] = cls.truncating
        return (db
                .loc[:, [cls.gene, cls.mapped_alt, cls.feature_type]]
                .dropna()
                .drop_duplicates()
                .reset_index(drop=True))

    @classmethod
    def prepare_observed_alterations(cls, inputs):
        observed_copy_number_alterations = cls.prepare_observed_copy_number_alterations(inputs[cls.cnas])
        observed_fusions = cls.prepare_observed_fusions(inputs[cls.fusions_gene1], inputs[cls.fusions_gene2])
        observed_variants = cls.prepare_observed_variants(inputs[cls.variants])
        observed_multi_index = (observed_copy_number_alterations.index
                                .union(observed_fusions.index)
                                .union(observed_variants.index))
        return pd.DataFrame(0, index=observed_multi_index, columns=['score_contribution'], dtype=float)

    @classmethod
    def prepare_observed_copy_number_alterations(cls, df):
        df[cls.feature_type] = cls.copy_number
        return (df
                .rename(columns={cls.feature: cls.gene, cls.alteration_type: cls.mapped_alt})
                .loc[:, [cls.model_id, cls.gene, cls.feature_type, cls.mapped_alt]]
                .dropna()
                .drop_duplicates()
                .reset_index(drop=True)
                .set_index([cls.model_id, cls.gene, cls.feature_type, cls.mapped_alt]))

    @classmethod
    def prepare_observed_fusions(cls, df1, df2):
        df = (pd
              .concat([df1.loc[:, [cls.feature, cls.partner, cls.model_id, cls.feature_type]],
                       df2.loc[:, [cls.feature, cls.partner, cls.model_id, cls.feature_type]]],
                       ignore_index=True))
        df[cls.feature_type] = cls.fusion
        idx = df[cls.feature].le(df[cls.partner])
        df.loc[idx, cls.gene1] = df.loc[idx, cls.partner]
        df.loc[idx, cls.gene2] = df.loc[idx, cls.feature]
        df.loc[~idx, cls.gene1] = df.loc[idx, cls.feature]
        df.loc[~idx, cls.gene2] = df.loc[idx, cls.partner]
        df.drop_duplicates([cls.model_id, cls.gene1, cls.gene2], inplace=True)
        df[cls.mapped_alt] = (df
                              .loc[:, [cls.gene1, cls.gene2]]
                              .apply(lambda row: '--'.join(reversed(sorted(row.values.astype(str)))), axis=1))
        df[cls.gene] = df[cls.gene1]
        return (df
                .loc[:, [cls.model_id, cls.gene, cls.feature_type, cls.mapped_alt]]
                .dropna()
                .drop_duplicates()
                .reset_index(drop=True)
                .set_index([cls.model_id, cls.gene, cls.feature_type, cls.mapped_alt]))

    @classmethod
    def prepare_observed_variants(cls, df):
        df[cls.feature_type] = cls.variant
        idx_missense = cls.return_observed_variants_missense_index(df)
        idx_truncating = cls.return_observed_variants_truncating_index(df)
        df.loc[idx_missense, cls.mapped_alt] = df.loc[idx_missense, cls.alteration]
        df.loc[idx_truncating, cls.mapped_alt] = cls.truncating
        return (df
                .rename(columns={cls.feature: cls.gene})
                .loc[:, [cls.model_id, cls.gene, cls.feature_type, cls.mapped_alt]]
                .dropna()
                .drop_duplicates()
                .reset_index(drop=True)
                .set_index([cls.model_id, cls.gene, cls.feature_type, cls.mapped_alt]))

    @classmethod
    def return_almanac_variants_missense_index(cls, db):
        return db[(db[cls.variant_annotation].eq(cls.missense)
                   & ~db[cls.protein_change].eq('')
                   & ~db[cls.protein_change].isin(['p.G12', 'p.G13', 'p.Q61'])
                   )].index

    @classmethod
    def return_almanac_variants_truncating_index(cls, db):
        return db[db[cls.variant_annotation].isin(cls.truncating_types)].index

    @classmethod
    def return_observed_variants_missense_index(cls, df):
        return df[df[cls.alteration_type].eq(cls.missense)].index

    @classmethod
    def return_observed_variants_truncating_index(cls, df):
        return df[df[cls.alteration_type].isin(cls.truncating_types)].index


class NonsynVariantCount(Models):
    label = 'nonsynonymous-variant-count'
    label_output = Models.format_label_for_output(label)
    description = 'We assign neighbors based on the absolute value of the difference of the number of coding somatic ' \
                  'variants. This is a proxy for mutational burden, because we do not have the number of somatic ' \
                  'bases considered when calling variants to use a denominator.'

    n_nonsyn_variants = 'n_nonsyn_variants'

    @classmethod
    def calculate(cls, inputs, samples):
        counts_series = cls.create_counts_series(inputs, samples)
        distance_dataframe = cls.create_difference_dataframe(counts_series)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        counts_series.to_csv('tables/features/{}.txt'.format(cls.label_output), sep='\t')
        return stacked_dataframe

    @staticmethod
    def calculate_pairwise_difference_series(series):
        return abs(series.values - series.values[:, None])

    @classmethod
    def create_counts_series(cls, inputs, samples):
        variants = inputs[cls.variants]
        series = variants[cls.model_id].value_counts().loc[samples].sort_values()
        series.name = cls.label
        return series

    @classmethod
    def create_difference_dataframe(cls, series):
        return pd.DataFrame(cls.calculate_pairwise_difference_series(series), index=series.index, columns=series.index)


class PCA(Models):
    @classmethod
    def calculate_euclidean_distance(cls, df, samples_list):
        euclidean = sklearn.metrics.pairwise.euclidean_distances(df)
        return pd.DataFrame(euclidean, index=samples_list, columns=samples_list)

    @classmethod
    def run_pca(cls, df, samples_list):
        pca = sklearn.decomposition.PCA()
        pca.fit(df)
        transformed = pca.transform(df)
        return pd.DataFrame(transformed, index=samples_list)


class PCAonAlmanac(PCA):
    label = 'pca-almanac-genes'
    label_output = Models.format_label_for_output(label)
    description = 'We run PCA and then nearest neighbors for the vectorization of Almanac genes, with mutants being ' \
                  'without consideration of feature type. For example, there is one feature called "TP53" and both ' \
                  'TP53 nonsense variants and copy number deletions can populate the element.'

    @classmethod
    def create_boolean_table(cls, inputs, samples):
        almanac = Almanac.import_dbs(inputs)
        return AlmanacGenes.create_boolean_table(inputs, samples, almanac)

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        boolean_dataframe = cls.create_boolean_table(input_dtypes, samples_list)
        pca_dataframe = cls.run_pca(boolean_dataframe, boolean_dataframe.index.tolist())
        distance_dataframe = cls.calculate_euclidean_distance(pca_dataframe, boolean_dataframe.index.tolist())
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)
        return stacked_dataframe


class PCAonCGC(PCA):
    label = 'pca-cgc-genes'
    description = 'We run PCA and then nearest neighbors for the vectorization of CGC genes, with mutants being ' \
                  'without consideration of feature type. For example, there is one feature called "TP53" and both ' \
                  'TP53 nonsense variants and copy number deletions can populate the element.'

    @classmethod
    def create_boolean_table(cls, inputs, samples):
        return CGC.create_boolean_table(inputs, samples)

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        boolean_dataframe = cls.create_boolean_table(input_dtypes, samples_list)
        pca_dataframe = cls.run_pca(boolean_dataframe, boolean_dataframe.index.tolist())
        distance_dataframe = cls.calculate_euclidean_distance(pca_dataframe, boolean_dataframe.index.tolist())
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)
        return stacked_dataframe


class RankedSort(Models):
    @classmethod
    def sort(cls, dataframe, models_to_sort, label):
        df = (dataframe
              .reset_index()
              .sort_values(models_to_sort, ascending=True, kind='mergesort')
              .reset_index(drop=True)
              )
        df[label] = (df
                     .reset_index()
                     .groupby('case')['index']
                     .rank(ascending=True, method='first')  # na_option='bottom', method='first')
                     .astype(int)
                     )
        return df.sort_values(['case', label]).set_index(['case', 'comparison']).loc[:, label]


class RankedSortAlmanacEvidenceCGC(RankedSort):
    label = 'multi-pass-sort_fda-cgc'
    label_output = Models.format_label_for_output(label)
    description = 'A weakness of agreement based measure is that there will be tied values. We tie break similarity ' \
                  'based on Molecular Oncology Almanac features associated with FDA evidence by ' \
                  'using similarity based on CGC genes.'

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        almanac_evidence = AlmanacEvidence.calculate(input_dtypes, samples_list)
        cgc = CGCFeatureTypes.calculate(input_dtypes, samples_list)
        combined = cls.sort(
            pd.concat([almanac_evidence, cgc], axis=1),
            [AlmanacEvidence.label, CGCFeatureTypes.label],
            cls.label
        )
        return combined


class RelativeSubstitutionRates(Models):
    label = 'dNdS'
    label_output = label
    description = 'We assign neighbors based on the absolute value of the difference of the number of coding somatic ' \
                  'variants divided by the number of silent variants.'

    n = 'nonsyn_variant_count'
    s = 'syn_variant_count'

    @classmethod
    def calculate(cls, inputs, samples):
        counts_series = cls.create_counts_series(inputs, samples)
        distance_dataframe = cls.create_difference_dataframe(counts_series)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        counts_series.to_csv('tables/features/{}.txt'.format(cls.label_output), sep='\t')
        return stacked_dataframe

    @staticmethod
    def calculate_pairwise_difference_series(series):
        return abs(series.values - series.values[:, None])

    @classmethod
    def create_counts_series(cls, inputs, samples):
        variants = inputs[cls.variants]
        series_ns = variants[cls.model_id].value_counts().loc[samples].sort_values()
        series_ns.name = cls.n

        summary = inputs[cls.summary]
        summary = summary[summary[cls.model_id].isin(samples)]
        series_s = summary.set_index(cls.model_id)[cls.s].sort_values()
        series_s.name = cls.s

        series = series_ns.divide(series_s)
        series.name = cls.label
        return series

    @classmethod
    def create_difference_dataframe(cls, series):
        return pd.DataFrame(cls.calculate_pairwise_difference_series(series), index=series.index, columns=series.index)


class SNFbyEvidenceCGC(Models):
    label = 'snf_fda-cgc-genes'
    label_output = Models.format_label_for_output(label)
    description = 'Rather than perform a two-pass heuristic sort, we use the python implementation of Similarity ' \
                  'Network Fusion (SNF, Bo Wang in the Goldenberg lab(https://www.nature.com/articles/nmeth.2810) by ' \
                  'Ross Markello (https://github.com/rmarkello/snfpy). We fuse networks that describe agreement ' \
                  'based measure of (1) almanac features associated with FDA evidence and (2) any ' \
                  'variant occurring in a Cancer Gene Census gene.'

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        almanac = AlmanacEvidence.import_dbs(input_dtypes)
        almanac_subset = AlmanacEvidence.subset_almanac_by_evidence(almanac, ['FDA-Approved'])
        boolean_dataframe_1 = AlmanacFeatures.create_boolean_table(input_dtypes, samples_list, almanac_subset)
        boolean_dataframe_2 = CGC.create_boolean_table(input_dtypes, samples_list)

        data = [boolean_dataframe_1.loc[samples_list, :],
                boolean_dataframe_2.loc[samples_list, :],
                ]
        affinity_networks = snf.make_affinity(data,
                                              metric='jaccard',
                                              normalize=False,
                                              K=20, mu=0.5)
        fused_network = snf.snf(affinity_networks, K=20)
        fused_dataframe = pd.DataFrame(fused_network, index=samples_list, columns=samples_list)
        distance_dataframe = pd.DataFrame(1, index=samples_list, columns=samples_list).subtract(fused_dataframe)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        # boolean_dataframe.to_csv('tables/distances/{}.boolean.txt'.format(cls.label_output), sep='\t')
        return stacked_dataframe


class SNFTypesAlmanac(Models):
    label = 'snf_almanac'
    description = 'Rather than collapse all data types into a single similarity matrix (e.g. with columns such as ' \
                  'CDKN2A somatic variant, CDKN2A copy number alteration), we use the python implementation of ' \
                  'Similarity Network Fusion (SNF, Bo Wang in the Goldenberg lab ' \
                  '(https://www.nature.com/articles/nmeth.2810) by Ross Markello(https://github.com/rmarkello/snfpy).' \
                  'We fuse networks that describe agreement based of variants in almanac genes in ' \
                  '(1) somatic variants, (2) copy number alterations, and (3) rearrangements.'

    @classmethod
    def create_boolean_table(cls, inputs, dtype, samples):
        almanac = Almanac.import_dbs(inputs)
        features = Almanac.generate_gene_features(almanac)

        if dtype == cls.fusions:
            columns = [cls.fusions_gene1, cls.fusions_gene2]
        else:
            columns = [dtype]

        index = pd.MultiIndex.from_product([features, samples])
        df = pd.DataFrame(columns=columns, index=index)
        for column in columns:
            df[column] = AlmanacGenes.create_bool(inputs[column], Almanac.feature_match_1, column)
        df = cls.reset_multi_indexed_dataframe(df.fillna(0), {'level_0': cls.feature, 'level_1': cls.model_id})

        sum_columns = columns
        df[cls.label] = df.loc[:, sum_columns].sum(axis=1).astype(bool).astype(int)
        return df.pivot_table(index=cls.model_id, columns=cls.feature, values=cls.label)

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        boolean_dataframe_variants = cls.create_boolean_table(input_dtypes, cls.variants, samples_list)
        boolean_dataframe_copy_numbers = cls.create_boolean_table(input_dtypes, cls.cnas, samples_list)
        boolean_dataframe_fusions = cls.create_boolean_table(input_dtypes, cls.fusions, samples_list)

        data = [boolean_dataframe_variants.loc[samples_list, :],
                boolean_dataframe_copy_numbers.loc[samples_list, :],
                boolean_dataframe_fusions.loc[samples_list, :],
                ]
        affinity_networks = snf.make_affinity(data,
                                              metric='jaccard',
                                              normalize=False,
                                              K=20, mu=0.5)
        fused_network = snf.snf(affinity_networks, K=20)
        fused_dataframe = pd.DataFrame(fused_network, index=samples_list, columns=samples_list)
        distance_dataframe = pd.DataFrame(1, index=samples_list, columns=samples_list).subtract(fused_dataframe)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        # boolean_dataframe.to_csv('tables/distances/{}.boolean.txt'.format(cls.label_output), sep='\t')
        return stacked_dataframe


class SNFTypesCGC(Models):
    label = 'snf_cgc'
    label_output = Models.format_label_for_output(label)
    description = 'Rather than collapse all data types into a single similarity matrix (e.g. with columns such as ' \
                  'CDKN2A somatic variant, CDKN2A copy number alteration), we use the python implementation of ' \
                  'Similarity Network Fusion (SNF, Bo Wang in the Goldenberg lab ' \
                  '(https://www.nature.com/articles/nmeth.2810) by Ross Markello(https://github.com/rmarkello/snfpy).' \
                  'We fuse networks that describe agreement based of variants in CGC genes in ' \
                  '(1) somatic variants, (2) copy number alterations, and (3) rearrangements.'

    @classmethod
    def create_boolean_table(cls, inputs, dtype, samples):
        cgc = inputs[CGC.input_label]
        features = CGC.create_features_list(cgc)

        if dtype == cls.fusions:
            columns = [cls.fusions_gene1, cls.fusions_gene2]
        else:
            columns = [dtype]

        index = pd.MultiIndex.from_product([features, samples])
        df = pd.DataFrame(columns=columns, index=index)
        for column in columns:
            df[column] = CGC.create_bool(inputs[column], CGC.cgc_bin, column)
        df = cls.reset_multi_indexed_dataframe(df.fillna(0), {'level_0': cls.feature, 'level_1': cls.model_id})

        sum_columns = columns
        df[cls.label] = df.loc[:, sum_columns].sum(axis=1).astype(bool).astype(int)
        return df.pivot_table(index=cls.model_id, columns=cls.feature, values=cls.label)

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        boolean_dataframe_variants = cls.create_boolean_table(input_dtypes, cls.variants, samples_list)
        boolean_dataframe_copy_numbers = cls.create_boolean_table(input_dtypes, cls.cnas, samples_list)
        boolean_dataframe_fusions = cls.create_boolean_table(input_dtypes, cls.fusions, samples_list)

        data = [boolean_dataframe_variants.loc[samples_list, :],
                boolean_dataframe_copy_numbers.loc[samples_list, :],
                boolean_dataframe_fusions.loc[samples_list, :],
                ]
        affinity_networks = snf.make_affinity(data,
                                              metric='jaccard',
                                              normalize=False,
                                              K=20, mu=0.5)
        fused_network = snf.snf(affinity_networks, K=20)
        fused_dataframe = pd.DataFrame(fused_network, index=samples_list, columns=samples_list)
        distance_dataframe = pd.DataFrame(1, index=samples_list, columns=samples_list).subtract(fused_dataframe)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        # boolean_dataframe.to_csv('tables/distances/{}.boolean.txt'.format(cls.label_output), sep='\t')
        return stacked_dataframe


class SNFTypesCGCwithEvidence(Models):
    label = 'snf_fda-cgc'
    label_output = Models.format_label_for_output(label)
    description = 'We perform similarity network fusion (SNF, Bo Wang in the Goldenberg lab ' \
                  '(https://www.nature.com/articles/nmeth.2810)) using the python implementation by Ross Markello ' \
                  '(https://github.com/rmarkello/snfpy) to fuse networks that contain: ' \
                  '(1) CGC genes that contain a somatic variant, ' \
                  '(2) CGC genes that contain a copy number alteration, ' \
                  '(3) CGC genes that contain a rearrangement, ' \
                  '(4) Almanac features associated with FDA evidence.'

    @classmethod
    def create_boolean_table(cls, inputs, dtype, samples):
        cgc = inputs[CGC.input_label]
        features = CGC.create_features_list(cgc)

        if dtype == cls.fusions:
            columns = [cls.fusions_gene1, cls.fusions_gene2]
        else:
            columns = [dtype]

        index = pd.MultiIndex.from_product([features, samples])
        df = pd.DataFrame(columns=columns, index=index)
        for column in columns:
            df[column] = CGC.create_bool(inputs[column], CGC.cgc_bin, column)
        df = cls.reset_multi_indexed_dataframe(df.fillna(0), {'level_0': cls.feature, 'level_1': cls.model_id})

        sum_columns = columns
        df[cls.label] = df.loc[:, sum_columns].sum(axis=1).astype(bool).astype(int)
        return df.pivot_table(index=cls.model_id, columns=cls.feature, values=cls.label)

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        boolean_dataframe_variants = cls.create_boolean_table(input_dtypes, cls.variants, samples_list)
        boolean_dataframe_copy_numbers = cls.create_boolean_table(input_dtypes, cls.cnas, samples_list)
        boolean_dataframe_fusions = cls.create_boolean_table(input_dtypes, cls.fusions, samples_list)

        almanac = AlmanacEvidence.import_dbs(input_dtypes)
        almanac_subset = AlmanacEvidence.subset_almanac_by_evidence(almanac, ['FDA-Approved'])
        boolean_dataframe_1 = AlmanacFeatures.create_boolean_table(input_dtypes, samples_list, almanac_subset)

        data = [boolean_dataframe_variants.loc[samples_list, :],
                boolean_dataframe_copy_numbers.loc[samples_list, :],
                boolean_dataframe_fusions.loc[samples_list, :],
                boolean_dataframe_1.loc[samples_list, :],
                ]
        affinity_networks = snf.make_affinity(data,
                                              metric='jaccard',
                                              normalize=False,
                                              K=20, mu=0.5)
        fused_network = snf.snf(affinity_networks, K=20)
        fused_dataframe = pd.DataFrame(fused_network, index=samples_list, columns=samples_list)
        distance_dataframe = pd.DataFrame(1, index=samples_list, columns=samples_list).subtract(fused_dataframe)
        stacked_dataframe = cls.stack_distances(distance_dataframe, cls.label)

        # boolean_dataframe.to_csv('tables/distances/{}.boolean.txt'.format(cls.label_output), sep='\t')
        return stacked_dataframe


class Tree(Models):
    label = 'somatic-tree'
    label_output = Models.format_label_for_output(label)
    description = 'This is somewhat inspired by CELLector by Hanna Najgebauer and Euan Stronach ' \
                  '(doi: 10.1016/j.cels.2020.04.007). One issue with agreement based measure is that each feature is ' \
                  'weighted the same. CELLector has a sorted list of genes/variants based on cancer type and will ' \
                  'report similar cell lines based on mutant / wild type status of each gene. While not exactly the ' \
                  'same, we use the annotations from various datasources appended to variants by the Molecular ' \
                  'Oncology Almanac to create a priority list for variants (hotspots ranked the highest, etc.). ' \
                  'For each case sample, we consider the genes which are observed to be mutated and preserve the ' \
                  'order that they would appear in the somatic.scored.txt output of Molecular Oncology Almanac. All ' \
                  'other samples are then sorted by their mutant / wild type status of these genes.'

    @classmethod
    def calculate(cls, input_dtypes, samples_list):
        sorted_dataframe = cls.create_sorted_table(input_dtypes, samples_list)
        sorted_dataframe = sorted_dataframe[sorted_dataframe['model_id'].isin(samples_list)]
        features = sorted_dataframe['feature'].drop_duplicates().sort_values().tolist()
        boolean_dataframe = cls.create_boolean_table(sorted_dataframe, samples_list, features)
        stacked_dataframe = cls.calculate_tree_distance(samples_list, sorted_dataframe, boolean_dataframe)
        return stacked_dataframe

    @classmethod
    def calculate_tree_distance(cls, samples, dataframe_features, dataframe_bool):
        dataframe_features.reset_index(inplace=True)
        distances = []
        empty_counter = 0
        for sample in samples:
            tmp_features = dataframe_features.copy(deep=True)
            tmp_booleans = dataframe_bool.copy(deep=True)
            features = tmp_features[tmp_features['model_id'].eq(sample)]['feature'].tolist()
            if not features:
                empty_counter += 1
                print(sample, empty_counter)
                continue
            relative_booleans = tmp_booleans.loc[:, features]
            relative_booleans = (relative_booleans
                                 .sort_values(by=relative_booleans.columns.tolist(), ascending=False)
                                 .reset_index()
                                 .reset_index()
                                 .rename(columns={'model_id': cls.comparison, 'index': cls.label})
                                 )
            relative_booleans[cls.case] = sample
            distances.append(relative_booleans.loc[:, [cls.case, cls.comparison, cls.label]])
        series = pd.concat(distances, ignore_index=True).set_index([cls.case, cls.comparison])
        series.name = cls.label
        return series

    @classmethod
    def create_boolean_table(cls, dataframe, samples, features):
        dataframe.set_index(['feature', 'model_id'], inplace=True)
        dataframe['value'] = 1

        index = pd.MultiIndex.from_product([features, samples])
        df = pd.DataFrame(0, columns=['value'], index=index)
        df.loc[dataframe.index, 'value'] = dataframe.loc[dataframe.index, 'value']
        return (df
                .reset_index().rename(columns={'level_0': 'feature', 'level_1': 'model_id'})
                .pivot_table(index='model_id', columns='feature', values='value')
                )

    @classmethod
    def create_sorted_table(cls, inputs, samples):
        variants = inputs[cls.variants]
        copy_number_alterations = inputs[cls.cnas]
        fusions_gene1 = inputs[cls.fusions_gene1]
        fusions_gene2 = inputs[cls.fusions_gene2]

        sort_columns = [  # 'feature_match_4', 'feature_match_3', 'feature_match_2', 'feature_match_1', 'evidence_map',
            'cancerhotspots_bin', 'cancerhotspots3D_bin', 'cgc_bin', 'gsea_pathways_bin',
            'gsea_modules_bin', 'cosmic_bin']

        variants['feature'] = variants['feature'] + ' Mut'
        copy_number_alterations['feature'] = copy_number_alterations['feature'] + ' CN'
        fusions_gene1['feature'] = fusions_gene1['feature'] + ' Fusion'
        fusions_gene2['feature'] = fusions_gene2['feature'] + ' Fusion'

        df = pd.concat([
            variants.loc[:, ['model_id', 'feature'] + sort_columns],  # + ' mut'
            copy_number_alterations.loc[:, ['model_id', 'feature'] + sort_columns],  # + '  CN'
            fusions_gene1.loc[:, ['model_id', 'feature'] + sort_columns],  # + '  Fusion'
            fusions_gene2.loc[:, ['model_id', 'feature'] + sort_columns]  # + '  Fusion'
        ])

        df = df.sort_values(sort_columns, ascending=False).drop_duplicates(subset=['model_id', 'feature'], keep='first')
        return df[df['model_id'].isin(samples)].reset_index(drop=True)
