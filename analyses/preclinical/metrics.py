import numpy as np
import pandas as pd


class Metrics:
    k = 'k'
    relevant = 'relevant'
    avg_precision = 'average_precision'
    mean_avg_precision = 'mean_average_precision'
    true_positives = 'tps'
    avg_precision_at_k = 'ap@k'
    precision_at_k = 'p@k'
    recall_at_k = 'r@k'
    true_positives_at_k = 'tps@k'

    calculated = 'calculated'
    description = 'description'

    case = 'case'
    comparison = 'comparison'
    n_shared = 'n_shared'

    @classmethod
    def calculate_average_precision(cls, samples, df):
        average_precision = pd.Series(index=samples, name=cls.avg_precision)
        true_positives = df[df[cls.relevant].astype(bool)]
        for name, group in true_positives.reset_index().groupby(cls.case):
            average_precision.loc[name] = group[cls.precision_at_k].mean()
        return average_precision

    @classmethod
    def calculate_average_precision_at_k(cls, df, k):
        return df[df[cls.k].eq(k)][cls.precision_at_k].mean()

    @classmethod
    def calculate_precision_recall_at_k(cls, df, true_positives):
        df[cls.precision_at_k] = 0
        df[cls.recall_at_k] = 0
        df[cls.true_positives_at_k] = 0

        for name, group in df.groupby(cls.case):
            idx = group.index
            df.loc[idx, cls.true_positives_at_k] = df.loc[idx, cls.relevant].astype(int).cumsum()
            df.loc[idx, cls.precision_at_k] = df.loc[idx, cls.true_positives_at_k].divide(df.loc[idx, cls.k])
            df.loc[idx, cls.recall_at_k] = df.loc[idx, cls.true_positives_at_k].divide(true_positives.loc[name])
        return df

    @classmethod
    def create_true_positives_series(cls, samples, dataframe):
        series = pd.Series(0, index=samples, name=cls.true_positives)
        for label, group in dataframe.groupby(cls.case)[cls.n_shared]:
            series.loc[label] = group.astype(bool).sum()
        return series

    @classmethod
    def evaluate_model(cls, df, model, true_positives_per_sample):
        df[cls.k] = cls.rank_comparisons(df[model], model)
        ranked = cls.sort_by_rank(df)
        return cls.calculate_precision_recall_at_k(ranked, true_positives_per_sample)

    @classmethod
    def evaluate_models(cls, samples, df, models, descriptions):
        df = cls.remove_identical_comparisons(df)
        df = df.reset_index()
        df = df[df['case'].isin(samples) & df['comparison'].isin(samples)]
        df = df.set_index(['case', 'comparison'])
        df[cls.relevant] = df[cls.n_shared].astype(bool).astype(int)
        tps_series = cls.create_true_positives_series(samples, df)

        models_dictionary = {}
        for model in models:
            description = descriptions[model]
            calculated = cls.evaluate_model(df, model, tps_series)
            model_avg_precision = cls.calculate_average_precision(samples, calculated)
            model_mean_avg_precision = model_avg_precision.mean()

            models_dictionary[model] = {}
            models_dictionary[model][cls.description] = description
            models_dictionary[model][cls.calculated] = calculated
            models_dictionary[model][cls.avg_precision] = model_avg_precision
            models_dictionary[model][cls.mean_avg_precision] = model_mean_avg_precision
            models_dictionary[model][cls.avg_precision_at_k] = {}
            for k in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 25, 50, 100]:
                models_dictionary[model][cls.avg_precision_at_k][k] = (
                    cls.calculate_average_precision_at_k(calculated, k))
            cls.print_model(models_dictionary, model)
        return models_dictionary

    @classmethod
    def print_model(cls, models_dictionary, model_label):
        print(model_label,
              models_dictionary[model_label][cls.avg_precision_at_k][1],
              models_dictionary[model_label][cls.avg_precision_at_k][2],
              models_dictionary[model_label][cls.avg_precision_at_k][3],
              models_dictionary[model_label][cls.mean_avg_precision])

    @classmethod
    def rank_comparisons(cls, series, model):
        reset = series.reset_index()
        reset[cls.k] = reset.groupby(cls.case)[model].rank(ascending=True, method='max')  # first
        return reset.set_index([cls.case, cls.comparison]).loc[:, cls.k]

    @classmethod
    def remove_identical_comparisons(cls, dataframe):
        df = dataframe.reset_index()
        df = df[df[cls.case].ne(df[cls.comparison])]
        return df.set_index([cls.case, cls.comparison])

    @classmethod
    def sort_by_rank(cls, dataframe):
        return dataframe.reset_index().sort_values([cls.case, cls.k]).set_index([cls.case, cls.comparison])


class Random(Metrics):
    @classmethod
    def evaluate_random_model(cls, samples, df, number_of_iterations):
        df = cls.remove_identical_comparisons(df)
        df[cls.relevant] = df[cls.n_shared].astype(bool).astype(int)
        tps_series = cls.create_true_positives_series(samples, df)

        results = {}
        for seed in range(0, number_of_iterations):
            dataframe = df.copy(deep=True)
            dataframe = dataframe.reset_index()

            rng = np.random.default_rng(seed=seed)
            dataframe['model'] = rng.permutation(dataframe.index)
            dataframe = dataframe.set_index(['case', 'comparison'])
            calculated = cls.evaluate_model(dataframe, 'model', tps_series)
            model_avg_precision = cls.calculate_average_precision(samples, calculated)
            model_mean_avg_precision = model_avg_precision.mean()
            results[seed] = model_mean_avg_precision

            del dataframe
        return results

    @classmethod
    def evaluate_random_models(cls, samples, df, seeds_dictionary):
        df = cls.remove_identical_comparisons(df)
        df[cls.relevant] = df[cls.n_shared].astype(bool).astype(int)
        tps_series = cls.create_true_positives_series(samples, df)

        models_dictionary = {}
        for label, seed in seeds_dictionary.items():
            dataframe = df.copy(deep=True)
            dataframe = dataframe.reset_index()
            rng = np.random.default_rng(seed=seed)
            dataframe[label] = rng.permutation(dataframe.index)
            dataframe = dataframe.set_index(['case', 'comparison'])

            calculated = cls.evaluate_model(dataframe, label, tps_series)
            model_avg_precision = cls.calculate_average_precision(samples, calculated)
            model_mean_avg_precision = model_avg_precision.mean()

            models_dictionary[label] = {}
            models_dictionary[label][cls.calculated] = calculated
            models_dictionary[label][cls.avg_precision] = model_avg_precision
            models_dictionary[label][cls.mean_avg_precision] = model_mean_avg_precision
            models_dictionary[label][cls.avg_precision_at_k] = {}
            for k in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 25, 50, 100]:
                models_dictionary[label][cls.avg_precision_at_k][k] = (
                    cls.calculate_average_precision_at_k(calculated, k))
            cls.print_model(models_dictionary, label)
        return models_dictionary
