import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


from metrics import Metrics


class Plots:
    @classmethod
    def holding(cls):
        return ''


class AveragePrecision(Plots):
    @classmethod
    def plot(cls, models_dictionary, models_list, outname='models'):
        list_avg_ps = []
        for model in models_list:
            series = pd.Series(models_dictionary[model][Metrics.avg_precision], name=model)
            list_avg_ps.append(series)
        avg_ps = pd.concat(list_avg_ps, axis=1)

        mean_sorted = avg_ps.mean(axis=0).sort_values(ascending=False).index.tolist()
        avg_ps = avg_ps.loc[:, mean_sorted]

        fig, ax = plt.subplots(figsize=(15, 7.5))
        ax = sns.boxplot(data=avg_ps)
        ax = sns.stripplot(data=avg_ps, color="black", jitter=0.2, size=2.5)

        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

        ax.set_ylabel('Average Precision', fontsize=16)

        # plt.title("Boxplot with jitter", loc="left")
        plt.xticks(rotation=45, fontsize=14, ha='right')
        plt.savefig(f'img/{outname}.avg_precision.png', bbox_inches='tight', dpi=300)


class AveragePrecisionK(Plots):
    @classmethod
    def plot(cls, models_dictionary, models_list, outname='models', ylim=0.20):
        columns = ['model', 'k', 'ap@k']
        list_ = []
        for model in models_list:
            apk = models_dictionary[model]['ap@k']
            for key in [1, 2, 3, 4, 5]:
                # for key in apk.keys():
                list_.append([model, key, apk[key]])
        data = pd.DataFrame(list_, columns=columns)

        fig, ax = plt.subplots(figsize=(7.5, 10))
        ax = sns.pointplot(x="k", y="ap@k", hue="model", data=data, ax=ax, fontsize=12)

        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

        ax.set_ylim([0, ylim])
        ax.set_ylabel('Average Precision @ K', fontsize=16)
        ax.set_xlabel('Rank (k)', fontsize=16)
        ax.set_title('Average Precision @ K performance', fontsize=18)
        ax.legend(fontsize=14, frameon=False)
        plt.savefig(f'img/{outname}.avg_precision_at_k.png', bbox_inches='tight', dpi=300)
