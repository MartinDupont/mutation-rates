import os
from scipy.stats import poisson
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from common import THIS_DIRECTORY


def calculate_mutation_rates(df, anomaly_threshold=0, h_null_mutation_rate=None):
    df = df[df['n_anomalies'] <= anomaly_threshold].copy()

    df = df.drop(columns=['n_anomalies'])

    # In the math explainer we showed that we can just count the total mutations over the total length
    df['total_len'] = df['sample_count'] * df['ref_length']

    df['mutation_rate'] = df['n_mutations'] / (df['total_len'] * df['ref_length'])

    if h_null_mutation_rate is None:
        h_null_mutation_rate = df['n_mutations'].sum() / df['total_len'].sum()

    df['expected_mutations'] = df['total_len'] * h_null_mutation_rate

    df['p_value'] = df.apply(lambda x: poisson.cdf(x['n_mutations'], x['total_len'] * h_null_mutation_rate), axis=1)

    df['predicted_essential'] = df['p_value'].apply(lambda x: x < 0.05)

    df['unique_ratio'] = df['n_sequences'] / df['sample_count']

    return df


def join_on_essential(df):
    df_essential = pd.read_csv('essential_genes.csv', sep='\t')[['Gene Name', 'Product']]

    df_matched = df.merge(df_essential, left_on=['identifier'], right_on=['Gene Name'], how='left')

    df_matched['essential'] = ~df_matched['Product'].isnull()

    return df_matched

def confusion_matrix(df_matched):
    tp = sum(df_matched['predicted_essential'] & df_matched['essential'])
    fp = sum(df_matched['predicted_essential'] & (~df_matched['essential']))
    fn = sum((~df_matched['predicted_essential']) & df_matched['essential'])
    tn = sum((~df_matched['predicted_essential']) & (~df_matched['essential']))

    matrix = np.array([[tp, fp], [fn, tn]])

    confusion = pd.DataFrame(matrix, columns = ['essential', 'non-essential'], index=['predicted_essential', 'predicted_non-essential'])
    print(confusion)

    print()
    confusion['essential'] = confusion['essential'] / confusion['essential'].sum()
    confusion['non-essential'] = confusion['non-essential'] / confusion['non-essential'].sum()

    print(confusion)



def make_summary_plots(df, df_matched):

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5))

    data = df['predicted_essential'].apply(lambda x: 0 if x else 1)
    counts = data.value_counts().sort_index()

    ax1.bar(counts.index, counts.values, color=['b', 'r'], alpha=0.5)
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels([True, False])

    ax1.set_xlabel('p-value < 0.05')
    ax1.set_ylabel('count')
    ax1.set_title("p-values all genes")

    data = df_matched['predicted_essential'].apply(lambda x: 0 if x else 1)
    counts = data.value_counts().sort_index()
    ax2.bar(counts.index, counts.values, color=['b', 'r'], alpha=0.5)
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels([True, False])

    ax2.set_xlabel('p-value < 0.05')
    ax2.set_ylabel('count')
    ax2.set_title("p-values essential genes")

    plt.tight_layout()
    plt.show()


    data = df['mutation_rate']
    plt.hist(data, 40, alpha=0.5, color="b", log=True)
    plt.xlabel('Mutations per base pair')
    plt.ylabel('Log count of genes')
    plt.title('Distribution of mutation frequency across genes')
    plt.show()
    plt.close()


if __name__ == '__main__':

    outfile = file_path = os.path.join(THIS_DIRECTORY, "df_mutations.pkl")
    df = pickle.load(open(outfile, 'rb'))

    # For the given data, this results in about 5% of pvals being less than 5%./
    low_mutation_rate = 0.00005

    df = calculate_mutation_rates(df, h_null_mutation_rate=low_mutation_rate)

    df_matched = join_on_essential(df)
    confusion_matrix(df_matched)

    df_essential = df_matched[df_matched['essential']]
    make_summary_plots(df, df_essential)
