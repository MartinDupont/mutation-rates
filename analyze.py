import pdb
from scipy.stats import poisson
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


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

    df['neg_log_p'] = - np.log10(df['p_value'])

    df['unique_ratio'] = df['n_sequences'] / df['sample_count']

    df['gene_id'] = df['identifier'].str.split(':').str[0]

    return df


def join_on_essential(df):
    df_essential = pd.read_csv('essential_genes.csv', sep='\t')[['Gene Name', 'Product']]

    df_matched = df.merge(df_essential, left_on=['gene_id'], right_on=['Gene Name'], how='inner')

    return df_matched


def make_summary_plots(df, df_matched):

    # Create a 1x2 grid layout for two plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    data = df['p_value']
    median = df['p_value'].median()
    ax1.hist(data, 100, alpha=0.5, color="b")
    ax1.axvline(median, label='median', color='red')
    ax1.set_xlabel('p-value')
    ax1.set_ylabel('count')
    ax1.set_title("p-value distribution across genes")
    ax1.legend()

    data = df_matched['p_value']
    median = df_matched['p_value'].median()
    ax2.hist(data, 100, alpha=0.5, color="b")
    ax2.axvline(median, label='median', color='red')
    ax2.set_xlabel('p-value')
    ax2.set_ylabel('count')
    ax2.set_title("p-value distribution across essential genes")
    ax2.legend()

    # Adjust layout to avoid overlap
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
    import pickle

    df = pickle.load(open('df', 'rb'))

    low_mutation_rate = 0.0001
    df_low = calculate_mutation_rates(df, h_null_mutation_rate=low_mutation_rate)
    zeros = df_low['p_value'] < 1e-8
    df_significantly_low = df_low[zeros]

    print('------------------------------------------')
    print('Genes with anomalously low mutation rates')
    print('------------------------------------------')
    print(df_significantly_low)

    high_mutation_rate = 0.005
    df_high = calculate_mutation_rates(df, h_null_mutation_rate=high_mutation_rate)
    ones = df_high['p_value'] > 1 - 1e-8
    df_significantly_high = df_high[ones]

    print('------------------------------------------')
    print('Genes with anomalously high mutation rates')
    print('------------------------------------------')
    print(df_significantly_high)

    df_middle = calculate_mutation_rates(df)

    df_matched = join_on_essential(df_middle)

    make_summary_plots(df_middle, df_matched)
