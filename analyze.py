import pdb
from scipy.stats import poisson
import numpy as np
import matplotlib.pyplot as plt


def calculate_mutation_rates(df):
    df = df[df['n_anomalies'] < 1].copy()

    df = df.drop(columns = ['n_anomalies'])

    # In the math explainer we showed that we can just count the total mutations over the total length
    df['total_len'] = df['sample_count'] * df['ref_length']

    df['mutation_rate'] = df['n_mutations'] / df['total_len']

    total_mutation_rate = df['n_mutations'].sum() / df['total_len'].sum()

    df['expected_mutations'] = df['total_len'] * total_mutation_rate

    df['p_value'] = df.apply(lambda x: poisson.cdf(x['n_mutations'], x['total_len'] * total_mutation_rate), axis=1)

    df['neg_log_p'] = - np.log10(df['p_value'])

    df['unique_ratio'] = df['n_sequences'] / df['sample_count']

    return df


def make_summary_plots(df):
    y = df['neg_log_p'].sort_values()

    x = np.log10(np.arange(len(y)) + 1)

    plt.scatter(x, y)
    plt.xlabel('-log ( expected pvals )')
    plt.ylabel('-log ( actual pvals )')

    plt.show()
    plt.savefig('qq_plot')
    plt.close()


    data = df['p_value']
    plt.hist(data, 100, alpha=0.5, color="black")
    plt.xlabel('p-value')
    plt.ylabel('count')
    plt.show()
    plt.close()


    data = df[df['p_value'] == 0]['mutation_rate']
    plt.hist(data, 100, alpha=0.5, color="black")
    plt.xlabel('Mutations per base pair')
    plt.show()
    plt.close()





if __name__ == '__main__':

    import pickle
    df = pickle.load(open('df', 'rb'))

    print("--------------------------")
    print(len(df[df['n_anomalies'] == 0]))
    print(len(df[df['n_anomalies'] < 10]))
    print(df.sort_values('n_anomalies'))

    df = calculate_mutation_rates(df)

    #print("-----------------------------")
    zeros = df['p_value'] == 0
    ones = df['p_value'] == 1
    df_significant = df[zeros | ones]
    pdb.set_trace()
    #print(df_significant.sort_values('mutation_rate').iloc[-40:])


    make_summary_plots(df)