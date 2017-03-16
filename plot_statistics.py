import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sp
import sys
import random
import seaborn as sns

def parser_matrix(line):
    s = line.strip().split("\t")

    name = s[0]
    chrom = s[1]
    start = s[2]
    expr = list(map(float, s[3:]))

    return (name, chrom, start, expr)

def main():
    if len(sys.argv) < 2:
        usage()
        sys.exit(2)

    normal_file = "TCGA-" + sys.argv[1] + "-normal/TCGA-" + sys.argv[1] + \
                  "-expression-normal.matrix"
    tumoral_file = "TCGA-" + sys.argv[1] + "-normalAssociated/TCGA-" + \
                   sys.argv[1] + "-expression-normalAssociated.matrix"

    normal_data = [parser_matrix(x) for x in open(normal_file, 'r').readlines()]
    tumoral_data = [parser_matrix(x) for x in open(tumoral_file, 'r').readlines()]

    for i in range(0, len(tumoral_data)):
        assert tumoral_data[i][0] == normal_data[i][0], "genes are different!"


    fig1, ax = plt.subplots()
    ax.boxplot([[np.mean(x[3]) for x in normal_data],
                [np.mean(x[3]) for x in tumoral_data]])
    ax.set_title(sys.argv[1])
    ax.set_yscale('log')
    ax.set_xticklabels(["normal", "tumoral"])
    ax.set_ylabel("Expression mean distribution")
    plt.savefig("expression-mean-distribution-" + sys.argv[1] +'.png')

    fig1, ax = plt.subplots()
    ax.scatter(list(map(np.log, [np.mean(x[3]) for x in normal_data])),
               list(map(np.log, [np.mean(x[3]) for x in tumoral_data])), s=1)
    ax.plot([0,20], [0,20], c='r')
    ax.set_title(sys.argv[1])
    ax.set_xlabel("Log normal mean expression")
    ax.set_ylabel("Log tumoral mean expression")
    plt.savefig("differential-mean-expression-" + sys.argv[1] +'.png')

    fig1, ax = plt.subplots()
    ax.boxplot([[np.median(x[3]) for x in normal_data],
                [np.median(x[3]) for x in tumoral_data]])
    ax.set_title(sys.argv[1])
    ax.set_yscale('log')
    ax.set_xticklabels(["normal", "tumoral"])
    ax.set_ylabel("Expression median distribution")
    plt.savefig("expression-median-distribution-" + sys.argv[1] +'.png')

    fig1, ax = plt.subplots()
    ax.scatter(list(map(np.log, [np.median(x[3]) for x in normal_data])),
               list(map(np.log, [np.median(x[3]) for x in tumoral_data])), s=1)
    ax.plot([0,20], [0,20], c='r')
    ax.set_title(sys.argv[1])
    ax.set_xlabel("Log normal median expression")
    ax.set_ylabel("Log tumoral median expression")
    plt.savefig("differential-median-expression-" + sys.argv[1] +'.png')

    ig1, ax = plt.subplots()
    ax.scatter(list(map(np.log, [np.mean(x[3]) for x in normal_data])),
               list(map(np.log, [np.median(x[3]) for x in normal_data])), s=1)
    ax.plot([0,20], [0,20], c='r')
    ax.set_title(sys.argv[1])
    ax.set_xlabel("Log normal mean expression")
    ax.set_ylabel("Log normal median expression")
    plt.savefig("normal-mean-median-" + sys.argv[1] +'.png')

    fig1, ax = plt.subplots()
    ax.scatter(list(map(np.log, [np.mean(x[3]) for x in tumoral_data])),
               list(map(np.log, [np.median(x[3]) for x in tumoral_data])), s=1)
    ax.plot([0,20], [0,20], c='r')
    ax.set_title(sys.argv[1])
    ax.set_xlabel("Log normal mean expression")
    ax.set_ylabel("Log normal median expression")
    plt.savefig("tumoral-mean-median-" + sys.argv[1] +'.png')

    #boxplot of the correlation ratio
    fig1, ax = plt.subplots()
    ratio = list(map(lambda x: (np.mean(x[1][3])+1.0) / (np.mean(x[0][3])+1.0),
                     zip(normal_data, tumoral_data)))
    ax.boxplot(ratio)
    ax.set_yscale('log')
    plt.savefig("ratio-distribution-" + sys.argv[1] +'.png')

    fig1, ax = plt.subplots()
    delta = list(map(lambda x: (np.mean(x[1][3])+1.0) - (np.mean(x[0][3])+1.0),
                     zip(normal_data, tumoral_data)))
    ax.boxplot(delta)
    plt.savefig("delta-distribution-" + sys.argv[1] +'.png')

    ratio_mean_distribution = []
    for x in range(1000):
        ratio_mean_distribution.append(np.mean(random.sample(ratio, 200)))

    #boxplot of the ratio mean
    fig1, ax = plt.subplots()
    ax.boxplot([ratio, ratio_mean_distribution])
    ax.set_yscale('log')
    plt.savefig("ratio-mean-boxplot-" + sys.argv[1] +'.png')

    fig1, ax = plt.subplots(figsize=(10,10))
    sns.distplot([np.log(x) for x in ratio], ax=ax)
    sns.distplot([np.log(x) for x in ratio_mean_distribution], ax=ax)
    plt.savefig("ratio-mean-distplot-" + sys.argv[1] +'.png')

    normal_patient_sum = []
    tumor_patient_sum = []
    for i in range(len(normal_data[0][3])):
        s = 0.0
        for g in normal_data:
            s += g[3][i]
        normal_patient_sum.append(s)

    for i in range(len(tumoral_data[0][3])):
        s = 0.0
        for g in tumoral_data:
            s += g[3][i]
        tumor_patient_sum.append(s)

    fig1, ax = plt.subplots(figsize=(10,10))
    ax.boxplot([normal_patient_sum, tumor_patient_sum])
    plt.savefig("patient-wise-sum-" + sys.argv[1] +'.png')


def usage():
    print("script.py TUMOR_TAG")

if __name__ == "__main__":

    main()
