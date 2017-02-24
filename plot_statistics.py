import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sp
import sys

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

    ig1, ax = plt.subplots()
    ax.scatter(list(map(np.log, [np.mean(x[3]) for x in tumoral_data])),
               list(map(np.log, [np.median(x[3]) for x in tumoral_data])), s=1)
    ax.plot([0,20], [0,20], c='r')
    ax.set_title(sys.argv[1])
    ax.set_xlabel("Log normal mean expression")
    ax.set_ylabel("Log normal median expression")
    plt.savefig("tumoral-mean-median-" + sys.argv[1] +'.png')

    pvals = sorted(map(lambda z : (z[0][0], sp.ttest_ind(z[0][3],z[1][3])[1]),
                       zip(normal_data, tumoral_data)),
                   key = lambda x :x[1])

    print(list(pvals)[0:10])

    print([x for x in pvals if x[0] == 'ENSG00000136997'])


def usage():
    print("script.py TUMOR_TAG")

if __name__ == "__main__":

    main()
