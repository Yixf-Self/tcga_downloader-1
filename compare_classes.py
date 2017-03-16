import matplotlib.pyplot as plt
import numpy as np
import random
import seaborn as sns
import scipy



tumor = "PRAD"


class IntermediateTuple():
    def __init__(self, line):
        s = line.strip().split()
        self.name = s[0]
        self.junction_start = int(s[1])
        self.junction_end = int(s[2])
        self.tss = int(s[3])
        self.distance = int(s[4])
        self.avg_normal = float(s[5])
        self.avg_tumor = float(s[6])
        self.std_normal = float(s[7])
        self.std_tumor = float(s[8])
        self.delta = float(s[9])
        self.ratio = float(s[10])
        self.mutation_count = int(s[11])
        self.ratio1 = float(s[12])
        self.cluster = s[13]


def main():

    data = [IntermediateTuple(l) for l in open("eirini/" + tumor + ".tsv")]
    print("all genes: ", len(data))

    data_mutation = list(filter(lambda x : x.cluster == "mutations", data))
    print("mutations: ", len(data_mutation))

    data_nomutation = list(filter(lambda x : x.cluster == "no_mutations", data))
    print("no_mutations: ", len(data_nomutation))

    mut_normal_fpkm = [x.avg_normal for x in data_mutation]
    mut_tumor_fpkm = [x.avg_tumor for x in data_mutation]
    nomut_normal_fpkm = [x.avg_normal for x in data_nomutation]
    nomut_tumor_fpkm = [x.avg_tumor for x in data_nomutation]

    mut_ratio = [x.ratio1 for x in data_mutation]
    nomut_ratio = [x.ratio1 for x in data_nomutation]

    fig1, ax = plt.subplots()
    ax.boxplot([mut_normal_fpkm,
                mut_tumor_fpkm])
    ax.set_title(tumor + ": Genes close to mutated junction (n=" + str(len(data_mutation)) + ")")
    ax.set_yscale('log')
    ax.set_xticklabels(["normal", "tumoral"])
    ax.set_ylabel("Mean (FPKM-UQ)")
    plt.savefig("c-mutations-mean-expression-" + tumor +'.png')

    fig1, ax = plt.subplots()
    ax.boxplot([nomut_normal_fpkm,
                nomut_tumor_fpkm])
    ax.set_title(tumor + ": Genes close to non mutated junction (n=" + str(len(data_nomutation)) + ")")
    ax.set_yscale('log')
    ax.set_xticklabels(["normal", "tumoral"])
    ax.set_ylabel("Mean (FPKM-UQ)")
    plt.savefig("c-no_mutations-mean-expression-" + tumor +'.png')

    pval_mean_diff = scipy.stats.ttest_ind(np.array(mut_ratio), np.array(nomut_ratio))[1]
    print("mean difference", pval_mean_diff)
    fig1, ax = plt.subplots()
    ax.boxplot([mut_ratio,
                nomut_ratio])
    ax.set_title(tumor + ": Expression ratio (p="+str(pval_mean_diff)+")")
    ax.set_yscale('log')
    ax.set_xticklabels(["mutations", "no_mutations"])
    ax.set_ylabel("Ratio of Mean (FPKM-UQ)")
    plt.savefig("c-expression-ratio-" + tumor +'.png')




if __name__ == "__main__":

    main()