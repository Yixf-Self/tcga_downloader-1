import matplotlib.pyplot as plt
import numpy as np
import random
import seaborn as sns
import scipy


tumor = "BLCA"
junction_file = "stefano/ctcf_biostring_lucilla.tsv"
max_distance = 100000
number_null_simulations = 1
list_genes = []
list_junctions = []
list_mutations = []


list_chroms = ["chr" + str(x) for x in range(1,23)]

class GeneAnnotation:
    def __init__(self, l):
        s = l.strip().split()
        self.name = s[0]
        self.chrom = s[1]
        self.position = int(s[2])


class Junction:
    def __init__(self, l):
        s = l.strip().split()
        self.chrom = "chr" + s[0]
        self.start = int(s[1])
        self.stop = int(s[2])
        self.center = self.start + 9
        self.left_junction = None
        self.right_junction = None
        self.null_score = []
        self.left_genes = []
        self.right_genes = []

    def add_left_genes(self, list_left):
        self.left_genes = sorted(list_left, key = lambda x : abs(x.position-self.center))

    def add_right_genes(self, list_right):
        self.right_genes = sorted(list_right, key = lambda x : abs(x.position-self.center))

    def add_left_junction(self, left_position):
        self.left_junction = left_position

    def add_right_junction(self, right_position):
        self.right_junction = right_position

    def set_score(self, s):
        self.score = s

    def add_null_score(self, s):
        self.null_score.append(s)

class Mutation:
    def __init__(self, l):
        s = l.strip().split()
        self.chrom = s[0]
        self.position = int(s[1])


def map_junction_genes(ljunctions, lgenes):
    for c in list_chroms:
        print(c)
        fj = [x for x in ljunctions if x.chrom == c]
        fg = [x for x in lgenes if x.chrom == c]
        for j in fj:
            right_genes = list(filter(lambda x: x.position > j.center and
                                                abs(x.position - j.center) < max_distance,
                                      fg))
            left_genes = list(filter(lambda x: x.position < j.center and
                                                abs(x.position - j.center) < max_distance,
                                      fg))
            j.add_right_genes(right_genes)
            j.add_left_genes(left_genes)

def map_left_and_right_junction(ljunctions):
    for c in list_chroms:
        print(c)
        fj = [x for x in ljunctions if x.chrom == c]
        for j in fj:
            left_js = list(filter(lambda x: x.center < j.center,
                                  fj))
            right_js = list(filter(lambda x: x.center > j.center,
                                   fj))
            if len(left_js) > 0:
                j.add_left_junction(list(sorted(left_js, key = lambda x : abs(x.center-j.center)))[0])

            if len(right_js) > 0:
                j.add_right_junction(list(sorted(right_js, key = lambda x : abs(x.center-j.center)))[0])


def parse_expr_data(l):
    s = l.strip().split()
    name = s[0]
    exprs = list(map(float, s[3:]))
    return (name, exprs)

#junctions with at least a gene on the left and a gene on the right
def point_D_filter(j):
    if len(j.left_genes) > 0 and len(j.right_genes) > 0:
        return True
    else:
        return False


class Point_E_ratio:
    def __init__(self):
        self.msg = "mean of ratio"

    def score(self, j, normal_expr, tumor_expr):
        list_ratio = []
        for x in j.left_genes + j.right_genes:
            mt = np.mean(tumor_expr[x.name])+1.0
            mn = np.mean(normal_expr[x.name])+1.0
            list_ratio.append(mt/mn)
        j.set_score(np.mean(list_ratio))

    def null_score(self, j, normal_expr, tumor_expr):

        random_ratio = []
        for x in random.sample(normal_expr.keys(), len(j.left_genes + j.right_genes)):
            random_ratio.append((np.mean(tumor_expr[x])+1.0) /
                                (np.mean(normal_expr[x])+1.0))
        j.add_null_score(np.mean(random_ratio))

def get_outliers_slope(distribution):

    for i in range(len(distribution)-5):
        if ((distribution[i+5].score - distribution[i].score) / 5.0) > 5.0:
            return i+4
    return len(distribution)

def check_mutation_juntion_intersection(list_mutations, list_junctions):

    count_mutations = 0
    count_intersecting = 0
    for c in list_chroms:
        fm = [x for x in list_mutations if x.chrom == c]
        fj = [x for x in list_junctions if x.chrom == c]
        count_mutations += len(fm)
        count_intersecting += len([m for m in list_mutations if len([j for j in list_junctions if j.start <= m.position <= j.stop]) > 0])

    print(str(count_mutations)," mutations; ", str(count_intersecting), " overlapping")


def main():

    normal_expr = dict([parse_expr_data(x) for x in open("TCGA-" + tumor + "-normal/TCGA-" + tumor + "-expression-normal.matrix")])
    tumor_expr = dict([parse_expr_data(x) for x in open("TCGA-" + tumor + "-normalAssociated/TCGA-" + tumor + "-expression-normalAssociated.matrix")])
    print("normal genes: ", len(normal_expr.keys()))
    print("tumor genes: ", len(tumor_expr.keys()))


    list_genes = [GeneAnnotation(l) for l in open("TCGA-" + tumor + "-normal/TCGA-" + tumor + "-expression-normal.matrix")]
    print("considered genes: ", len(list_genes))

    list_junctions = [Junction(l) for l in open(junction_file)]
    print("considered junctions: ", len(list_junctions))

    map_junction_genes(list_junctions, list_genes)

    map_left_and_right_junction(list_junctions)

    list_junctions = list(filter(point_D_filter, list_junctions))
    print("Junctions after filtering: ", len(list_junctions))

    list_mutations = [Mutation(x) for x in open("stefano/somatic_mutations.txt")]
    check_mutation_juntion_intersection(list_mutations, list_junctions)

    fig1, ax = plt.subplots()
    ax.boxplot([[len(x.left_genes) for x in list_junctions],
                [len(x.right_genes) for x in list_junctions]])
    ax.set_title("Junctions (" + str(len(list_junctions)) + ") - Genes (" + str(len(list_genes)) + ")" + ", threshold = " + str(max_distance))
    #ax.set_yscale('log')
    ax.set_xticklabels(["Left", "Right"])
    ax.set_ylabel("Genes Count")
    plt.savefig("s-"+ tumor +"left-right-distribution.png")


    method_E = Point_E_ratio()

    for j in list_junctions:
        method_E.score(j, normal_expr, tumor_expr)

    for i in range(number_null_simulations):
        print(str(i))
        for j in list_junctions:
            method_E.null_score(j, normal_expr, tumor_expr)

    score_values = sorted([x.score for x in list_junctions])

    sorted_junctions = list(sorted(list_junctions, key = lambda x: x.score))
    #top_outliers_position = get_outliers_slope(sorted_junctions)
    top_outliers_position = int(0.95*len(sorted_junctions))
    bottom_outliers_position = int(0.05*len(sorted_junctions))

    fig1, ax = plt.subplots()
    ax.plot(
        range(len(score_values)),
        score_values,
        c='r'
    )
    for i in range(number_null_simulations):
        print(i)
        null_score_values = sorted([j.null_score[i] for j in list_junctions])
        ax.plot(
            range(len(null_score_values)),
            null_score_values,
            c='b'
        )
    ax.axvline(top_outliers_position)
    ax.axvline(bottom_outliers_position)
    ax.set_title(tumor + ": sorted score " + method_E.msg)
    ax.set_yscale('log')
    ax.set_ylabel("Score")
    plt.savefig("s-"+ tumor +"score-sorted.png")

    count_bottom_mutations = []
    count_junctions_mutations = []
    count_outliers_mutations = []
    for c in list_chroms:
        fb = [x for x in sorted_junctions[0:bottom_outliers_position] if x.chrom == c]
        fj = [x for x in sorted_junctions[bottom_outliers_position:top_outliers_position] if x.chrom == c]
        fo = [x for x in sorted_junctions[top_outliers_position:] if x.chrom == c]
        fm = [x for x in list_mutations if x.chrom == c]

        for j in fb:
            count_bottom_mutations.append(len([m for m in fm if j.start < m.position < j.stop]))

        for j in fj:
            count_junctions_mutations.append(len([m for m in fm if j.start < m.position < j.stop]))

        for j in fo:
            count_outliers_mutations.append(len([m for m in fm if j.start < m.position < j.stop]))

    print("Mutation per junction: ",  top_outliers_position-bottom_outliers_position, np.mean(count_junctions_mutations), np.std(count_junctions_mutations))
    print("Mutation per top outlier: ", len(sorted_junctions) - top_outliers_position, np.mean(count_outliers_mutations), np.std(count_outliers_mutations))
    print("Mutation per low outlier: ", bottom_outliers_position, np.mean(count_bottom_mutations), np.std(count_bottom_mutations))


    mutations_count = []
    for i in range(100):
        print(i)
        count_mut = []
        samp = random.sample(sorted_junctions, bottom_outliers_position)
        for c in list_chroms:
            fj = [x for x in samp if x.chrom == c]
            fm = [x for x in list_mutations if x.chrom == c]
            for j in fj:
                count_mut.append(len([m for m in fm if j.start < m.position < j.stop]))
        mutations_count.append(np.mean(count_mut))


    fig1, ax = plt.subplots()
    sns.distplot(mutations_count, ax=ax)
    plt.savefig("s-"+ tumor +"-distribution.png")

    print(np.mean(mutations_count), np.std(mutations_count))


if __name__ == "__main__":

    main()