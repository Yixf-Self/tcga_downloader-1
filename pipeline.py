import numpy as np
import random
import os
import pylatex
from scipy import stats


junction_file = "additional_data/ctcf_biostring_lucilla.tsv"
number_null_simulations = 1
list_genes = []
list_junctions = []
list_mutations = []
#tumors = ["PRAD", "BLCA", "BRCA", "LIHC", "LUAD", "KIRC", "KIRP", "UCEC", "THCA", "HNSC", "LUSC"]
tumors = ["PRAD", "BLCA", "LIHC", "LUAD", "KIRC", "KIRP", "UCEC", "THCA", "HNSC", "LUSC"]
max_distances = [10000, 50000, 100000,200000]
thresholds = [0.01, 0.05, 0.1]


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
        self.tag = "MeanOfRatio"

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

class Point_E_maxLR:
    def __init__(self):
        self.msg = "max left/right mean ratio"
        self.tag = "MaxLeftRightMeanRatio"

    def score(self, j, normal_expr, tumor_expr):
        list_ratio_left = []
        list_ratio_right = []
        for x in j.left_genes:
            mt = np.mean(tumor_expr[x.name])+1.0
            mn = np.mean(normal_expr[x.name])+1.0
            list_ratio_left.append(mt/mn)
        left_ratio = np.mean(list_ratio_left)
        for x in j.right_genes:
            mt = np.mean(tumor_expr[x.name])+1.0
            mn = np.mean(normal_expr[x.name])+1.0
            list_ratio_right.append(mt/mn)
        right_ratio = np.mean(list_ratio_right)

        if left_ratio > right_ratio:
            j.set_score(left_ratio)
        else:
            j.set_score(right_ratio)

    def null_score(self, j, normal_expr, tumor_expr):

        j.add_null_score(-1.0)

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

def get_pval(mu, sigma, v):
    if sigma == 0:
        return 0
    else:
        P = stats.norm.cdf((v-mu)/ sigma)
        if P > 0.5:
            return 1.0 - P
        else:
            return P

def do_analysis(doc):

    normal_expr = dict([parse_expr_data(x) for x in open("TCGA-" + tumor + "-normal/TCGA-" + tumor + "-expression-normal.matrix")])
    tumor_expr = dict([parse_expr_data(x) for x in open("TCGA-" + tumor + "-normalAssociated/TCGA-" + tumor + "-expression-normalAssociated.matrix")])

    list_genes = [GeneAnnotation(l) for l in open("TCGA-" + tumor + "-normal/TCGA-" + tumor + "-expression-normal.matrix")]
    print("considered genes: ", len(list_genes))

    list_junctions = [Junction(l) for l in open(junction_file)]
    print("considered junctions: ", len(list_junctions))

    map_junction_genes(list_junctions, list_genes)

    map_left_and_right_junction(list_junctions)

    list_junctions = list(filter(point_D_filter, list_junctions))
    print("Junctions after filtering: ", len(list_junctions))

    list_mutations = [Mutation(x) for x in open("additional_data/somatic_mutations.txt")]
    check_mutation_juntion_intersection(list_mutations, list_junctions)

    for j in list_junctions:
        e_method.score(j, normal_expr, tumor_expr)

    #for i in range(number_null_simulations):
    #    for j in list_junctions:
    #        e_method.null_score(j, normal_expr, tumor_expr)

    sorted_junctions = list(sorted(list_junctions, key = lambda x: x.score))

    with doc.create(pylatex.Tabular("l|r|r|r|r|r",
                             row_height=1.5)) as data_table:
        data_table.add_row(["class",
                            "thresh.",
                            "num.",
                            "mean",
                            "std",
                            "pval"],
                           mapper=pylatex.utils.bold,
                           color="lightgray")
        data_table.add_hline()
        for px in thresholds:
            #top_outliers_position = get_outliers_slope(sorted_junctions)
            top_outliers_position = int((1-px)*len(sorted_junctions))
            bottom_outliers_position = int(px*len(sorted_junctions))
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

            null_mean_distribution = []
            for i in range(250):
                samp = random.sample(count_junctions_mutations, bottom_outliers_position)
                null_mean_distribution.append(np.mean(samp))

            null_mean = np.mean(null_mean_distribution)
            null_sigma = np.mean(null_mean_distribution)

            data_table.add_row(("center",
                                str(px),
                                str(top_outliers_position-bottom_outliers_position),
                                '%.5f' % np.mean(count_junctions_mutations),
                                '%.5f' % np.std(count_junctions_mutations),
                                get_pval(null_mean, null_sigma, np.mean(count_junctions_mutations))
                                ))
            if np.mean(count_bottom_mutations) > null_mean:
                bottom_name = "bottom UP"
            else:
                bottom_name = "bottom DOWN"
            data_table.add_row((bottom_name,
                                str(px),
                                str(bottom_outliers_position),
                                '%.5f' % np.mean(count_bottom_mutations),
                                '%.5f' % np.std(count_bottom_mutations),
                                get_pval(null_mean, null_sigma, np.mean(count_bottom_mutations))
                                ))
            if np.mean(count_outliers_mutations) > null_mean:
                top_name = "top UP"
            else:
                top_name = "top DOWN"
            data_table.add_row((top_name,
                                str(px),
                                str(len(sorted_junctions) - top_outliers_position),
                                '%.5f' % np.mean(count_outliers_mutations),
                                '%.5f' % np.std(count_outliers_mutations),
                                get_pval(null_mean, null_sigma, np.mean(count_outliers_mutations))
                                ))
            data_table.add_hline()


e_methods = [Point_E_ratio(), Point_E_maxLR()]

def main():

    if not os.path.exists("report"):
        os.mkdir("report")

    global tumor
    global max_distance
    global e_method

    for t in tumors:
        print(t)
        doc = pylatex.Document('report/results'+t)
        doc.append(pylatex.NewPage())
        with doc.create(pylatex.Section(t)):
            for m in e_methods:
                for d in max_distances:
                    tumor = t
                    max_distance = d
                    e_method = m
                    with doc.create(pylatex.Subsection("Score = " + m.tag + ", Max dist. = " + str(d))):
                        do_analysis(doc)

        doc.generate_pdf(clean_tex=False)



if __name__ == "__main__":

    main()
