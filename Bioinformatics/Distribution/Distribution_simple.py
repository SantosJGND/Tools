import os

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, kstest

from utilities import draw_comp_ks, geno_Lwind_split

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", type=str, default="filtered_HQ_INV.bed")

    parser.add_argument("--kwindow", type=int, default=2e5)

    parser.add_argument("--ksteps", type=int, default=1e5)

    parser.add_argument("--alpha", type=float, default=1e-2)

    parser.add_argument("--cmet", type=str, default="fdr_bh")

    parser.add_argument("--centro", type=str, default="")

    args = parser.parse_args()

    kwindow = args.kwindow
    ksteps = args.ksteps
    alpha = args.alpha
    corr_method = args.cmet

    file = args.input

    ###
    if args.centro:
        centromeres = pd.read_csv(args.centro, header=None, sep="\t")
        centromeres.columns = ["ID", "chrom", "IN", "OUT", "L"]

    ###
    invbed = pd.read_csv(file, sep="\t")
    invbed = invbed.rename(
        columns={"CHROM": "chrom", "clustStart": "IN", "clustEnd": "OUT"}
    )
    print(invbed.head())

    # invbed.columns = ["chrom", "IN", "OUT", "L"]

    ## replace names in bed column
    # invbed.chrom= invbed.chrom.str.split('.',n= 1, expand= True)[1]

    invbed.head()

    out_dir = (
        os.path.splitext(os.path.basename(file))[0]
        + "_{}_{}".format(alpha, corr_method)
        + "/"
    )
    print(out_dir)

    os.makedirs(out_dir, exist_ok=True)

    #######################
    #######################
    pdout = []
    kt_out = []
    bed_filtered = []
    for chrom in invbed.chrom.unique():

        bedsel = invbed.loc[invbed.chrom == chrom].reset_index(drop=True)
        bedsel = bedsel.sort_values(by="IN")

        bedsel["IN"] = pd.to_numeric(bedsel["IN"])
        bedsel["OUT"] = pd.to_numeric(bedsel["OUT"])

        print(bedsel.head(30))
        bed_filtered.append(bedsel)

        chrom_start = 0
        chrom_end = bedsel.OUT.max() + bedsel.IN.min()

        #################################################
        ################################################# Kolmogrov-Smirnov test
        ## using means.
        vals = (bedsel.IN + bedsel.OUT) / 2
        props = list(vals / chrom_end)
        vals = list(vals)

        kt = kstest(props, "uniform")

        ### power of Kolmogorov-smirnov test
        ###

        stat_store, count_store = draw_comp_ks(
            n=bedsel.shape[0],
            rangepos=[chrom_start, chrom_end],
            comp="uniform",
            nrep=10000,
            geno_bornes=[chrom_start, chrom_end],
            ksteps=ksteps,
            kwindow=kwindow,
            chrom_start=chrom_start,
            chrom_end=chrom_end,
        )

        cthist, ctidx = np.histogram(
            count_store, bins=len(set(count_store)), density=True
        )
        ct_cdf = [
            1 - (sum(cthist[:x]) / sum(cthist)) for x in range(1, len(cthist) + 1)
        ]

        #
        stat_store = pd.DataFrame(
            stat_store, columns=["KS", "PVAL", "DF", "DMEAN", "DSTD"]
        )
        stat_store.head()

        pvalktest = kstest(list(stat_store.PVAL), "uniform")

        power_here = stat_store[stat_store.PVAL < alpha].shape[0] / stat_store.shape[0]
        report = "power (KS | H0): {} for alpha = {}".format(
            stat_store[stat_store.PVAL < alpha].shape[0] / stat_store.shape[0], alpha
        )

        kt_out.append([chrom, bedsel.shape[0], *list(kt), power_here])

    kt_out = np.array(kt_out)
    kt_out = pd.DataFrame(kt_out, columns=["CHROM", "NINV", "KS", "PVAL", "FDR"])

    kt_out.to_csv(out_dir + "KStest.txt", index=False, sep="\t")
