import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, kstest


def draw_comp(
    n=50,
    rangepos=[0, 1],
    rug=200,
    comp="uniform",
    nrep=1000,
    chrom_start=0,
    chrom_end=1e6,
):

    stat_store = []

    for idx in range(nrep):
        pos = np.random.randint(*rangepos, size=n)

        xs = np.linspace(chrom_start, chrom_end, rug)
        density = gaussian_kde(pos)

        Y = density(xs)

        kt = kstest((pos - rangepos[0]) / rangepos[1], comp)

        stat_store.append([*kt, density.factor, np.mean(Y), np.std(Y)])

    return np.array(stat_store)


def draw_comp_ks(
    n=50,
    rangepos=[0, 1],
    rug=200,
    comp="uniform",
    nrep=1000,
    geno_bornes=[0, 1000],
    ksteps=5e5,
    kwindow=1e6,
    chrom_start=0,
    chrom_end=1e6,
):

    stat_store = []
    count_store = []

    for idx in range(nrep):
        pos = np.random.randint(*rangepos, size=n)

        xs = np.linspace(chrom_start, chrom_end, rug)
        density = gaussian_kde(pos)

        Y = density(xs)

        kt = kstest((pos - rangepos[0]) / rangepos[1], comp)

        stat_store.append([*kt, density.factor, np.mean(Y), np.std(Y)])

        ##
        ##
        pos = np.random.randint(*geno_bornes, size=n)
        pos = sorted(pos)
        #
        td = np.array(pos, dtype=int)
        td = pd.DataFrame(td, columns=["POS"])

        Windows, Out = geno_Lwind_split(
            td, geno_bornes=[chrom_start, chrom_end], Steps=ksteps, window_size=kwindow
        )
        #

        kp = [len(set(x)) for x in Windows.values()]

        count_store.extend(kp)

    return np.array(stat_store), count_store


def geno_Lwind_split(summary, geno_bornes=[], Steps=25e3, window_size=5e4):
    """
    split genotype array into windows by length, steps.
    assumes genotype has a single chrom.
    """

    POS = np.array(summary.POS, dtype=int)
    if not geno_bornes:
        ## assume that chromosome does not end at last INV position;
        ## without genome sizes, best bet is that INVs are uniformily distributed.
        geno_bornes = [0, (max(POS) + min(POS))]

    window_starts = np.array(
        np.arange(geno_bornes[0], geno_bornes[1], Steps), dtype=int
    )

    Windows = {x: [] for x in window_starts}
    Out = {z: z + window_size for z in window_starts}

    current_winds = []

    for idx in range(len(POS)):
        posh = int(POS[idx])

        if len(window_starts):
            if window_starts[0] <= posh:
                current_winds.append(window_starts[0])

                window_starts = window_starts[1:]

            d = 1 - int(len(window_starts) > 0)
            ids = 0
            while d == 0:
                if Out[window_starts[ids]] < posh:
                    ids += 1
                else:
                    current_winds.append(window_starts[0])

                    window_starts = window_starts[1:]
                    d += 1

            window_starts = window_starts[ids:]

        current_rm = []

        for windx in current_winds:
            if posh > Out[windx]:
                if idx == len(POS) - 1:
                    current_rm.append(windx)
                else:
                    if POS[idx + 1] > Out[windx]:
                        current_rm.append(windx)

                continue

            Windows[windx].append(idx)

        current_winds = [x for x in current_winds if x not in current_rm]

    return Windows, Out
