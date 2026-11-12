#!/usr/bin/env python3
"""Find reciprocal crossovers: pairs of strains whose ancestry switches in the same SNP
interval in opposite directions (one A->B, the other B->A), the signature of the two
reciprocal products of one crossover, i.e. of meiotic siblings (or their descendants).

Input, either
  --anc  <prefix>.anc.tsv.gz from haploswarm.py / split_generations.py (true ancestry per SNP), or
  --vcf  a haploid or homozygous-diploid VCF plus --parent-a / --parent-b: sample names whose
         alleles define ancestry (the real parents, or proxies for the two lineages). Only sites
         where the two differ are used; a strain's ancestry at a site is the parent it matches.

For each group of strains (generation, from --samples; otherwise all strains together):
  * breakpoints: every ancestry switch between two adjacent informative SNPs
  * reciprocal events: SNP intervals where at least one strain switches A->B and another B->A
  * expected_by_chance: reciprocal events expected if the distinct A->B and B->A breakpoint
    intervals of each chromosome were placed independently of each other (so high SNP
    density makes chance coincidences rare)
  * shared_same_direction: intervals where two or more strains switch the same way (shared
    descent, e.g. clones or descendants of one recombinant)
Strains with identical ancestry (clones) are counted once (--keep-clones to count each).

Usage:
  python reciprocal_crossovers.py --anc bench3/tetrad.anc.tsv.gz --samples bench3/tetrad.samples.tsv --out bench3/tetrad_recip
  python reciprocal_crossovers.py --vcf isolates.vcf.gz --parent-a P1 --parent-b P2 --out recip

Writes <out>.summary.tsv (one row per group) and <out>.events.tsv (one row per reciprocal event).
"""
import argparse, collections, gzip, sys
import numpy as np


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_anc(path):
    f = opener(path)
    hdr = f.readline().rstrip('\n').split('\t')
    if hdr[:2] != ['chrom', 'pos']:
        sys.exit('Error: ancestry table must start with chrom and pos columns (made by haploswarm.py)')
    chrom, pos, rows = [], [], []
    for l in f:
        c = l.rstrip('\n').split('\t')
        chrom.append(c[0]); pos.append(int(c[1]))
        rows.append([1 if x == '1' else 0 if x == '0' else -1 for x in c[2:]])
    return hdr[2:], np.array(chrom), np.array(pos), np.array(rows, dtype=np.int8)   # rows: SNPs x samples


def read_vcf(path, pa, pb):
    ids, chrom, pos, rows = None, [], [], []
    for l in opener(path):
        if l.startswith('##'):
            continue
        c = l.rstrip('\n').split('\t')
        if l.startswith('#'):
            ids = c[9:]
            for p in (pa, pb):
                if p not in ids:
                    sys.exit(f'Error: {p} is not a sample in {path}')
            ia, ib = ids.index(pa), ids.index(pb)
            continue
        gt = [x.split(':')[0].replace('|', '/').split('/')[0] for x in c[9:]]
        a, b = gt[ia], gt[ib]
        if a in ('.', '') or b in ('.', '') or a == b:
            continue                                   # not informative between the parents
        rows.append([1 if g == a else 0 if g == b else -1 for g in gt])
        chrom.append(c[0]); pos.append(int(c[1]))
    return ids, np.array(chrom), np.array(pos), np.array(rows, dtype=np.int8)


def analyse(cols, names, chrom, pos, A, keep_clones):
    """Breakpoints and reciprocal events among the strains in cols."""
    X = A[:, cols]
    if not keep_clones:                                # one representative per identical ancestry
        _, first = np.unique(X.T, axis=0, return_index=True)
        first = np.sort(first)
        X = X[:, first]; cols = [cols[i] for i in first]
    n_strains = X.shape[1]
    events, n_bp, n_dist, n_shared, expected = [], 0, 0, 0, 0.0
    for ch in dict.fromkeys(chrom):                    # chromosome order as in the file
        idx = np.where(chrom == ch)[0]
        Y = X[idx]
        # fill missing calls from the left so they don't create breakpoints
        for j in range(Y.shape[1]):
            col = Y[:, j]
            if not (col < 0).any():
                continue
            for k in range(1, len(col)):
                if col[k] < 0:
                    col[k] = col[k - 1]
        d = np.diff(Y.astype(np.int16), axis=0)       # interval k: between SNP k and k+1
        ab = (d == -1)                                # A (1) -> B (0)
        ba = (d == 1)
        n_bp += int(ab.sum() + ba.sum())
        m = len(idx) - 1
        cab, cba = ab.sum(1), ba.sum(1)
        if m > 0:   # chance: the distinct A->B and B->A intervals placed independently among the m intervals
            expected += (cab > 0).sum() * (cba > 0).sum() / m
        n_dist += int(((cab + cba) > 0).sum())
        n_shared += int(((cab >= 2) | (cba >= 2)).sum())
        for k in np.where((cab > 0) & (cba > 0))[0]:
            events.append((ch, int(pos[idx[k]]), int(pos[idx[k + 1]]),
                           ','.join(names[cols[j]] for j in np.where(ab[k])[0]),
                           ','.join(names[cols[j]] for j in np.where(ba[k])[0])))
    return n_strains, n_bp, n_dist, events, expected, n_shared


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--anc', help='ancestry table (.anc.tsv.gz) from haploswarm.py or split_generations.py')
    ap.add_argument('--vcf', help='VCF of real or simulated isolates (use with --parent-a and --parent-b)')
    ap.add_argument('--parent-a'); ap.add_argument('--parent-b')
    ap.add_argument('--samples', help='samples.tsv: analyse each group (generation) separately; references are skipped')
    ap.add_argument('--keep-clones', action='store_true', help='count strains with identical ancestry separately')
    ap.add_argument('--out', required=True, help='output prefix')
    a = ap.parse_args()

    if a.anc:
        names, chrom, pos, A = read_anc(a.anc)
    elif a.vcf and a.parent_a and a.parent_b:
        names, chrom, pos, A = read_vcf(a.vcf, a.parent_a, a.parent_b)
    else:
        sys.exit('Error: give --anc, or --vcf with --parent-a and --parent-b')
    print(f'{len(names)} samples, {len(pos)} informative SNPs on {len(set(chrom))} chromosomes', flush=True)

    skip = {a.parent_a, a.parent_b}
    groups = collections.OrderedDict()
    if a.samples:
        info = {l.split('\t')[1]: l.rstrip('\n').split('\t') for l in list(open(a.samples))[1:]}
        for i, n in enumerate(names):
            r = info.get(n)
            if r is None or r[3] == 'ref' or n in skip:
                continue
            groups.setdefault(r[2], []).append(i)
        gen = lambda g: int(g.split('_g')[1]) if '_g' in g else 10**9
        groups = collections.OrderedDict(sorted(groups.items(), key=lambda kv: (gen(kv[0]), kv[0])))
    else:
        groups['all'] = [i for i, n in enumerate(names) if n not in skip]

    with open(a.out + '.summary.tsv', 'w') as fs, open(a.out + '.events.tsv', 'w') as fe:
        fs.write('group\tstrains\tdistinct_genotypes\tbreakpoints\tbreakpoints_per_genotype\tdistinct_breakpoint_intervals'
                 '\treciprocal_events\texpected_by_chance\tobserved_over_expected\tshared_same_direction\n')
        fe.write('group\tchrom\tleft_snp\tright_snp\tstrains_A_to_B\tstrains_B_to_A\n')
        print(f"{'group':12s} {'strains':>7s} {'genotypes':>9s} {'breakpts':>8s} {'reciprocal':>10s} {'by chance':>9s} {'shared':>7s}")
        for g, cols in groups.items():
            n_all = len(cols)
            ns, nbp, nd, ev, ex, sh = analyse(cols, names, chrom, pos, A, a.keep_clones)
            oe = len(ev) / ex if ex > 0 else float('nan')
            fs.write(f'{g}\t{n_all}\t{ns}\t{nbp}\t{nbp / max(ns, 1):.2f}\t{nd}\t{len(ev)}\t{ex:.2f}\t{oe:.2f}\t{sh}\n')
            for e in ev:
                fe.write(g + '\t' + '\t'.join(map(str, e)) + '\n')
            print(f'{g:12s} {n_all:7d} {ns:9d} {nbp:8d} {len(ev):10d} {ex:9.1f} {sh:7d}', flush=True)
    print(f'-> {a.out}.summary.tsv, {a.out}.events.tsv')


if __name__ == '__main__':
    main()
