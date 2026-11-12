#!/usr/bin/env python3
"""Positive control for ADMIXTURE: data built exactly as its model assumes, from the same parents.

From a haploswarm.py run made with --include-parents, for every sampled generation this writes
two input sets, each with the same clonal source populations (copies of parent A and parent B,
zero SNP differences from the parents):

  <outdir>/swarm_sources/g<gg>/   the real swarm strains of that generation + the clonal sources
  <outdir>/model/g<gg>/           model-conforming strains + the clonal sources: one per swarm
                                  strain, with the same true A fraction q, but every SNP drawn
                                  independently (A allele with probability q), so there is no
                                  linkage and no relatedness, as ADMIXTURE assumes

Both use the same SNPs (every site where the parents differ). The only difference between the
two sets is the swarm's linkage and relatedness.

Files per generation (same layout as split_generations.py, so run_admixture_by_gen.sh,
plot_admixture.R and compare_generations.R work unchanged):
  <tag>.dip.vcf.gz, .vcf.gz, .geno, .pop (A / B for the sources, so supervised runs happen too),
  .samples.tsv, .truth.tsv, .anc.tsv.gz, .refs.txt

For the model set, A_frac_bp is the strain's target q and A_frac_snp the realised share of
A alleles over the SNPs.

Usage:
  python control_admixture.py bench5/single_tetrad --sources 20 --outdir bench5/control
"""
import argparse, gzip, os, sys
import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('prefix', help='haploswarm.py output prefix (run with --include-parents)')
    ap.add_argument('--sources', type=int, default=20, help='clonal copies of each parent per set (default 20)')
    ap.add_argument('--gens', help='generations to write (default: all sampled)')
    ap.add_argument('--outdir', help='output folder (default: <prefix>_control)')
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    name = os.path.basename(a.prefix)
    outdir = a.outdir or a.prefix + '_control'

    # genotypes (0/1) of every sample at every SNP, from the haploid VCF
    ids, rows = None, []
    for l in gzip.open(a.prefix + '.vcf.gz', 'rt'):
        if l.startswith('##'):
            continue
        c = l.rstrip('\n').split('\t')
        if l.startswith('#'):
            ids = c[9:]; hdr = c[:9]; continue
        rows.append((c[:9], [int(x) for x in c[9:]]))
    fixed = [r[0] for r in rows]
    G = np.array([r[1] for r in rows], dtype=np.uint8)          # SNPs x samples
    if 'parentA' not in ids or 'parentB' not in ids:
        sys.exit('Error: run haploswarm.py with --include-parents (needs parentA and parentB)')
    pa, pb = G[:, ids.index('parentA')], G[:, ids.index('parentB')]
    keep = np.where(pa != pb)[0]                                  # SNPs where the parents differ
    fixed = [fixed[j] for j in keep]; G = G[keep]; pa, pb = pa[keep], pb[keep]
    anc_rows = [l.rstrip('\n').split('\t') for l in gzip.open(a.prefix + '.anc.tsv.gz', 'rt')]
    anc_ids = anc_rows[0][2:]
    anc_pos = {(r[0], r[1]): i for i, r in enumerate(anc_rows[1:])}
    ANC = np.array([anc_rows[1 + anc_pos[(f[0], f[1])]][2:] for f in fixed]).astype(np.uint8)   # 1 = A
    ANC = ANC[:, [anc_ids.index(s) for s in ids]]

    samples = [l.rstrip('\n').split('\t') for l in open(a.prefix + '.samples.tsv')][1:]
    group = {s[1]: s[2] for s in samples}
    truth = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in list(open(a.prefix + '.truth.tsv'))[1:]}
    gens = sorted({int(group[s].rsplit('_g', 1)[1]) for s in ids if s in group and group[s].startswith('swarm_g')})
    if a.gens:
        gens = [g for g in gens if g in {int(x) for x in a.gens.split(',')}]
    nS = a.sources
    print(f'{len(keep)} SNPs where the parents differ; {nS} clonal copies of each parent per set')

    def write(d, tag, names, groups, kinds, labels, geno, anc, qbp, mats):
        os.makedirs(d, exist_ok=True)
        p = os.path.join(d, tag)
        m = geno.shape[0]
        with open(p + '.samples.tsv', 'w') as f:
            f.write('order\tsample\tgroup\tkind\n')
            for i, (s, g, k) in enumerate(zip(names, groups, kinds), 1):
                f.write(f'{i}\t{s}\t{g}\t{k}\n')
        with open(p + '.truth.tsv', 'w') as f:
            f.write('sample\tgroup\tkind\tMAT\tA_frac_snp\tA_frac_bp\tn_tracts\n')
            for i, s in enumerate(names):
                f.write(f'{s}\t{groups[i]}\t{kinds[i]}\t{mats[i]}\t{anc[:, i].mean():.5f}\t{qbp[i]:.5f}\t-\n')
        with open(p + '.pop', 'w') as f:
            f.write(''.join((l if k == 'ref' else '-') + '\n' for l, k in zip(labels, kinds)))
        with open(p + '.refs.txt', 'w') as f:
            f.write(''.join(s + '\n' for s, k in zip(names, kinds) if k == 'ref'))
        with open(p + '.geno', 'w') as f:
            for j in range(m):
                f.write(''.join(map(str, geno[j])) + '\n')
        with gzip.open(p + '.anc.tsv.gz', 'wt', compresslevel=4) as f:
            f.write('chrom\tpos\t' + '\t'.join(names) + '\n')
            for j in range(m):
                f.write(f'{fixed[j][0]}\t{fixed[j][1]}\t' + '\t'.join(map(str, anc[j])) + '\n')
        for suffix, fmt in (('.dip.vcf.gz', lambda x: f'{x}/{x}'), ('.vcf.gz', lambda x: str(x))):
            with gzip.open(p + suffix, 'wt', compresslevel=4) as f:
                f.write('##fileformat=VCFv4.2\n##source=control_admixture.py\n')
                f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                f.write('\t'.join(hdr) + '\t' + '\t'.join(names) + '\n')
                for j in range(m):
                    f.write('\t'.join(fixed[j]) + '\t' + '\t'.join(fmt(x) for x in geno[j]) + '\n')

    src_names = [f'srcA_{i}' for i in range(1, nS + 1)] + [f'srcB_{i}' for i in range(1, nS + 1)]
    src_geno = np.hstack([np.tile(pa[:, None], nS), np.tile(pb[:, None], nS)])
    src_anc = np.hstack([np.ones((len(pa), nS), np.uint8), np.zeros((len(pa), nS), np.uint8)])
    src_groups = ['ref_A'] * nS + ['ref_B'] * nS
    for g in gens:
        cols = [i for i, s in enumerate(ids) if group.get(s) == f'swarm_g{g}']
        sw = [ids[i] for i in cols]
        q = np.array([float(truth[s][5]) for s in sw])
        mats = [truth[s][3] for s in sw]
        # 1. the real swarm strains + clonal sources
        write(os.path.join(outdir, 'swarm_sources', f'g{g:02d}'), f'{name}_g{g:02d}',
              sw + src_names, [f'swarm_g{g}'] * len(sw) + src_groups, ['sim'] * len(sw) + ['ref'] * 2 * nS,
              ['-'] * len(sw) + ['A'] * nS + ['B'] * nS,
              np.hstack([G[:, cols], src_geno]), np.hstack([ANC[:, cols], src_anc]),
              np.concatenate([q, np.ones(nS), np.zeros(nS)]), mats + ['-'] * 2 * nS)
        # 2. model-conforming strains (independent SNPs, same q) + the same clonal sources
        fromA = (rng.random((len(pa), len(sw))) < q[None, :]).astype(np.uint8)
        geno = np.where(fromA == 1, pa[:, None], pb[:, None]).astype(np.uint8)
        mnames = [f'model_g{g}_{s.split("_", 1)[1]}' for s in sw]
        write(os.path.join(outdir, 'model', f'g{g:02d}'), f'{name}_model_g{g:02d}',
              mnames + src_names, [f'model_g{g}'] * len(sw) + src_groups, ['sim'] * len(sw) + ['ref'] * 2 * nS,
              ['-'] * len(sw) + ['A'] * nS + ['B'] * nS,
              np.hstack([geno, src_geno]), np.hstack([fromA, src_anc]),
              np.concatenate([q, np.ones(nS), np.zeros(nS)]), ['-'] * len(sw) + ['-'] * 2 * nS)
        print(f'generation {g:2d}: {len(sw)} strains (q mean {q.mean():.3f}, sd {q.std():.3f}) + {2 * nS} clonal sources '
              f'-> {outdir}/{{swarm_sources,model}}/g{g:02d}/', flush=True)


if __name__ == '__main__':
    main()
