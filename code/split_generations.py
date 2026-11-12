#!/usr/bin/env python3
"""Split a haploswarm.py data set into one input set per generation.

Each generation's swarm strains are kept together with all reference isolates,
so every generation is analysed on its own with the same references. SNPs that
are monomorphic within a subset are dropped from that subset.

For generation g the files are written to <outdir>/g<gg>/<name>_g<gg>.*:
  dip.vcf.gz, vcf.gz, geno, str (+ mainparams, extraparams), pop,
  samples.tsv, truth.tsv, anc.tsv.gz, refs.txt (reference IDs, for plink --keep)

Usage:
  python split_generations.py bench2/tetrad                    # all generations
  python split_generations.py bench2/tetrad --gens 1,2,3,5,8,11 --outdir bench2/tetrad_bygen
  python split_generations.py bench2/tetrad --max-strains 40    # also subsample each generation
"""
import argparse, gzip, os, random, re, sys
import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('prefix', help='haploswarm output prefix, e.g. bench2/tetrad')
    ap.add_argument('--gens', help='generations to write (default: all sampled generations)')
    ap.add_argument('--outdir', help='output folder (default: <prefix>_bygen)')
    ap.add_argument('--max-strains', type=int, default=0, help='randomly keep at most this many strains per generation')
    ap.add_argument('--no-refs', action='store_true', help='leave the reference isolates out')
    ap.add_argument('--pop', help='with a multi-population model: only this population (default: all populations of a generation together)')
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    name = os.path.basename(a.prefix)
    outdir = a.outdir or a.prefix + '_bygen'

    samples = [l.rstrip('\n').split('\t') for l in open(a.prefix + '.samples.tsv')][1:]   # order sample group kind
    ids = [s[1] for s in samples]
    gen_of = {s[1]: int(s[2].rsplit('_g', 1)[1]) for s in samples
              if s[3] == 'sim' and (not a.pop or s[2].rsplit('_g', 1)[0] == a.pop)}
    refs = [i for i, s in enumerate(samples) if s[3] == 'ref']
    gens = sorted(set(gen_of.values())) if not a.gens else [int(x) for x in a.gens.split(',')]

    # read everything once
    print('reading', a.prefix, '...', flush=True)
    vcf_hdr, vcf_rows = [], []
    for l in gzip.open(a.prefix + '.dip.vcf.gz', 'rt'):
        if l.startswith('##'):
            vcf_hdr.append(l)
        elif l.startswith('#'):
            assert l.rstrip('\n').split('\t')[9:] == ids, 'VCF sample order differs from samples.tsv'
        else:
            vcf_rows.append(l.rstrip('\n').split('\t'))
    fixed = [r[:9] for r in vcf_rows]
    G = np.array([[0 if x == '0/0' else 1 if x == '1/1' else 9 for x in r[9:]] for r in vcf_rows], dtype=np.uint8)  # snps x samples
    anc_lines = [l.rstrip('\n').split('\t') for l in gzip.open(a.prefix + '.anc.tsv.gz', 'rt')]
    assert anc_lines[0][2:] == ids, 'ancestry table sample order differs from samples.tsv'
    anc = np.array([r[2:] for r in anc_lines[1:]])
    truth = {l.split('\t')[0]: l for l in list(open(a.prefix + '.truth.tsv'))[1:]}
    truth_hdr = open(a.prefix + '.truth.tsv').readline()
    pop = [l.strip() for l in open(a.prefix + '.pop')]
    str_lines = open(a.prefix + '.str').read().split('\n')
    has_map = 'LINKAGE 1' in open(a.prefix + '.extraparams').read()
    extraparams = open(a.prefix + '.extraparams').read()
    mainparams = open(a.prefix + '.mainparams').read()
    cm = None
    if has_map:          # rebuild cumulative cM per chromosome so distances can be recomputed after dropping SNPs
        d = str_lines[1].split()
        cm, c = [], 0.0
        for x in d:
            c = 0.0 if x == '-1' else c + float(x)
            cm.append(c)
    chrom = [f[0] for f in fixed]

    for g in gens:
        sel = [i for i, s in enumerate(ids) if gen_of.get(s) == g]
        if not sel:
            print(f'generation {g}: no sampled strains, skipped'); continue
        if a.max_strains and len(sel) > a.max_strains:
            sel = sorted(rng.sample(sel, a.max_strains))
        cols = sel + ([] if a.no_refs else refs)
        sub = G[:, cols]
        keep = np.where((sub.min(1) != sub.max(1)))[0]
        tag = f'{name}_g{g:02d}'
        d = os.path.join(outdir, f'g{g:02d}')
        os.makedirs(d, exist_ok=True)
        p = os.path.join(d, tag)
        sid = [ids[i] for i in cols]

        with open(p + '.samples.tsv', 'w') as f:
            f.write('order\tsample\tgroup\tkind\n')
            for n, i in enumerate(cols, 1):
                f.write(f'{n}\t' + '\t'.join(samples[i][1:]) + '\n')
        with open(p + '.truth.tsv', 'w') as f:
            f.write(truth_hdr + ''.join(truth[s] for s in sid))
        with open(p + '.pop', 'w') as f:
            f.write(''.join(pop[i] + '\n' for i in cols))
        with open(p + '.refs.txt', 'w') as f:
            f.write(''.join(ids[i] + '\n' for i in cols if i in refs))
        with gzip.open(p + '.anc.tsv.gz', 'wt', compresslevel=4) as f:
            f.write('chrom\tpos\t' + '\t'.join(sid) + '\n')
            for j in keep:
                f.write(f'{anc_lines[j + 1][0]}\t{anc_lines[j + 1][1]}\t' + '\t'.join(anc[j, cols]) + '\n')
        with open(p + '.geno', 'w') as f:
            for j in keep:
                f.write(''.join(map(str, sub[j])) + '\n')
        for suffix, fmt in (('.dip.vcf.gz', lambda x: f'{x}/{x}'), ('.vcf.gz', lambda x: str(x))):
            with gzip.open(p + suffix, 'wt', compresslevel=4) as f:
                f.write(''.join(vcf_hdr) + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(sid) + '\n')
                for j in keep:
                    f.write('\t'.join(fixed[j]) + '\t' + '\t'.join(fmt(x) for x in sub[j]) + '\n')
        # STRUCTURE: rows are samples; data rows start after the marker-name row (and map row)
        start = 2 if has_map else 1
        rows = {r.split(' ', 1)[0]: r.split(' ') for r in str_lines[start:] if r.strip()}
        with open(p + '.str', 'w') as f:
            names = str_lines[0].split()
            f.write(' '.join(names[j] for j in keep) + '\n')
            if has_map:
                f.write(' '.join('-1' if k == 0 or chrom[keep[k]] != chrom[keep[k - 1]] else f'{cm[keep[k]] - cm[keep[k - 1]]:.6f}'
                                 for k in range(len(keep))) + '\n')
            for s in sid:
                r = rows[s]
                f.write(' '.join(r[:3] + [r[3 + j] for j in keep]) + '\n')
        mp = re.sub(r'#define NUMINDS \d+', f'#define NUMINDS {len(cols)}', mainparams)
        mp = re.sub(r'#define NUMLOCI \d+', f'#define NUMLOCI {len(keep)}', mp)
        mp = re.sub(r'#define INFILE \S+', f'#define INFILE {tag}.str', mp)
        mp = re.sub(r'#define OUTFILE \S+', f'#define OUTFILE {tag}.structure', mp)
        open(p + '.mainparams', 'w').write(mp)
        open(p + '.extraparams', 'w').write(extraparams)
        print(f'generation {g:2d}: {len(sel)} strains + {0 if a.no_refs else len(refs)} references, '
              f'{len(keep)} polymorphic SNPs -> {p}.*', flush=True)


if __name__ == '__main__':
    main()
