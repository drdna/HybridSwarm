#!/bin/bash
#SBATCH --job-name=admix_bygen
#SBATCH --array=0-11
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --mem=24GB
#SBATCH -A coa_farman_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user farman@uky.edu
#
# Run PLINK2 + ADMIXTURE on every generation folder made by split_generations.py.
# As a SLURM array, each task does one generation (task 0 = first g* folder, ...);
# set --array to 0-(number of generations - 1). Without SLURM it loops over all.
#
#   sbatch run_admixture_by_gen.sh tetrad_swarmonly                # defaults below
#   sbatch run_admixture_by_gen.sh tetrad_swarmonly 5 10 allmaf    # KMAX RUNS PRUNE
#
# PRUNE decides how SNPs are chosen (MAF >= 0.05; LD pruning is --indep-pairwise 50 10 0.1):
#   allmaf  : MAF on all samples of the generation, no LD pruning (default)
#   all     : MAF + LD pruning on all samples (a standard pipeline)
#   none    : every SNP polymorphic in the generation
#   refsmaf / refs : MAF (+ LD pruning) on the reference isolates only; only for folders
#             made with references. Not recommended (see README).
#
# Folders made with split_generations.py --no-refs (swarm only) get unsupervised runs only;
# supervised runs need A / B references in the .pop file (e.g. the *_parents folders).
#
# Output, in each g<gg>/ folder:
#   <tag>_ADMIX.bed/.bim/.fam                  the SNPs ADMIXTURE used
#   unsup/<tag>_ADMIX_K<K>_run<r>.Q/.P/.log    unsupervised, K = 1..KMAX, RUNS seeds each
#   sup/<tag>_ADMIX_K2_run<r>.Q/.P/.log        supervised, only if the .pop file has A / B references

set -euo pipefail
DIR=${1:?usage: run_admixture_by_gen.sh <folder from split_generations.py> [KMAX] [RUNS] [PRUNE]}
KMAX=${2:-5}
RUNS=${3:-10}
PRUNE=${4:-allmaf}
THREADS=${SLURM_CPUS_PER_TASK:-1}

PLINK_IMG=/share/singularity/images/ccs/conda/amd-conda25-rocky9.sinf
ADMIXTURE_IMG=/share/singularity/images/ccs/conda/amd-conda15-rocky8.sinf
plink2()    { singularity run --app plink200a612 "$PLINK_IMG" plink2 "$@"; }
admixture() { singularity run --app admixture130 "$ADMIXTURE_IMG" admixture "$@"; }

GENS=( $(ls -d "$DIR"/g*/ | sort) )
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  (( SLURM_ARRAY_TASK_ID < ${#GENS[@]} )) || { echo "no generation folder for task $SLURM_ARRAY_TASK_ID"; exit 0; }
  GENS=( "${GENS[$SLURM_ARRAY_TASK_ID]}" )
fi

for g in "${GENS[@]}"; do
  pushd "$g" > /dev/null
  tag=$(basename "$(ls *.dip.vcf.gz)" .dip.vcf.gz)
  echo "== $tag  (PRUNE=$PRUNE, K=1..$KMAX, $RUNS runs)"

  plink2 --vcf "$tag.dip.vcf.gz" --make-bed --out "$tag" --threads "$THREADS" > /dev/null
  if [[ "$PRUNE" == refs* && ! -s "$tag.refs.txt" ]]; then
    echo "PRUNE=$PRUNE needs reference isolates, but this folder has none; use allmaf"; exit 1
  fi
  case "$PRUNE" in
    allmaf)  plink2 --bfile "$tag" --maf 0.05 --write-snplist --out "${tag}_sel" --threads "$THREADS" > /dev/null
             plink2 --bfile "$tag" --extract "${tag}_sel.snplist" --make-bed --out "${tag}_ADMIX" > /dev/null ;;
    all)     plink2 --bfile "$tag" --maf 0.05 --indep-pairwise 50 10 0.1 --bad-ld --out "${tag}_sel" --threads "$THREADS" > /dev/null
             plink2 --bfile "$tag" --extract "${tag}_sel.prune.in" --make-bed --out "${tag}_ADMIX" > /dev/null ;;
    none)    plink2 --bfile "$tag" --make-bed --out "${tag}_ADMIX" > /dev/null ;;
    refsmaf) plink2 --bfile "$tag" --keep "$tag.refs.txt" --maf 0.05 --write-snplist --out "${tag}_sel" --threads "$THREADS" > /dev/null
             plink2 --bfile "$tag" --extract "${tag}_sel.snplist" --make-bed --out "${tag}_ADMIX" > /dev/null ;;
    refs)    plink2 --bfile "$tag" --keep "$tag.refs.txt" --maf 0.05 --indep-pairwise 50 10 0.1 --bad-ld --out "${tag}_sel" --threads "$THREADS" > /dev/null
             plink2 --bfile "$tag" --extract "${tag}_sel.prune.in" --make-bed --out "${tag}_ADMIX" > /dev/null ;;
    *)       echo "PRUNE must be allmaf, all, none, refsmaf or refs"; exit 1 ;;
  esac
  echo "   $(wc -l < "${tag}_ADMIX.bim") SNPs, $(wc -l < "${tag}_ADMIX.fam") samples"

  # unsupervised
  mkdir -p unsup && cd unsup
  ln -sf "../${tag}_ADMIX.bed" "../${tag}_ADMIX.bim" "../${tag}_ADMIX.fam" .
  for K in $(seq 1 "$KMAX"); do
    for r in $(seq 1 "$RUNS"); do
      # a failed run is logged and skipped; without this, set -e would end the whole task
      if admixture -j"$THREADS" -s "$r" --cv "${tag}_ADMIX.bed" "$K" > "${tag}_ADMIX_K${K}_run${r}.log" 2>&1; then
        mv "${tag}_ADMIX.$K.Q" "${tag}_ADMIX_K${K}_run${r}.Q"; mv "${tag}_ADMIX.$K.P" "${tag}_ADMIX_K${K}_run${r}.P"
      else
        echo "   K=$K run $r FAILED (exit $?): see unsup/${tag}_ADMIX_K${K}_run${r}.log"
        rm -f "${tag}_ADMIX.$K.Q" "${tag}_ADMIX.$K.P"
      fi
    done
  done
  cd ..

  # supervised (K = 2; the .pop file must sit next to the .bed with the same name)
  if ! grep -qE '^(A|B)$' "$tag.pop"; then
    echo "   no reference isolates: supervised runs skipped"
    popd > /dev/null; continue
  fi
  mkdir -p sup && cd sup
  ln -sf "../${tag}_ADMIX.bed" "../${tag}_ADMIX.bim" "../${tag}_ADMIX.fam" .
  cp "../$tag.pop" "${tag}_ADMIX.pop"
  for r in $(seq 1 "$RUNS"); do
    if admixture -j"$THREADS" -s "$r" --supervised "${tag}_ADMIX.bed" 2 > "${tag}_ADMIX_K2_run${r}.log" 2>&1; then
      mv "${tag}_ADMIX.2.Q" "${tag}_ADMIX_K2_run${r}.Q"; mv "${tag}_ADMIX.2.P" "${tag}_ADMIX_K2_run${r}.P"
    else
      echo "   supervised run $r FAILED (exit $?): see sup/${tag}_ADMIX_K2_run${r}.log"
      rm -f "${tag}_ADMIX.2.Q" "${tag}_ADMIX.2.P"
    fi
  done
  cd ..
  popd > /dev/null
done
