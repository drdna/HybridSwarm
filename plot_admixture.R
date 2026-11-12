#!/usr/bin/env Rscript
# Plot ADMIXTURE results against the simulated truth from make_admix_inputs.py
#
# Works with any file naming, e.g. bench_ADMIX_K2_run1.Q or bench_thin.2.Q:
# every <prefix>*.Q file in --dir is read, K is taken from its number of
# columns, and its log is the file with the same name ending .log or .out
# (ADMIXTURE writes the CV error and log-likelihood only to the screen, so
# save it: admixture --cv x.bed 2 | tee x_K2_run1.log). With several runs
# per K, the run with the highest log-likelihood is plotted.
#
#   Rscript plot_admixture.R --prefix bench_ADMIX --truth bench_thin.truth.tsv \
#       --samples bench_thin.samples.tsv --anc bench_thin.anc.tsv --label "unsupervised, LD-pruned"
#
# Options:
#   --prefix   start of the .Q file names (required)
#   --dir      folder with the .Q and log files (default .)
#   --fam      PLINK .fam used for the runs; gives the sample order (default <dir>/<prefix>.fam)
#   --bim      PLINK .bim used for the runs; the SNPs ADMIXTURE saw (default: .fam name with .bim)
#   --samples  sample table from make_admix_inputs.py
#   --truth    truth table from make_admix_inputs.py
#   --anc      per-SNP ancestry from make_admix_inputs.py (.anc.tsv). With --bim, the truth is
#              recomputed over only the SNPs in the .bim, e.g. after LD pruning. Recommended.
#   --out      output prefix (default <dir>/<prefix>_plots)
#   --label    text added to plot titles
#   --basis    bp (default): truth = share of the genome inherited from A; snp: share of the SNPs used
#   --kmin, --kmax  only use runs with K in this range (default: every K found)
#
# Writes <out>.pdf (one plot per page), <out>_comparison.tsv (per strain, K = 2)
# and <out>_runs.tsv (log-likelihood and CV error of every run).

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

# ---- arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opt <- list(dir = ".", label = "")
i <- 1
while (i <= length(args)) {
  if (i == length(args) || startsWith(args[i + 1], "--")) stop("Missing value for ", args[i])
  opt[[sub("^--", "", args[i])]] <- args[i + 1]
  i <- i + 2
}
if (is.null(opt$prefix)) stop("Usage: Rscript plot_admixture.R --prefix <start of .Q names> --samples x.samples.tsv --truth x.truth.tsv [--anc x.anc.tsv]")
for (k in c("samples", "truth")) if (is.null(opt[[k]])) stop("--", k, " is required")
if (is.null(opt$fam)) opt$fam <- file.path(opt$dir, paste0(opt$prefix, ".fam"))
if (is.null(opt$bim)) opt$bim <- sub("\\.fam$", ".bim", opt$fam)
if (is.null(opt$out)) opt$out <- file.path(opt$dir, paste0(opt$prefix, "_plots"))
ttl <- function(x) if (nzchar(opt$label)) paste0(x, " (", opt$label, ")") else x

# ---- samples in .fam order, joined to the truth ---------------------------------
if (!file.exists(opt$fam)) stop("Can't find ", opt$fam, "; give the .fam used for the runs with --fam")
fam <- read.table(opt$fam, stringsAsFactors = FALSE)[, 2]
samples <- read.delim(opt$samples, stringsAsFactors = FALSE)
truth <- read.delim(opt$truth, stringsAsFactors = FALSE)
info <- data.frame(sample = fam, stringsAsFactors = FALSE) %>%
  left_join(samples %>% select(sample, group, kind), by = "sample") %>%
  left_join(truth %>% select(sample, MAT, A_frac_snp, A_frac_bp, n_tracts), by = "sample") %>%
  mutate(kind = ifelse(is.na(kind), "other", kind),
         group = ifelse(is.na(group), "other", group),
         A_frac_snp = suppressWarnings(as.numeric(A_frac_snp)))
extra <- info$sample[info$kind == "other"]
if (length(extra)) message("Not in the sample table (plotted as 'other'): ", paste(extra, collapse = ", "))

truth_note <- "true A fraction at the SNPs in the truth table"
nsnp_used <- NA
if (!is.null(opt$anc)) {
  if (!file.exists(opt$bim)) stop("--anc needs the .bim of the runs; give it with --bim")
  bim <- read.table(opt$bim, stringsAsFactors = FALSE)
  # The .P files have one row per SNP ADMIXTURE used: make sure the .bim is the same SNP set
  pfile <- list.files(opt$dir, pattern = paste0("^", gsub("([.()+^$|\\[\\]])", "\\\\\\1", opt$prefix), ".*\\.P$"),
                      full.names = TRUE)[1]
  if (!is.na(pfile)) {
    np <- length(readLines(pfile, warn = FALSE))
    if (np != nrow(bim)) {
      cand <- list.files(unique(c(dirname(opt$bim), dirname(opt$fam), opt$dir)), pattern = "\\.bim$", full.names = TRUE)
      cand <- cand[sapply(cand, function(f) length(readLines(f, warn = FALSE)) == np)]
      if (!length(cand)) stop(sprintf("%s has %d SNPs but ADMIXTURE used %d (rows of %s). Give the .fam/.bim of the ADMIXTURE input with --fam / --bim.",
                                      opt$bim, nrow(bim), np, basename(pfile)))
      message(sprintf("Note: %s has %d SNPs but ADMIXTURE used %d; using %s instead (pass --fam %s next time).",
                      basename(opt$bim), nrow(bim), np, basename(cand[1]), sub("\\.bim$", ".fam", basename(cand[1]))))
      opt$bim <- cand[1]
      bim <- read.table(opt$bim, stringsAsFactors = FALSE)
    }
  }
  anc <- read.delim(opt$anc, check.names = FALSE, colClasses = "character")
  # match SNPs on chromosome + position when the ancestry table has a chrom column, else on position
  key_anc <- if ("chrom" %in% names(anc)) paste(anc$chrom, anc$pos) else anc$pos
  key_bim <- if ("chrom" %in% names(anc)) paste(bim$V1, bim$V4) else as.character(bim$V4)
  used <- anc[key_anc %in% key_bim, , drop = FALSE]
  nsnp_used <- nrow(bim)
  if (nrow(used) < nrow(bim)) warning(nrow(bim) - nrow(used), " SNPs in the .bim are not in ", opt$anc)
  info$A_frac_snp <- sapply(info$sample, function(s)
    if (s %in% names(used)) { v <- used[[s]]; v <- v[v != "."]; if (length(v)) mean(v == "1") else NA } else NA)
  truth_note <- sprintf("true A fraction at the %d SNPs ADMIXTURE used", nrow(used))
}
snp_note <- truth_note
info$A_frac_bp <- suppressWarnings(as.numeric(info$A_frac_bp))
if (is.null(opt$basis)) opt$basis <- "bp"
if (!opt$basis %in% c("bp", "snp")) stop("--basis must be bp or snp")
if (opt$basis == "bp") {
  info$truth <- info$A_frac_bp
  truth_note <- "true A fraction of the whole genome (bp)"
} else info$truth <- info$A_frac_snp
message("Truth: ", truth_note, if (opt$basis == "bp") paste0(" (also reported: ", snp_note, ")") else "")

# ---- runs ------------------------------------------------------------------------
esc <- gsub("([.()+^$|\\[\\]])", "\\\\\\1", opt$prefix)
qfiles <- list.files(opt$dir, pattern = paste0("^", esc, ".*\\.Q$"), full.names = TRUE)
if (!length(qfiles)) stop("No ", opt$prefix, "*.Q files in ", opt$dir)

refA <- which(info$kind == "ref" & info$group == "ref_A")
orient <- function(q) {
  # cluster numbers are arbitrary: put the cluster that best matches ancestry A first
  score <- if (length(refA)) colMeans(q[refA, , drop = FALSE]) else
    apply(q, 2, function(x) suppressWarnings(cor(x, info$truth, use = "complete.obs")))
  score[is.na(score)] <- -Inf
  first <- which.max(score)
  rest <- setdiff(seq_len(ncol(q)), first)
  # then the cluster that best matches the B references (if any), then the rest by overall size
  refB <- which(info$kind == "ref" & info$group == "ref_B")
  if (length(rest) > 1) {
    sb <- if (length(refB)) colMeans(q[refB, rest, drop = FALSE]) else rep(0, length(rest))
    rest <- rest[order(-sb, -colMeans(q[, rest, drop = FALSE]))]
  }
  q <- q[, c(first, rest), drop = FALSE]
  colnames(q) <- paste0("C", seq_len(ncol(q)))
  q
}
runs <- bind_rows(lapply(qfiles, function(f) {
  stem <- sub("\\.Q$", "", f)
  log <- c(paste0(stem, ".log"), paste0(stem, ".out"))
  log <- log[file.exists(log)][1]
  l <- if (!is.na(log)) readLines(log, warn = FALSE) else character()
  num <- function(re) { x <- grep(re, l, value = TRUE); if (length(x)) as.numeric(sub(".*: *", "", x[length(x)])) else NA }
  data.frame(file = basename(f), K = ncol(read.table(f, nrows = 1)),
             loglik = num("^Loglikelihood:"), cv = num("CV error \\(K="), stringsAsFactors = FALSE)
}))
kmin <- if (is.null(opt$kmin)) -Inf else as.numeric(opt$kmin)
kmax <- if (is.null(opt$kmax)) Inf else as.numeric(opt$kmax)
runs <- runs %>% filter(K >= kmin, K <= kmax)
if (!nrow(runs)) stop("No runs with K between ", kmin, " and ", kmax, " in ", opt$dir)
runs <- runs %>% group_by(K) %>% mutate(best = row_number() == which.max(replace(loglik, is.na(loglik), -Inf))) %>%
  ungroup() %>% arrange(K, file)
write.table(runs, paste0(opt$out, "_runs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
Ks <- sort(unique(runs$K))
read_q <- function(file) {
  q <- as.matrix(read.table(file.path(opt$dir, file)))
  if (nrow(q) != nrow(info)) stop(file, " has ", nrow(q), " rows but ", opt$fam, " has ", nrow(info), " samples")
  orient(q)
}
Q <- setNames(lapply(Ks, function(K) read_q(runs$file[runs$K == K & runs$best])), Ks)
# proportions of the best run for every K (used by compare_generations.R)
write.table(bind_rows(lapply(Ks, function(K)
  data.frame(sample = info$sample, group = as.character(info$group), kind = info$kind,
             truth = info$truth, A_frac_snp = info$A_frac_snp, K = K, Q[[as.character(K)]], check.names = FALSE) %>%
    pivot_longer(starts_with("C"), names_to = "cluster", values_to = "q"))),
  paste0(opt$out, "_q.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(nrow(runs), " runs for K = ", paste(Ks, collapse = ","), "; plotting the highest-likelihood run for each K")

plots <- list()
sub_txt <- if (is.na(nsnp_used)) "" else sprintf("ADMIXTURE used %d SNPs. ", nsnp_used)

# 1. Cross-validation error (and log-likelihood) for every run
if (any(!is.na(runs$cv))) {
  cvm <- runs %>% filter(!is.na(cv)) %>% group_by(K) %>% summarise(cv = mean(cv), .groups = "drop")
  plots$cv <- ggplot(runs %>% filter(!is.na(cv)), aes(K, cv)) +
    geom_line(data = cvm) +
    geom_jitter(width = 0.08, height = 0, alpha = 0.6) +
    geom_point(data = cvm[which.min(cvm$cv), ], colour = "firebrick", size = 4, shape = 1, stroke = 1.2) +
    scale_x_continuous(breaks = Ks) +
    labs(x = "K", y = "Cross-validation error", title = ttl("ADMIXTURE cross-validation"),
         subtitle = paste0(sub_txt, "Points = runs, line = mean. Lowest mean CV error at K = ",
                           cvm$K[which.min(cvm$cv)], " (circled).")) +
    theme_bw()
} else message("No 'CV error' lines in the logs; skipping the CV plot")
if (any(!is.na(runs$loglik))) {
  plots$loglik <- ggplot(runs %>% filter(!is.na(loglik)), aes(K, loglik, colour = best)) +
    geom_jitter(width = 0.08, height = 0) +
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "firebrick"), labels = c("other runs", "best run (plotted)")) +
    scale_x_continuous(breaks = Ks) +
    labs(x = "K", y = "Log-likelihood", colour = NULL, title = ttl("Log-likelihood of each run"),
         subtitle = "Runs that differ at the same K got stuck in different solutions") +
    theme_bw() + theme(legend.position = "bottom")
}

# 2. Admixture proportions for every K, with the truth on top
gen_of <- function(g) suppressWarnings(as.numeric(sub(".*_g", "", g)))   # swarm_g10 -> 10
simg <- unique(info$group[info$kind == "sim"])
grp_levels <- unique(c("ref_A", "ref_B", simg[order(gen_of(simg), simg)], "other"))
info$group <- factor(info$group, levels = grp_levels[grp_levels %in% info$group])
ord <- order(info$group, -info$truth)
info$x <- match(seq_len(nrow(info)), ord)
pal <- function(K) c("#2166AC", "#D6604D", "#4DAC26", "#8073AC", "#F1A340", "#35978F", "#C51B7D", "#878787",
                     "#B8E186", "#FDB863", "#5AAE61", "#E7D4E8")[seq_len(K)]
truth_long <- info %>% filter(!is.na(truth)) %>%
  transmute(sample, group, x, track = "Truth", C1 = truth, C2 = 1 - truth) %>%
  pivot_longer(c(C1, C2), names_to = "cluster", values_to = "q")
q_long <- bind_rows(lapply(Ks, function(K)
  data.frame(info[, c("sample", "group", "x")], Q[[as.character(K)]], check.names = FALSE) %>%
    pivot_longer(starts_with("C"), names_to = "cluster", values_to = "q") %>%
    mutate(track = paste0("K = ", K))))
bars <- bind_rows(truth_long, q_long) %>%
  mutate(track = factor(track, levels = c("Truth", paste0("K = ", Ks))),
         cluster = factor(cluster, levels = paste0("C", seq_len(max(Ks, 2)))))
plots$bars <- ggplot(bars, aes(x, q, fill = cluster)) +
  geom_col(width = 1) +
  facet_grid(track ~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal(max(Ks, 2)), labels = function(l) ifelse(l == "C1", "C1 (= A)", l)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Strains (within each group, sorted by true A fraction)", y = "Ancestry proportion", fill = "Cluster",
       title = ttl("Admixture proportions (best run per K)"),
       subtitle = paste0(sub_txt, "Top row: ", truth_note, " (blue = A). C1 is the cluster matching ancestry A.")) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = unit(1, "pt"),
        strip.text.y = element_text(angle = 0), strip.text.x = element_text(size = 7), legend.position = "bottom")

# 3. Cluster membership: each strain assigned to its largest cluster
memb <- bind_rows(lapply(Ks, function(K) {
  q <- Q[[as.character(K)]]
  data.frame(group = info$group, K = factor(paste0("K = ", K), levels = paste0("K = ", Ks)),
             cluster = factor(colnames(q)[max.col(q, ties.method = "first")], levels = paste0("C", seq_len(max(Ks)))),
             call = ifelse(apply(q, 1, max) >= 0.9, "Q >= 0.9 (\"pure\")", "Q < 0.9 (admixed)"))
}))
plots$membership <- ggplot(memb, aes(group, fill = cluster, alpha = call)) +
  geom_bar(colour = "grey30", linewidth = 0.2) +
  facet_wrap(~K) +
  scale_fill_manual(values = pal(max(Ks, 2)), drop = FALSE) +
  scale_alpha_manual(values = c(`Q < 0.9 (admixed)` = 0.4, `Q >= 0.9 ("pure")` = 1)) +
  labs(x = NULL, y = "Number of strains", fill = "Largest cluster", alpha = NULL,
       title = ttl("Cluster membership (best run per K)"),
       subtitle = "Each strain assigned to its largest cluster; solid = Q >= 0.9 in that cluster") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

# 4. Estimated vs true A fraction at K = 2, with the spread across runs
if (2 %in% Ks) {
  allq2 <- sapply(runs$file[runs$K == 2], function(f) read_q(f)[, "C1"])
  cmp <- info %>%
    mutate(est_A = Q[["2"]][, "C1"], est_A_runs_sd = apply(as.matrix(allq2), 1, sd),
           error = est_A - truth) %>%
    mutate(error_vs_snp_truth = est_A - A_frac_snp) %>%
    select(sample, group, kind, MAT, A_frac_bp, A_frac_snp, truth, n_tracts, est_A, est_A_runs_sd, error, error_vs_snp_truth)
  write.table(cmp, paste0(opt$out, "_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  sim <- cmp %>% filter(kind == "sim", !is.na(truth))
  # Do the two clusters actually separate the A and B references? If not, "C1" is not ancestry A.
  refq <- cmp %>% filter(kind == "ref") %>% group_by(group) %>% summarise(C1 = mean(est_A), .groups = "drop")
  if (all(c("ref_A", "ref_B") %in% refq$group) && refq$C1[refq$group == "ref_B"] > 0.5) {
    msg <- sprintf("WARNING: at K = 2 the B references are also mostly in C1 (mean Q %.2f), so the clusters do not separate A from B; C1 is not ancestry A.",
                   refq$C1[refq$group == "ref_B"])
    message(msg)
    opt$label <- paste0(opt$label, if (nzchar(opt$label)) "; " else "", "K = 2 clusters do NOT separate A and B references")
  }
  stats <- sim %>% group_by(group) %>%
    summarise(n = n(), r = suppressWarnings(cor(est_A, truth)), bias = mean(error),
              RMSE = sqrt(mean(error^2)), .groups = "drop") %>%
    mutate(lab = sprintf("n = %d\nr = %.3f\nbias = %+.3f\nRMSE = %.3f", n, r, bias, RMSE))
  cat("\nK = 2, best run, simulated strains:\n"); print(as.data.frame(select(stats, -lab)), digits = 3)
  plots$truth <- ggplot(sim, aes(truth, est_A)) +
    geom_abline(linetype = 2, colour = "grey50") +
    # few SNPs -> truth takes few values, so size points by count; otherwise plain points
    (if (length(unique(sim$truth)) < 40) geom_count(alpha = 0.6, colour = "#2166AC")
     else geom_point(alpha = 0.5, size = 1.2, colour = "#2166AC")) +
    geom_text(data = stats, aes(x = 0.02, y = 0.98, label = lab), hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
    facet_wrap(~group) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = paste0("Truth: ", truth_note), y = "ADMIXTURE estimate of A (K = 2)", size = "Strains",
         title = ttl("ADMIXTURE estimate vs truth"), subtitle = "Dashed line: perfect estimate") +
    theme_bw()
  # accuracy across generations (drift): one point per sampled generation
  gstats <- stats %>% mutate(gen = gen_of(as.character(group))) %>% filter(!is.na(gen))
  if (nrow(gstats) > 1) {
    means <- sim %>% group_by(group) %>%
      summarise(true = mean(truth), estimated = mean(est_A), sd_true = sd(truth), .groups = "drop") %>%
      mutate(gen = gen_of(as.character(group)))
    acc <- gstats %>% select(gen, n, bias, RMSE, r) %>%
      left_join(means %>% select(gen, sd_true), by = "gen") %>%
      pivot_longer(c(bias, RMSE, r, sd_true), names_to = "measure", values_to = "value") %>%
      mutate(measure = factor(measure, levels = c("bias", "RMSE", "r", "sd_true"),
                              labels = c("Bias (mean estimate - truth)", "RMSE per strain",
                                         "Correlation with truth (r)", "SD of true A fraction among strains")))
    plots$by_gen <- ggplot(acc, aes(gen, value)) +
      geom_hline(data = data.frame(measure = factor("Bias (mean estimate - truth)", levels = levels(acc$measure)), y = 0),
                 aes(yintercept = y), linetype = 2, colour = "grey50") +
      geom_line() + geom_point() +
      facet_wrap(~measure, scales = "free_y") +
      scale_x_continuous(breaks = sort(unique(acc$gen))) +
      labs(x = "Generation", y = NULL, title = ttl("Accuracy of ADMIXTURE (K = 2) by generation"),
           subtitle = "Bottom right: how much true ancestry varies among strains, which shrinks as the swarm mixes") +
      theme_bw()
    plots$means <- ggplot(means %>% pivot_longer(c(true, estimated), names_to = "what", values_to = "A"),
                          aes(gen, A, colour = what)) +
      geom_line() + geom_point() +
      scale_colour_manual(values = c(true = "black", estimated = "#D6604D")) +
      scale_x_continuous(breaks = sort(unique(means$gen))) +
      labs(x = "Generation", y = "Mean A fraction of sampled strains", colour = NULL,
           title = ttl("Mean ancestry by generation: truth vs ADMIXTURE (K = 2)"),
           subtitle = "The true mean drifts away from 0.5 because the swarm descends from very few founders") +
      theme_bw() + theme(legend.position = "bottom")
  }
  # error decomposition: markers (SNP truth vs genome truth) and the tool (estimate vs SNP truth)
  dec <- sim %>% filter(!is.na(A_frac_bp), !is.na(A_frac_snp)) %>% group_by(group) %>%
    summarise(`Markers: SNP truth vs genome truth` = sqrt(mean((A_frac_snp - A_frac_bp)^2)),
              `Tool: estimate vs SNP truth` = sqrt(mean((est_A - A_frac_snp)^2)),
              `Total: estimate vs genome truth` = sqrt(mean((est_A - A_frac_bp)^2)),
              `True spread among strains (SD, genome)` = sd(A_frac_bp), .groups = "drop") %>%
    mutate(gen = gen_of(as.character(group))) %>% filter(!is.na(gen))
  if (nrow(dec) > 1) {
    plots$decomp <- ggplot(dec %>% pivot_longer(-c(group, gen), names_to = "what", values_to = "value"),
                           aes(gen, value, colour = what)) +
      geom_line() + geom_point() +
      scale_colour_manual(values = c("#4DAC26", "#D6604D", "black", "grey60")) +
      scale_x_continuous(breaks = sort(unique(dec$gen))) +
      labs(x = "Generation", y = "RMSE (or SD) of A fraction", colour = NULL,
           title = ttl("Where the error comes from (K = 2)"),
           subtitle = paste0("Markers: error even a perfect method would make with these SNPs (", snp_note,
                             "). Tool: ADMIXTURE's own error.\nWhen an error line is above the grey line, the estimates rank strains worse than guessing the mean.")) +
      theme_bw() + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
  }
  plots$error <- ggplot(sim, aes(truth, error)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_jitter(width = 0.005, height = 0, alpha = 0.6, colour = "#D6604D") +
    facet_wrap(~group) +
    labs(x = paste0("Truth: ", truth_note), y = "Estimate - truth",
         title = ttl("Error in ADMIXTURE estimates (K = 2)"),
         subtitle = "Systematic error shows as a trend away from zero") +
    theme_bw()
} else message("No K = 2 runs, so no comparison with the truth")

pdf(paste0(opt$out, ".pdf"), width = 11, height = 7)
for (p in plots) print(p)
invisible(dev.off())
message("Wrote ", opt$out, ".pdf, ", opt$out, "_runs.tsv", if (2 %in% Ks) paste0(" and ", opt$out, "_comparison.tsv"))
