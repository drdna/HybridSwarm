#!/usr/bin/env Rscript
# Compare ADMIXTURE accuracy across generations (and scenarios), after
# run_admixture_by_gen.sh has run on folders made by split_generations.py.
#
#   Rscript compare_generations.R --set tetrad=bench2/tetrad_bygen --set perithecium=bench2/perithecium_bygen
#
# For every g<gg>/unsup and g<gg>/sup folder that has .Q files, it first runs
# plot_admixture.R (one PDF per generation and mode, skipped if already made;
# add --redo yes to remake them), then summarises K = 2 accuracy against the
# whole-genome truth by generation.
#
# Options:
#   --set      label=folder, repeatable (required)
#   --out      output prefix (default: generations_summary)
#   --basis    bp (default) or snp, passed to plot_admixture.R
#   --redo     yes to remake the per-generation plots
#   --kmin, --kmax  only use runs with K in this range, e.g. --kmin 2 --kmax 4 (default: every K)
#
# Writes <out>.pdf and <out>.tsv.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

args <- commandArgs(trailingOnly = TRUE)
sets <- c(); opt <- list(out = "generations_summary", basis = "bp", redo = "no", kmin = "-Inf", kmax = "Inf")
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); val <- args[i + 1]
  if (is.na(val)) stop("Missing value for ", args[i])
  if (key == "set") sets <- c(sets, val) else opt[[key]] <- val
  i <- i + 2
}
kmin <- as.numeric(opt$kmin); kmax <- as.numeric(opt$kmax)
if (!length(sets)) stop("Usage: Rscript compare_generations.R --set label=folder [--set label=folder ...]")

self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
plotter <- file.path(dirname(normalizePath(self)), "plot_admixture.R")
if (!file.exists(plotter)) stop("plot_admixture.R must be in the same folder as this script")

rows <- list(); allq <- list(); allruns <- list(); allcmp <- list()
for (s in sets) {
  label <- sub("=.*", "", s); dir <- sub("^[^=]*=", "", s)
  for (g in sort(list.dirs(dir, recursive = FALSE))) {
    if (!grepl("/g[0-9]+$", g)) next
    gen <- as.integer(sub(".*/g", "", g))
    tag <- sub("\\.samples\\.tsv$", "", list.files(g, pattern = "\\.samples\\.tsv$")[1])
    for (mode in c("unsup", "sup")) {
      d <- file.path(g, mode)
      if (!length(list.files(d, pattern = "\\.Q$"))) next
      out <- file.path(d, paste0(tag, "_", mode))
      cmpf <- paste0(out, "_comparison.tsv")
      if (!file.exists(cmpf) || !file.exists(paste0(out, "_q.tsv")) || opt$redo == "yes") {
        message(sprintf("%s generation %d (%s): plotting", label, gen, mode))
        fam <- file.path(g, paste0(tag, "_ADMIX.fam"))
        res <- system2("Rscript", c(plotter, "--prefix", paste0(tag, "_ADMIX_K"), "--dir", d, "--fam", fam,
                                    "--samples", file.path(g, paste0(tag, ".samples.tsv")),
                                    "--truth", file.path(g, paste0(tag, ".truth.tsv")),
                                    "--anc", file.path(g, paste0(tag, ".anc.tsv.gz")),
                                    "--basis", opt$basis, "--out", out, "--kmin", opt$kmin, "--kmax", opt$kmax,
                                    "--label", shQuote(sprintf("%s, generation %d, %s", label, gen,
                                                               ifelse(mode == "sup", "supervised", "unsupervised")))),
                       stdout = TRUE, stderr = TRUE)
        if (!file.exists(cmpf)) { message(paste(res, collapse = "\n")); next }
      }
      cmp <- read.delim(cmpf)
      runs <- read.delim(paste0(out, "_runs.tsv")) %>% filter(K >= kmin, K <= kmax)
      mlab <- ifelse(mode == "sup", "supervised", "unsupervised")
      key <- paste(label, mlab, sep = " | ")
      allq[[key]] <- bind_rows(allq[[key]], read.delim(paste0(out, "_q.tsv")) %>% filter(K >= kmin, K <= kmax) %>% mutate(gen = gen))
      allruns[[key]] <- bind_rows(allruns[[key]], runs %>% mutate(gen = gen))
      allcmp[[key]] <- bind_rows(allcmp[[key]], cmp %>% mutate(gen = gen))
      sim <- cmp %>% filter(kind == "sim", !is.na(truth))
      refB <- mean(cmp$est_A[cmp$group == "ref_B"]); refA <- mean(cmp$est_A[cmp$group == "ref_A"])
      cvm <- runs %>% filter(!is.na(cv)) %>% group_by(K) %>% summarise(cv = mean(cv), .groups = "drop")
      nsnp <- tryCatch(length(readLines(file.path(g, paste0(tag, "_ADMIX.bim")))), error = function(e) NA)
      rows[[length(rows) + 1]] <- data.frame(
        set = label, mode = ifelse(mode == "sup", "supervised", "unsupervised"), gen = gen, n = nrow(sim), snps = nsnp,
        true_mean = mean(sim$truth), est_mean = mean(sim$est_A),
        bias = mean(sim$est_A - sim$truth), RMSE = sqrt(mean((sim$est_A - sim$truth)^2)),
        r = suppressWarnings(cor(sim$est_A, sim$truth)), true_sd = sd(sim$truth),
        marker_RMSE = sqrt(mean((sim$A_frac_snp - sim$A_frac_bp)^2)),
        called_pure = mean(sim$est_A < 0.1 | sim$est_A > 0.9),
        refA_in_A = refA, refB_in_A = refB, refs_separated = refA > 0.5 & refB < 0.5,
        best_K_cv = if (nrow(cvm)) cvm$K[which.min(cvm$cv)] else NA)
    }
  }
}
if (!length(rows)) stop("No results found: are there .Q files in g*/unsup or g*/sup?")
summ <- bind_rows(rows) %>% arrange(set, mode, gen)
write.table(summ, paste0(opt$out, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(summ %>% select(set, mode, gen, snps, bias, RMSE, r, true_sd, called_pure, refs_separated, best_K_cv)), digits = 3)

theme_set(theme_bw() + theme(legend.position = "bottom"))
xb <- scale_x_continuous(breaks = sort(unique(summ$gen)))
truth_lab <- if (opt$basis == "bp") "whole-genome truth" else "SNP truth"
plots <- list()

long <- summ %>%
  transmute(set, mode, gen, `Bias (mean estimate - truth)` = bias, `RMSE per strain` = RMSE,
            `Correlation with truth (r)` = r, `Strains called pure (Q < 0.1 or > 0.9)` = called_pure) %>%
  pivot_longer(-c(set, mode, gen), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = unique(measure)))
plots$acc <- ggplot(long, aes(gen, value, colour = set, linetype = mode, shape = mode)) +
  geom_hline(data = data.frame(measure = factor("Bias (mean estimate - truth)", levels = levels(long$measure)), y = 0),
             aes(yintercept = y), colour = "grey60") +
  geom_line() + geom_point() + facet_wrap(~measure, scales = "free_y") + xb +
  labs(x = "Generation (rounds of intermating)", y = NULL, colour = "Scenario", linetype = "ADMIXTURE", shape = "ADMIXTURE",
       title = "ADMIXTURE accuracy (K = 2) by generation, each generation analysed separately",
       subtitle = paste0("Against ", truth_lab, "; generation 0 = founders"))

plots$noise <- ggplot(summ, aes(gen, RMSE, colour = set, linetype = mode, shape = mode)) +
  geom_line() + geom_point() +
  geom_line(aes(y = true_sd, linetype = NULL, shape = NULL), colour = "grey40", linetype = 3,
            data = summ %>% distinct(set, gen, true_sd)) +
  geom_line(aes(y = marker_RMSE, linetype = NULL, shape = NULL), colour = "#4DAC26", linetype = 2,
            data = summ %>% group_by(set, gen) %>% summarise(marker_RMSE = mean(marker_RMSE), .groups = "drop")) +
  facet_wrap(~set) + xb +
  labs(x = "Generation", y = "RMSE or SD of A fraction", colour = "Scenario", linetype = "ADMIXTURE", shape = "ADMIXTURE",
       title = "Estimation error vs how much strains truly differ",
       subtitle = "Dotted grey: true SD among strains (above it = worse than guessing the mean). Dashed green: error from the SNPs used alone.")

plots$means <- ggplot(summ, aes(gen, est_mean, colour = set, linetype = mode, shape = mode)) +
  geom_line() + geom_point() +
  geom_line(aes(y = true_mean, linetype = NULL, shape = NULL), colour = "black", data = summ %>% distinct(set, gen, true_mean)) +
  facet_wrap(~set) + xb +
  labs(x = "Generation", y = "Mean A fraction", colour = "Scenario", linetype = "ADMIXTURE", shape = "ADMIXTURE",
       title = "Mean ancestry of the sampled strains: truth (black) vs ADMIXTURE")

if (any(!is.na(summ$refB_in_A))) plots$refs <- ggplot(summ, aes(gen, refB_in_A, colour = set, linetype = mode, shape = mode)) +
  geom_hline(yintercept = 0.5, colour = "grey60") +
  geom_line() + geom_point() + xb + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Generation", y = "Mean Q of B references in the A cluster", colour = "Scenario",
       linetype = "ADMIXTURE", shape = "ADMIXTURE",
       title = "Do the K = 2 clusters separate the A and B references?",
       subtitle = "Near 0 = yes. Above 0.5 = the B references sit in the 'A' cluster, so the K = 2 split is not A vs B")

if (any(!is.na(summ$best_K_cv)))
  plots$k <- ggplot(summ %>% filter(!is.na(best_K_cv)), aes(gen, best_K_cv, colour = set)) +
    geom_line() + geom_point() + xb +
    labs(x = "Generation", y = "K with lowest CV error", colour = "Scenario",
         title = "K chosen by cross-validation (unsupervised)")

pdf(paste0(opt$out, ".pdf"), width = 11, height = 7)
for (p in plots) print(p)
invisible(dev.off())
message("Wrote ", opt$out, ".pdf and ", opt$out, ".tsv")

# ==============================================================================
# The full set of plot_admixture.R pages, for each scenario and mode, but with
# every generation analysed by its own ADMIXTURE runs: <out>_<set>_<mode>.pdf
# ==============================================================================
pal <- function(K) c("#2166AC", "#D6604D", "#4DAC26", "#8073AC", "#F1A340", "#35978F", "#C51B7D", "#878787",
                     "#B8E186", "#FDB863", "#5AAE61", "#E7D4E8")[seq_len(K)]
glab <- function(g) factor(paste("gen", g), levels = paste("gen", sort(unique(g))))
for (key in names(allq)) {
  qt <- allq[[key]]; rn <- allruns[[key]]; cm <- allcmp[[key]]
  ttl <- function(x) paste0(x, " (", key, "; each generation analysed separately)")
  gens <- sort(unique(qt$gen)); Ks <- sort(unique(qt$K))
  xg <- scale_x_continuous(breaks = gens)
  P <- list()

  # 1-2. cross-validation error and log-likelihood of every run, by generation
  if (any(!is.na(rn$cv))) {
    cvm <- rn %>% filter(!is.na(cv)) %>% group_by(gen, K) %>% summarise(cv = mean(cv), .groups = "drop")
    best <- cvm %>% group_by(gen) %>% slice_min(cv, n = 1)
    P$cv <- ggplot(rn %>% filter(!is.na(cv)), aes(K, cv)) +
      geom_line(data = cvm) + geom_jitter(width = 0.08, height = 0, alpha = 0.6, size = 1) +
      geom_point(data = best, colour = "firebrick", size = 3, shape = 1, stroke = 1.1) +
      facet_wrap(~glab(gen), scales = "free_y") + scale_x_continuous(breaks = Ks) +
      labs(y = "Cross-validation error", title = ttl("ADMIXTURE cross-validation"),
           subtitle = "Points = runs, line = mean; lowest mean CV error circled")
  }
  if (any(!is.na(rn$loglik)))
    P$loglik <- ggplot(rn %>% filter(!is.na(loglik)), aes(K, loglik, colour = best)) +
      geom_jitter(width = 0.08, height = 0, size = 1) +
      scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "firebrick"), labels = c("other runs", "best run (plotted)")) +
      facet_wrap(~glab(gen), scales = "free_y") + scale_x_continuous(breaks = Ks) +
      labs(y = "Log-likelihood", colour = NULL, title = ttl("Log-likelihood of each run"))

  # 3. admixture proportions: one column per generation, one row per K, truth on top
  ord <- qt %>% filter(K == min(Ks)) %>% distinct(gen, sample, kind, group, truth) %>%
    mutate(rank = case_when(group == "ref_A" ~ 1, group == "ref_B" ~ 2, kind == "sim" ~ 4, TRUE ~ 3)) %>%
    arrange(gen, rank, desc(truth)) %>% group_by(gen) %>% mutate(x = row_number()) %>% ungroup()
  tru <- ord %>% filter(!is.na(truth)) %>%
    transmute(gen, x, track = "Truth", C1 = truth, C2 = 1 - truth) %>%
    pivot_longer(c(C1, C2), names_to = "cluster", values_to = "q")
  bars <- bind_rows(tru, qt %>% left_join(ord %>% select(gen, sample, x), by = c("gen", "sample")) %>%
                      transmute(gen, x, track = paste0("K = ", K), cluster, q)) %>%
    mutate(track = factor(track, levels = c("Truth", paste0("K = ", Ks))),
           cluster = factor(cluster, levels = paste0("C", seq_len(max(Ks, 2)))))
  nref <- ord %>% group_by(gen) %>% summarise(xr = sum(rank < 4) + 0.5)
  P$bars <- ggplot(bars, aes(x, q, fill = cluster)) +
    geom_col(width = 1) +
    geom_vline(data = nref, aes(xintercept = xr), colour = "white", linewidth = 0.8) +
    facet_grid(track ~ glab(gen), scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal(max(Ks, 2)), labels = function(l) ifelse(l == "C1", "C1 (= A)", l)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Within each generation: A references, B references (left of the white line), then swarm strains sorted by true A fraction",
         y = "Ancestry proportion", fill = "Cluster", title = ttl("Admixture proportions (best run per K)"),
         subtitle = "Top row: true A fraction (blue = A). C1 = cluster of the A references, C2 = cluster of the B references (if any); other colours are not comparable between generations.") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = unit(1, "pt"),
          strip.text.y = element_text(angle = 0), strip.text.x = element_text(size = 7), axis.title.x = element_text(size = 9))

  # 4. cluster membership of the swarm strains
  memb <- qt %>% filter(kind == "sim") %>% group_by(gen, K, sample) %>%
    summarise(cluster = cluster[which.max(q)], qmax = max(q), .groups = "drop") %>%
    mutate(call = ifelse(qmax >= 0.9, "Q >= 0.9 (\"pure\")", "Q < 0.9 (admixed)"),
           K = factor(paste0("K = ", K), levels = paste0("K = ", Ks)),
           cluster = factor(cluster, levels = paste0("C", seq_len(max(Ks)))))
  P$membership <- ggplot(memb, aes(factor(gen), fill = cluster, alpha = call)) +
    geom_bar(colour = "grey30", linewidth = 0.2) + facet_wrap(~K) +
    scale_fill_manual(values = pal(max(Ks, 2)), drop = FALSE) +
    scale_alpha_manual(values = c(`Q < 0.9 (admixed)` = 0.4, `Q >= 0.9 ("pure")` = 1)) +
    labs(x = "Generation", y = "Swarm strains", fill = "Largest cluster", alpha = NULL,
         title = ttl("Cluster membership of the swarm strains (best run per K)"),
         subtitle = "Each strain assigned to its largest cluster; solid = Q >= 0.9 in that cluster")

  # 5-9. K = 2 against the truth
  if (nrow(cm) && "est_A" %in% names(cm)) {
    sim <- cm %>% filter(kind == "sim", !is.na(truth))
    refq <- cm %>% filter(kind == "ref") %>% group_by(gen, group) %>% summarise(C1 = mean(est_A), .groups = "drop")
    bad <- refq %>% filter(group == "ref_B", C1 > 0.5) %>% pull(gen)
    st <- sim %>% group_by(gen) %>%
      summarise(n = n(), r = suppressWarnings(cor(est_A, truth)), bias = mean(est_A - truth),
                RMSE = sqrt(mean((est_A - truth)^2)), .groups = "drop") %>%
      mutate(lab = sprintf("n = %d\nr = %.3f\nbias = %+.3f\nRMSE = %.3f%s", n, r, bias, RMSE,
                           ifelse(gen %in% bad, "\nclusters != A/B", "")))
    P$truth <- ggplot(sim, aes(truth, est_A)) +
      geom_abline(linetype = 2, colour = "grey50") + geom_point(alpha = 0.5, size = 1, colour = "#2166AC") +
      geom_text(data = st, aes(x = 0.02, y = 0.98, label = lab), hjust = 0, vjust = 1, size = 2.6, inherit.aes = FALSE) +
      facet_wrap(~glab(gen)) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = paste0("Truth (", truth_lab, ")"), y = "ADMIXTURE estimate of A (K = 2)",
           title = ttl("ADMIXTURE estimate vs truth"),
           subtitle = "Dashed line: perfect estimate. 'clusters != A/B': the B references sit in the A cluster in that generation")
    acc <- sim %>% group_by(gen) %>%
      summarise(`Bias (mean estimate - truth)` = mean(est_A - truth), `RMSE per strain` = sqrt(mean((est_A - truth)^2)),
                `Correlation with truth (r)` = suppressWarnings(cor(est_A, truth)),
                `SD of true A fraction among strains` = sd(truth), .groups = "drop") %>%
      pivot_longer(-gen, names_to = "measure", values_to = "value") %>%
      mutate(measure = factor(measure, levels = unique(measure)))
    P$by_gen <- ggplot(acc, aes(gen, value)) +
      geom_hline(data = data.frame(measure = factor("Bias (mean estimate - truth)", levels = levels(acc$measure)), y = 0),
                 aes(yintercept = y), linetype = 2, colour = "grey50") +
      geom_line() + geom_point() + facet_wrap(~measure, scales = "free_y") + xg +
      labs(x = "Generation", y = NULL, title = ttl("Accuracy of ADMIXTURE (K = 2) by generation"))
    mn <- sim %>% group_by(gen) %>% summarise(true = mean(truth), estimated = mean(est_A), .groups = "drop") %>%
      pivot_longer(-gen, names_to = "what", values_to = "A")
    P$means <- ggplot(mn, aes(gen, A, colour = what)) + geom_line() + geom_point() + xg +
      scale_colour_manual(values = c(true = "black", estimated = "#D6604D")) +
      labs(x = "Generation", y = "Mean A fraction of sampled strains", colour = NULL,
           title = ttl("Mean ancestry by generation: truth vs ADMIXTURE (K = 2)"))
    dec <- sim %>% group_by(gen) %>%
      summarise(`Markers: SNP truth vs genome truth` = sqrt(mean((A_frac_snp - A_frac_bp)^2)),
                `Tool: estimate vs SNP truth` = sqrt(mean((est_A - A_frac_snp)^2)),
                `Total: estimate vs genome truth` = sqrt(mean((est_A - A_frac_bp)^2)),
                `True spread among strains (SD, genome)` = sd(A_frac_bp), .groups = "drop") %>%
      pivot_longer(-gen, names_to = "what", values_to = "value")
    P$decomp <- ggplot(dec, aes(gen, value, colour = what)) + geom_line() + geom_point() + xg +
      scale_colour_manual(values = c("#4DAC26", "#D6604D", "black", "grey60")) +
      labs(x = "Generation", y = "RMSE (or SD) of A fraction", colour = NULL, title = ttl("Where the error comes from (K = 2)"),
           subtitle = "An error line above the grey line means the estimates rank strains worse than guessing the mean") +
      guides(colour = guide_legend(nrow = 2))
    if (nrow(refq)) P$refs <- ggplot(refq, aes(gen, C1, colour = group)) +
      geom_hline(yintercept = 0.5, colour = "grey60") + geom_line() + geom_point() + xg + coord_cartesian(ylim = c(0, 1)) +
      scale_colour_manual(values = c(ref_A = "#2166AC", ref_B = "#D6604D")) +
      labs(x = "Generation", y = "Mean Q in the A cluster (C1)", colour = NULL,
           title = ttl("Do the K = 2 clusters separate the A and B references?"),
           subtitle = "A references should be near 1 and B references near 0")
    P$error <- ggplot(sim, aes(truth, est_A - truth)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
      geom_point(alpha = 0.5, size = 1, colour = "#D6604D") + facet_wrap(~glab(gen)) +
      labs(x = paste0("Truth (", truth_lab, ")"), y = "Estimate - truth", title = ttl("Error in ADMIXTURE estimates (K = 2)"))
  }
  f <- paste0(opt$out, "_", gsub("[^A-Za-z0-9]+", "_", key), ".pdf")
  pdf(f, width = 13, height = 8)
  for (p in P) print(p)
  invisible(dev.off())
  message("Wrote ", f)
}
