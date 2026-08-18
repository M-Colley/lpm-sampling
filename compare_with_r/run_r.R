# Draw the same number of LPM samples with the reference R implementation.
#
#   Rscript run_r.R [replicates] [method]
#
# Requires:  install.packages("BalancedSampling")
# BalancedSampling is by the authors of the method (Grafstrom et al.) and is
# AGPL-3 licensed. It is used here only as an external reference to compare
# against; no code from it is copied into this package.

args <- commandArgs(trailingOnly = TRUE)
replicates <- if (length(args) >= 1) as.integer(args[1]) else 5000
method <- if (length(args) >= 2) args[2] else "lpm2"

suppressMessages(library(BalancedSampling))

frame_path <- file.path("frames", "frame.csv")
if (!file.exists(frame_path)) stop("run frames.py first to write frames/frame.csv")

frame <- read.csv(frame_path)
coords <- as.matrix(frame[, c("x", "y")])
pi <- frame$pi
N <- nrow(coords)

counts <- numeric(N)
balances <- numeric(replicates)

set.seed(1)
started <- Sys.time()

for (r in seq_len(replicates)) {
  s <- if (method == "lpm1") lpm1(pi, coords) else lpm2(pi, coords)
  counts[s] <- counts[s] + 1
  # sb() expects the same probabilities, coordinates and sample indices
  balances[r] <- sb(pi, coords, s)
}

elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))

dir.create("results", showWarnings = FALSE)
write.csv(data.frame(frequency = counts / replicates),
          file.path("results", paste0("r_", method, "_freq.csv")), row.names = FALSE)
write.csv(data.frame(balance = balances),
          file.path("results", paste0("r_", method, "_balance.csv")), row.names = FALSE)

cat(sprintf("R %s: %d replicates in %.1f s (%.2f ms per draw)\n",
            method, replicates, elapsed, 1000 * elapsed / replicates))
cat(sprintf("mean spatial balance %.5f\n", mean(balances)))
