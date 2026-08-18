# R half of the deterministic checks (see diagnostic.py).
#
#   Rscript diagnostic.R
#
# Scores the fixed inputs written by diagnostic.py with the reference
# implementation, so compare.py can check the two sides agree exactly:
#   (a) BalancedSampling::sb       vs  lpm_sampling.spatial_balance
#   (b) BalancedSampling::getPips  vs  lpm_sampling.pi_from_size
#
# Samples arrive as 0/1 masks, so the 0-based / 1-based index convention gap
# between Python and R cannot silently shift a sample.

suppressMessages(library(BalancedSampling))

frame_path <- file.path("frames", "frame.csv")
mask_path <- file.path("frames", "diag_masks.csv")
if (!file.exists(mask_path)) stop("run diagnostic.py first to write frames/diag_masks.csv")

frame <- read.csv(frame_path)
coords <- as.matrix(frame[, c("x", "y")])
pi <- frame$pi

dir.create("results", showWarnings = FALSE)

# --- (a) spatial balance ----------------------------------------------------
masks <- as.matrix(read.csv(mask_path, header = FALSE))
balance <- numeric(nrow(masks))

for (r in seq_len(nrow(masks))) {
  balance[r] <- sb(pi, coords, which(masks[r, ] == 1))
}

write.csv(data.frame(balance = balance),
          file.path("results", "r_diag_sb.csv"), row.names = FALSE)
cat(sprintf("balance: %d cases (range %.5f to %.5f)\n",
            length(balance), min(balance), max(balance)))

# --- (b) pps probabilities --------------------------------------------------
sizes <- as.matrix(read.csv(file.path("frames", "diag_pps_sizes.csv"), header = FALSE))
targets <- read.csv(file.path("frames", "diag_pps_n.csv"), header = FALSE)[[1]]

probabilities <- matrix(0.0, nrow = nrow(sizes), ncol = ncol(sizes))
for (k in seq_len(ncol(sizes))) {
  probabilities[, k] <- getPips(sizes[, k], targets[k])
}

# Header row on purpose -- compare.py skips exactly one line in every results
# file, so a headerless matrix would silently lose its first unit.
colnames(probabilities) <- paste0("case", seq_len(ncol(probabilities)) - 1L)
write.csv(probabilities, file.path("results", "r_diag_pps.csv"), row.names = FALSE)
cat(sprintf("pps:     %d cases over %d units\n", ncol(sizes), nrow(sizes)))
