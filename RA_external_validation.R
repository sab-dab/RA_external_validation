# ============================================================
# RA External Validation Script (Genus-level)
# Author: Sabira Dabeer
# ============================================================

# -----------------------------
# 1. Load libraries
# -----------------------------
library(dada2)
library(dplyr)
library(stringr)
library(tibble)
library(pROC)
library(xgboost)

# -----------------------------
# 2. Set path
# -----------------------------
path <- "C:/Users/sabir/RA_validation"

# -----------------------------
# 3. Load filtered files
# -----------------------------
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))

sample.names <- basename(fnFs) %>%
  str_remove("_1.fastq.gz")

# -----------------------------
# 4. Filtering
# -----------------------------
filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(240, 200),
  maxN = 0,
  maxEE = c(5,5),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE
)

print(out)

# -----------------------------
# 5. Learn errors (fast mode)
# -----------------------------
errF <- learnErrors(filtFs, nbases = 1e7)
errR <- learnErrors(filtRs, nbases = 1e7)

# -----------------------------
# 6. Dereplication
# -----------------------------
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

# -----------------------------
# 7. DADA
# -----------------------------
dadaFs <- dada(derepFs, err = errF)
dadaRs <- dada(derepRs, err = errR)

# -----------------------------
# 8. Merge
# -----------------------------
mergers <- mergePairs(
  dadaFs, derepFs,
  dadaRs, derepRs,
  minOverlap = 8,
  maxMismatch = 2
)

# -----------------------------
# 9. Sequence table
# -----------------------------
seqtab <- makeSequenceTable(mergers)
print(dim(seqtab))

# -----------------------------
# 10. Remove chimera
# -----------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus")
print(dim(seqtab.nochim))

# -----------------------------
# 11. Taxonomy assignment
# -----------------------------
taxa <- assignTaxonomy(
  seqtab.nochim,
  "silva_nr99_v138.1_train_set.fa.gz"
)

taxa_df <- as.data.frame(taxa)

# -----------------------------
# 12. Genus aggregation
# -----------------------------
taxa_df <- taxa_df[colnames(seqtab.nochim), ]

genus <- taxa_df$Genus
genus[is.na(genus)] <- "Unknown"

asv_counts <- as.data.frame(seqtab.nochim)

genus_table <- rowsum(t(asv_counts), group = genus)

genus_table_t <- t(genus_table)
genus_table_t <- as.data.frame(genus_table_t)

# -----------------------------
# 13. Add labels (EDIT if needed)
# -----------------------------
genus_table_t$Group <- c(
  "HC","HC","HC","HC","HC","HC","RA","HC","HC","RA",
  "RA","RA","RA","RA","RA","HC","RA","HC","RA","RA"
)

# -----------------------------
# 14. Load trained model
# -----------------------------
model <- xgb.load("model_xgb.model")
train_features <- readRDS("train_asvs.rds")

# -----------------------------
# 15. Align features
# -----------------------------
val_df <- genus_table_t
val_df <- val_df[, !colnames(val_df) %in% c("Group")]

common <- intersect(train_features, colnames(val_df))
val_df <- val_df[, common, drop = FALSE]

missing <- setdiff(train_features, colnames(val_df))
for (g in missing) {
  val_df[[g]] <- 0
}

val_df <- val_df[, train_features]

# -----------------------------
# 16. Predict
# -----------------------------
val_matrix <- as.matrix(val_df)
pred_prob <- predict(model, val_matrix)

# -----------------------------
# 17. Evaluate
# -----------------------------
true_labels <- genus_table_t$Group

roc_obj <- roc(true_labels, pred_prob)
print(auc(roc_obj))

pred_class <- ifelse(pred_prob > 0.5, "RA", "HC")
print(table(pred_class, true_labels))
