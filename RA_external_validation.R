 ===============================================================
# PART 1: GENUS-LEVEL TRAINING SCRIPT
# Train XGBoost model on Li et al. dataset aggregated to genus
# Save model + training genus feature names
# ===============================================================



library(dplyr)
library(stringr)
library(xgboost)

# -------------------------------------------------
# Load files
# -------------------------------------------------
asv_profile <- readRDS("1.ASV.profile.rds")
tax_info    <- readRDS("1.taxonomy.info.rds")

# -------------------------------------------------
# Convert ASV table
# rows = samples, cols = ASVs
# -------------------------------------------------
asv_mat <- as.data.frame(asv_profile)
asv_mat <- t(asv_mat)
asv_mat <- as.data.frame(asv_mat)

# labels
asv_mat$Group <- ifelse(grepl("^RA", rownames(asv_mat)), "RA", "HC")

# -------------------------------------------------
# Correct taxonomy table
# -------------------------------------------------
tax_df <- data.frame(
  ASV = rownames(tax_info),
  Taxon = tax_info$Taxon,
  stringsAsFactors = FALSE
)

# extract genus
tax_df$Genus <- str_extract(tax_df$Taxon, "g__[^;]+")

tax_df$Genus[is.na(tax_df$Genus)] <- "Unknown"
tax_df$Genus <- sub("g__", "", tax_df$Genus)

# -------------------------------------------------
# Keep ASVs present in matrix
# -------------------------------------------------
X <- asv_mat[, !colnames(asv_mat) %in% "Group"]

tax_df <- tax_df[match(colnames(X), tax_df$ASV), ]
tax_df$Genus[is.na(tax_df$Genus)] <- "Unknown"

# -------------------------------------------------
# Aggregate to genus
# -------------------------------------------------
X_t <- t(data.matrix(X))

genus_table <- rowsum(X_t, group = tax_df$Genus)

genus_df <- as.data.frame(t(genus_table), check.names = FALSE)

# add labels
genus_df$Group <- asv_mat$Group

# -------------------------------------------------
# Check structure
# -------------------------------------------------
dim(genus_df)
length(colnames(genus_df))
head(colnames(genus_df),20)

# -------------------------------------------------
# Train model
# -------------------------------------------------
train_x <- data.matrix(genus_df[, !colnames(genus_df) %in% "Group"])
train_y <- ifelse(genus_df$Group == "RA", 1, 0)

dtrain <- xgb.DMatrix(data=train_x,label=train_y)

params <- list(
  objective="binary:logistic",
  eval_metric="auc",
  max_depth=4,
  eta=0.05,
  subsample=0.8,
  colsample_bytree=0.8
)

model <- xgb.train(
  params=params,
  data=dtrain,
  nrounds=300
)

cat("Saved genus model.\n")
cat("Features:", length(features), "\n")

# ==========================================================
# Part 2: GENUS LEVEL EXTERNAL VALIDATION
# Uses trained genus XGBoost model
# ==========================================================

library(dada2)
library(dplyr)
library(stringr)
library(tibble)
library(pROC)
library(xgboost)

# ----------------------------------------------------------
# PATH
# ----------------------------------------------------------
setwd(path)

# ----------------------------------------------------------
# FASTQ FILES
# ----------------------------------------------------------
unlink(file.path(path,"filtered"), recursive=TRUE)
dir.create(file.path(path,"filtered"))

# forward files only
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz$", full.names = TRUE))

# extract sample names
sample.names <- sub("_1.fastq.gz", "", basename(fnFs))

# create reverse files in same order
fnRs <- file.path(path, paste0(sample.names, "_2.fastq.gz"))

fwd <- sub("_1.fastq.gz","",list.files(path, pattern="_1.fastq.gz$"))
rev <- sub("_2.fastq.gz","",list.files(path, pattern="_2.fastq.gz$"))

setdiff(fwd, rev)
setdiff(rev, fwd)

# check missing reverse files
fnRs[file.exists(fnRs) == FALSE]

length(fnFs)
length(fnRs)

# ----------------------------------------------------------
# FILTERING
# ----------------------------------------------------------
filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(240,200),
  maxN = 0,
  maxEE = c(5,5),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE
)

# ----------------------------------------------------------
# LEARN ERRORS
# ----------------------------------------------------------
errF <- learnErrors(filtFs, nbases = 1e7)
errR <- learnErrors(filtRs, nbases = 1e7)

# ----------------------------------------------------------
# DEREP + DADA
# ----------------------------------------------------------
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF)
dadaRs <- dada(derepRs, err = errR)

# ----------------------------------------------------------
# MERGE
# ----------------------------------------------------------
mergers <- mergePairs(
  dadaFs, derepFs,
  dadaRs, derepRs,
  minOverlap = 8,
  maxMismatch = 2
)

# ----------------------------------------------------------
# SEQUENCE TABLE
# ----------------------------------------------------------
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus")

# ----------------------------------------------------------
# TAXONOMY
# ----------------------------------------------------------
taxa <- assignTaxonomy(
  seqtab.nochim,
  "silva_nr99_v138.1_train_set.fa.gz"
)

taxa_df <- as.data.frame(taxa)
taxa_df <- taxa_df[colnames(seqtab.nochim), ]

# ----------------------------------------------------------
# GENUS AGGREGATION
# ----------------------------------------------------------
genus <- taxa_df$Genus
genus[is.na(genus)] <- "Unknown"

asv_counts <- as.data.frame(seqtab.nochim)

genus_table <- rowsum(t(asv_counts), group = genus)

genus_table_t <- as.data.frame(t(genus_table))

# ----------------------------------------------------------
# TRUE LABELS
# EDIT THIS TO MATCH YOUR SAMPLE ORDER
# ----------------------------------------------------------
genus_table_t$Group <- c(
  "HC","HC","HC","HC","HC","HC","RA","HC","HC","RA",
  "RA","RA","RA","RA","RA","HC","RA","HC","RA","RA",
  "HC","HC","HC","HC","HC","HC","HC","HC","HC","HC",
  "RA","RA","RA","RA","RA","RA","RA","RA","RA","RA"
)


# ----------------------------------------------------------
# LOAD TRAINED MODEL
# ----------------------------------------------------------
model <- xgb.load("model_xgb_genus.model")
train_features <- readRDS("train_genus_features.rds")

# ----------------------------------------------------------
# ALIGN FEATURES
# ----------------------------------------------------------
val_df <- genus_table_t[, !colnames(genus_table_t) %in% "Group"]

common <- intersect(train_features, colnames(val_df))
val_df <- val_df[, common, drop = FALSE]

missing <- setdiff(train_features, colnames(val_df))

for(g in missing){
  val_df[[g]] <- 0
}

val_df <- val_df[, train_features]

# ----------------------------------------------------------
# PREDICT
# ----------------------------------------------------------

val_df[] <- lapply(val_df, function(x) as.numeric(as.character(x)))

val_matrix <- data.matrix(val_df)

pred_prob <- predict(model, val_matrix)


# ----------------------------------------------------------
# ROC
# ----------------------------------------------------------
true_labels <- factor(genus_table_t$Group,
                      levels = c("HC","RA"))

roc_obj <- roc(true_labels, pred_prob)

print(auc(roc_obj))

# ----------------------------------------------------------
# PLOT ROC
# ----------------------------------------------------------
plot(
  roc_obj,
  col = "blue",
  lwd = 3,
  main = paste("External Validation ROC (AUC =", round(auc(roc_obj),3),")")
)

abline(a=0,b=1,lty=2,col="gray")

# Figure 12 - Save External Validation ROC Curve

png(
  filename = "Figure12_External_Validation_ROC.png",
  width = 2400,
  height = 2000,
  res = 300
)

plot(
  roc_obj,
  col = "blue",
  lwd = 3,
  main = paste(
    "External Validation ROC (AUC =",
    round(auc(roc_obj), 3),
    ")"
  ),
  cex.main = 1.4
)

abline(a = 0, b = 1, lty = 2, col = "gray")

dev.off()
# ----------------------------------------------------------
# CONFUSION MATRIX
# ----------------------------------------------------------
pred_class <- ifelse(pred_prob > 0.5, "RA", "HC")

print(table(Predicted = pred_class,
            Actual = genus_table_t$Group))

summary(pred_prob)
pred_prob

length(common)
length(train_features)

data.frame(Sample=rownames(genus_table_t), Label=genus_table_t$Group)

table(genus_table_t$Group) 

ci.auc(roc_obj)
coords(roc_obj, "best", ret=c("threshold","sensitivity","specificity")) 
###############################################################
#checking authenticity 
data.frame(
  Sample = rownames(genus_table_t),
  Prob = pred_prob,
  Label = genus_table_t$Group
)

tapply(pred_prob, genus_table_t$Group, summary)
