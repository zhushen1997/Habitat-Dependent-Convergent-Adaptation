###======================================================================================================================================###
# step 1: data clean
###======================================================================================================================================###


rm(list=ls())
library(pmartR)
library(impute)
library(mtExtra)
library(missForest)
library(parallel)


### 1. load data
peak_df <- read.csv("1_metabolome_matrix/1_metabolome_matrix.csv",  check.names = FALSE)
meta_df <- read.csv("1_metabolome_matrix/2_group.csv", check.names = FALSE)
meta_df <- meta_df[order(match(meta_df$Sample, colnames(peak_df))), ]
meta_df$Sample <- as.character(meta_df$Sample)


### 2. build pmartR subject
# The data matrix is e_data, the sample information is f_data, and it will automatically replace 0 with NA.
peak_se <- as.metabData(e_data = peak_df,
                        f_data = meta_df,
                        edata_cname = "CPD_ID",     # The column name of metabolite ID in e_data
                        fdata_cname = "Sample"      # The column name of sample ID in the f_data table
)


### 3. data clean
## 3.1 sample filter 
# Eliminate the samples with an excessively high proportion of missing data
sample_missing_ratio <- colSums(is.na(peak_se$e_data[, -1])) / nrow(peak_se$e_data)
bad_samples <- names(which(sample_missing_ratio > 0.5))

if (length(bad_samples) > 0) {
  filter_bad_samples <- custom_filter(peak_se, f_data_remove = bad_samples)
  peak_se_filtered <- applyFilt(filter_bad_samples, peak_se)
} else {
  peak_se_filtered <- peak_se
}

## 3.2 metabolites filter
# Remove the metabolites with an excessively high proportion of missing values
moleculeFilt <- molecule_filter(omicsData = peak_se_filtered)
summary(moleculeFilt, min_num = 6)
plot(moleculeFilt, min_num = 6)
peak_se_filtered <- applyFilt(filter_object = moleculeFilt, omicsData = peak_se_filtered, min_num = 6)


## 3.3 QC-RSD filter(CV filtering based on QC samples)
# Based on the filtering of QC samples using the relative standard deviation (RSD)
# Eliminating the metabolites with significant fluctuations in the instrument or pre-treatment process.
sam_class <- ifelse(meta_df$group == "QC", "qc", "sample")
# get log2(K+1) matrix
mat <- peak_se_filtered$e_data %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "CPD_ID") %>%
  as.matrix() %>%
  {. + 1} %>%
  log2() %>%
  t()

# RSD filtering based on QC samples
qc_rsd_filter <- mtExtra::qc_filter(
  x         = mat,
  y         = sam_class,      # An vector containing elements such as "sample", "qc" and "blank"
  thres_rsd = 30,             # QC-RSD threshold 30 %
  f_mv      = FALSE,          # Used to indicate whether to perform the missing value filtering operation on "sample" data or "qc" data
  f_mv_qc_sam = FALSE,        # The flag for filtering based on the percentage of missing values in "qc" or "sample", where TRUE indicates "qc" and FALSE indicates "sample"
  thres_mv  = 0.3)            # If the missing rate is greater than 30%, it will also be excluded.

table(qc_rsd_filter$idx)

# Reconvert the matrix filtered by QC-RSD back into a pmartR object
filtered_df <- as.data.frame(t(qc_rsd_filter$dat))
filtered_df$CPD_ID <- rownames(filtered_df)
filtered_df <- as.data.frame(filtered_df) %>%
  tibble::remove_rownames() %>%
  dplyr::relocate(CPD_ID, .before = everything())


peak_se_qcRSD <- as.metabData(e_data = filtered_df,
                              f_data = meta_df,
                              edata_cname = "CPD_ID",
                              fdata_cname = "Sample"
)


## 3.4 CV Filtering (CV Filtering Based on Sample Grouping) (Optional step, but sometimes the differences between samples are quite significant. In such cases, you can appropriately increase the cv_threshold or refrain from using CV filtering based on sample grouping)
# Group by group, remove the metabolites with large fluctuations within the group, and retain the metabolites that "may truly have biological significance"
cv_threshold <- 50
# Calculate CV
cvfilter <- cv_filter(omicsData = peak_se_qcRSD, use_groups = TRUE)
# Check the minimum and maximum CV values
min_cv <- min(cvfilter$CV)
max_cv <- max(cvfilter$CV)
# filtering
if (max_cv > cv_threshold) {
  peak_se_CV <- applyFilt(filter_object = cvfilter, omicsData = peak_se_qcRSD, cv_threshold = cv_threshold)
} else {
  peak_se_CV <- peak_se_qcRSD
}

# Restore the original scale matrix
keep_idx <- row.names(peak_se_CV$e_data)
CV_raw <- peak_se$e_data[keep_idx, ]
peak_se_CV_raw <- as.metabData(e_data = CV_raw,
                               f_data = meta_df,
                               edata_cname = "CPD_ID",
                               fdata_cname = "Sample"
)


### 4. Normalization (LOEES / PQN)
# LOESS mainly addresses the signal drift of a single feature over time (such as due to fluctuations in instrument performance)
# PQN mainly addresses the overall concentration differences among samples (such as those caused by different dilution levels)

# log2 transfer
peak_se_log_final <- edata_transform(
  omicsData  = peak_se_CV_raw,
  data_scale = "log2"           
)

# LOEES Normalization
loees_factors <-  normalize_loess(
  omicsData = peak_se_log_final,     # log2 transfer matrix
  method    = "fast",                # "fast" (default), "affy", or "pairs"
  span      = 0.4                    # Smoothness
)

# PQN Normalization
pqn_factors <- normalize_global(
  omicsData     = loees_factors,
  subset_fn     = "all",
  norm_fn       = "median",
  apply_norm    = TRUE,
  backtransform = TRUE
)


#### 5. Fill in the missing values
### In our reasearch, missForest was choosed to fill in the missing value

### 5.1 KNN
knn_raw <- pqn_factors$e_data %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "CPD_ID") %>%
  as.matrix()

knn_imp <- impute.knn(knn_raw, k = 5, rowmax = 0.5)$data %>%   
  as.data.frame() %>%
  tibble::rownames_to_column(var = "ID") %>%
  tibble::remove_rownames()

any(is.na(knn_imp))
any(knn_imp == 0)

knn_imp_noQC <- knn_imp %>% dplyr::select(-c(26:29))

### 5.2 missForest
# The columns correspond to the variables and the rows to the observations;
expr <- pqn_factors$e_data %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "CPD_ID") %>%
  as.matrix()
expr[expr == 0] <- NA

# multi-threading
cl <- makeCluster(detectCores() - 1)
setDefaultCluster(cl)

set.seed(1234)
imp  <- missForest(expr,
                   ntree   = 200,
                   mtry    = sqrt(ncol(expr)),
                   maxiter = 10,
                   verbose = TRUE)

expr_imp <- imp$ximp %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "CPD_ID") %>%
  tibble::remove_rownames()

imp$OOBerror
sum(is.na(expr_imp))

expr_imp_noQC <- expr_imp %>% dplyr::select(-c(26:29))

### 7. save results
write.table(knn_imp, file = "1_metabolome_matrix/3_metabolome_LOESS_PQN_KNN_matrix.tsv", sep = "\t", quote = F, row.names = FALSE)
write.table(knn_imp_noQC, file = "1_metabolome_matrix/4_metabolome_LOESS_PQN_KNN_matrix_noQC.tsv", sep = "\t", quote = F, row.names = FALSE)

write.table(expr_imp, file = "1_metabolome_matrix/5_metabolome_LOESS_PQN_missForest_matrix.tsv", sep = "\t", quote = F, row.names = FALSE)
write.table(expr_imp_noQC, file = "1_metabolome_matrix/6_metabolome_LOESS_PQN_missForest_matrix_noQC.tsv", sep = "\t", quote = F, row.names = FALSE)