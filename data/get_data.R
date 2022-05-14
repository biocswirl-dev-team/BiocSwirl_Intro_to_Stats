# PanNETs Gene Expression Dataset -----------------------------------------

# Homepage:
# https://www.nature.com/articles/s41467-018-06498-2

# Download RNA-seq gene expression data
pannets_expr_rnaseq <-
  paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118014&",
         "format=file&file=GSE118014_PanNETs_log2TPM_33_RSEM_STAR_",
         "Process_Samples.txt.gz") %>%
  read_tsv(col_types = cols()) %>%
  rename_with(toupper) %>%
  rename(Gene = GENE_ID)

# Download microarray gene expression data
pannets_expr_array <-
  paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE117851&",
         "format=file&file=GSE117851_Matrix_Gene_Expression_Collapsed_",
         "47_PanNETs_average.txt.gz") %>%
  read_tsv(col_types = cols()) %>%
  rename_with(toupper) %>%
  rename(Gene = X1) %>%
  select(Gene, one_of(colnames(pannets_expr_rnaseq)))

# Function for downloading and loading Excel supplemental data
download_supp_data <- function (link, ...) {
  temp <- tempfile()
  download.file(link, destfile = temp, quiet = TRUE)
  result <- read_excel(temp, na = c("", "none"), ...)
  unlink(temp)
  result
}

# Download clinical metadata
pannets_clinical <-
  paste0("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018",
         "-06498-2/MediaObjects/41467_2018_6498_MOESM4_ESM.xlsx") %>%
  download_supp_data() %>%
  transmute(Tumour = toupper(`Our ID`),
            Age, Sex = Gender,
            Genotype = Genotype...6,
            Metastasis = recode(`distant metastasis`,
                                "Y" = TRUE, "N" = FALSE),
            Subtype = ifelse(is.na(Genotype),
                             "A-D-M WT", "A-D-M Mutant")) %>%
  separate_rows(Genotype, sep = ", ") %>%
  pivot_wider(names_from = Genotype, values_from = Genotype,
              values_fn = list(Genotype = ~ TRUE),
              values_fill = list(Genotype = FALSE)) %>%
  rename_with(~ paste0(., "_mutant"), c(ATRX, DAXX, MEN1)) %>%
  select(-`NA`)

# Download ESTIMATE output
pannets_estimate <-
  paste0("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018",
         "-06498-2/MediaObjects/41467_2018_6498_MOESM6_ESM.xlsx") %>%
  download_supp_data(skip = 2) %>%
  transmute(Tumour = toupper(Description),
            Immune_score = round(ImmuneScore),
            Tumour_purity = round(`Tumor Purity`) / 100)

# Obtain RNA-seq and microarray gene expression for A-D-M genes
pannets_expr_adm <-
bind_rows(rnaseq = pannets_expr_rnaseq, array = pannets_expr_array,
          .id = "Method") %>%
  filter(Gene %in% c("ATRX", "DAXX", "MEN1")) %>%
  pivot_longer(cols = -c(Method, Gene), names_to = "Tumour",
               values_to = "Expr") %>%
  pivot_wider(id_cols = Tumour, names_from = c(Gene, Method),
              names_sep = "_", values_from = Expr)

# Merge all metadata and data together
pannets_metadata <-
  pannets_clinical %>%
  left_join(pannets_estimate, by = "Tumour") %>%
  left_join(pannets_expr_adm, by = "Tumour") %>%
  relocate(Tumour_purity, Immune_score, .after = Metastasis)

# Output data and metadata
write_csv(pannets_expr_rnaseq, "data/pannets_expr_rnaseq.csv.gz")
write_csv(pannets_expr_array, "data/pannets_expr_array.csv.gz")
write_csv(pannets_metadata, "data/pannets_metadata.csv")
