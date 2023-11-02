# Introduction ------------------------------------------------------------
# This script builds config files for kGWASflow.

library(dplyr)
library(googlesheets4)
library(stringr)
library(tidyr)

# Adding my Google service account credentials
gs4_auth(path = "~/.credentials/google_sheets_api/service_account.json")

# Test run: Varitome and flower measurements ------------------------------
# Here I will test kGWASflow with the anther / pistil ratio measurements and 
# the Varitome accessions. I'd like to know how the pipeline functions and how 
# it scales.
flower_phenotypes <- read_sheet("1Jyv0pj_SpE_9evfOt0Eb2nGFJN_Kfg4kMC3l1sJLsJU") %>%
  group_by(accession_id) %>%
  summarize(mean_anther_pistil_ratio = mean(anther_over_pistil)) %>%
  rename(accession = accession_id)

varitome_metadata <- read.table(
    file = file.path(getwd(), "R_data", "varitome_metadata.txt"),
    sep = ',',
    header = TRUE
  ) %>%
  select(Run, Sample.Name)

varitome_paths <- read.table(
    file = file.path(getwd(), "R_data", "varitome_paths.txt"),
    sep = ',',
    header = TRUE
  ) %>%
  mutate(
    library_name = str_extract(path, "SRR\\d+"),
    fq_read = case_when(
      str_detect(path, "_1\\.fastq\\.gz$") ~ "fq1",
      str_detect(path, "_2\\.fastq\\.gz$") ~ "fq2",
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_wider(
    names_from = fq_read,
    values_from = path
  )

cedar_metadata <- read_sheet("1V2kH8G4tfYsYqnYb6bVHpGjwwkjd9le0arGBJ2o4r8s") %>%
  select(name_CW, name_original_razifard)

# Joining in the CW identifiers
varitome_metadata <- varitome_metadata %>%
  left_join(cedar_metadata, join_by(Sample.Name == name_original_razifard)) %>%
  drop_na() %>%
  rename(accession = name_CW)

# There are multiple libraries for some of the accessions, figure out how to deal 
# with these. I'll choose a random one for the test.
qc <- varitome_metadata %>%
  group_by(accession) %>%
  summarize(n = n())

set.seed(13)
varitome_metadata <- varitome_metadata %>%
  group_by(accession) %>%
  slice_sample(n = 1)

# Joining with the phenotype data
flower_metadata <- flower_phenotypes %>%
  left_join(varitome_metadata, by = "accession") %>%
  drop_na()

# Now I need to get it in the right format, including the full fastq paths.
samplesheet_flowers <- flower_metadata %>%
  select(accession, Run) %>%
  rename(sample_name = accession, library_name = Run) %>%
  left_join(varitome_paths, by = "library_name") %>%
  mutate(sra = library_name)

# Writing out the samplesheet
samplesheet_flowers %>%
  write.table(
    file = file.path(getwd(), "kgwasflow", "config", "flowers_and_varitome_test", "samples.tsv"),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  )

# Making and saving the phenos.tsv file
phenos_flowers <- data.frame(
  pheno_name = "anther_pistil_ratio", 
  pheno_path = "/home/u16/cedar/git/pollen_quantitative_genetics/kgwasflow/config/flowers_and_varitome_test/anther_pistil_ratio.pheno"
)

phenos_flowers %>%
  write.table(
    file = file.path(getwd(), "kgwasflow", "config", "flowers_and_varitome_test", "phenos.tsv"),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  )

# Making and saving the anther / pistil ratio phenotype file.
phenotype_file_anther_pistil <- flower_phenotypes %>%
  rename(accession_id = accession, phenotype_value = mean_anther_pistil_ratio) %>%
  filter(accession_id %in% samplesheet_flowers$sample_name) # Only keeping the Varitome accessions

phenotype_file_anther_pistil %>%
  write.table(
    file = file.path(getwd(), "kgwasflow", "config", "flowers_and_varitome_test", "anther_pistil_ratio.pheno"),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  )
