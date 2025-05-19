#Making input tsv files for running the kmers gwas nextflow pipeline with both long and short read genomes included

# Creating the master sheet------------------------------------------------------------------

# loading necessary libraries
library(dplyr)
library(stringr)
library(tidyr)


# Set working directory
setwd('/Users/cperkins/Desktop/palanivelu_lab/kmers_gwas/')

# building a sheet that contains everything so that it's all aligned, and then I will pull the necessary columns from this sheet for the .tsv files I need for kgwasflow

# reading in accessions sheet containing the names for each accession, including CW
accessions <- read.csv('data/accessions.csv') %>%
  select(name_CW:name_original_razifard)
# 
# # reading in varitome metadata sheet containing all the seq info for the varitome div panel, including sras & library
varitome_metadata <- read.csv('data/solavar_metadata.txt')
var_sra <- varitome_metadata %>% # pulling out necessary columns, SRA#, sample name, SRX#
  select(Run, Sample.Name, Experiment) %>%
  mutate(Sample.Name = str_replace_all(Sample.Name, ' ', '')) #some of the samples have an unnecessary space

alonge_metadata <- read.csv('data/alonge_metadata.txt')

# alonge has a mix of short and long read sequencing. We only want to keep the long reads (duplicate sequencing methods for some accessions)
# Filter for rows where LibraryLayout is "SINGLE"
alonge_metadata_single <- alonge_metadata %>%
  filter(LibraryLayout == "SINGLE")

alonge_sra <- alonge_metadata_single %>%
select(Run, Sample.Name, Experiment)
# 
# #joining sras with cw
# 
# var_sheet <- merge(var_sra, accessions, by.x = "Sample.Name", by.y = "name_original_razifard", suffixes = c("_var_sra", "_accessions"))
# 
# 
# 
# alonge_sheet <- merge(alonge_sra, accessions, by.x = "Sample.Name", by.y = "packet_name_1", suffixes = "_alonge_sra", "_accessions")

# 
var_sheet <- var_sra %>%
 left_join(accessions, join_by(Sample.Name == name_original_razifard))

var_sheet <- var_sheet %>%
 arrange(desc(Run))

alonge_sheet <- alonge_sra %>%
  left_join(accessions, join_by(Sample.Name == packet_name_1))

alonge_sheet <- alonge_sheet %>%
  arrange(desc(Run))

# Add read_type to each sheet
var_sheet$read_type <- "paired"
alonge_sheet$read_type <- "single"

# Join the two sheets
dp_sheet <- full_join(var_sheet, alonge_sheet)

#finding the accessions that are not present in either varitome or alonge datasets
missing_in_dp_sheet <- setdiff(accessions$name_CW, dp_sheet$name_CW)
print(missing_in_dp_sheet)

# Get names in accessions that are missing in dp_sheet
missing_in_dp_sheet <- setdiff(accessions$name_CW, dp_sheet$name_CW)

# Extract the rows from accessions
accessions_missing <- accessions[accessions$name_CW %in% missing_in_dp_sheet, ]
print(accessions_missing$packet_name_1)

# There are multiple libraries for some of the accessions, figure out how to deal
# with these. I'll choose a random one for the test.
qc <- dp_sheet %>%
  group_by(name_CW) %>%
  summarize(n = n())

set.seed(13)
dp_sheet <- dp_sheet %>%
  group_by(name_CW) %>%
  slice_sample(n = 1)

#Keeping only necessary columns
dp_sheet <- dp_sheet %>%
  select(Run, Sample.Name, Experiment, name_CW, read_type) %>%
  drop_na()

# adding missing accesssions

# Create the manual data frame
manual_entries <- tribble(
  ~Run, ~Sample.Name, ~Experiment, ~name_CW, ~read_type,
  "SRR14191280", "LA4345", "SRX10558138", "CW0000", "long",
  "SRR29823272", "LA0490", "SRX25322143", "CW0001", "long",
  "SRR29540650", "LA1994", "SRX25048037", "CW0002", "paired",
  "SRR31986498", "LA2661", "SRX27341537", "CW0003", "paired",
  "SRR941559",   "LA2375", "SRX326415",   "CW0006", "paired",
  "SRR31986664", "LA2662", "SRX27341371", "CW0015", "paired",
  "SRR10208149", "LA3317", "SRX6928073",  "CW0030", "long",
  "SRR11093182", "PI647486", "SRX7732019", "CW0044", "long",
  "SRR10214511", "LA4026", "SRX6934191",  "CW0066", "long",
  "SRR10214512", "LA3840", "SRX6934190",  "CW0070", "long",
  "SRR10199004", "LA0502", "SRX6919144",  "CW0092", "long",
  "SRR10199003", "LA3465", "SRX6919145",  "CW0110", "long",
  "SRR10208148", "LA3242", "SRX6928074",  "CW0127", "long",
  "SRR9856933",  "LA3475", "SRX6610777",  "CW0128", "long"
)

# Replace 'long' with 'single'
manual_entries <- manual_entries %>%
  mutate(read_type = ifelse(read_type == "long", "single", read_type))

dp_sheet_full <- full_join(dp_sheet, manual_entries)

write.csv(dp_sheet_full, file = file.path(getwd(), "data", "dp_sheet.csv"), row.names = FALSE)


# fetchngs input tsvs--------------------------------------------------------------------------------------------------

#with dp_sheet created, I will now use it to build all my other inputs.
#First input will be for fetchngs

#read in constructed sheet from data

# loading necessary libraries
library(dplyr)
library(stringr)
library(tidyr)


# Set working directory
setwd('/Users/cperkins/Desktop/palanivelu_lab/kmers_gwas/')

#dp_sheet has all necessary info for pulling all the SRR IDs that match exactly with the SRRs we have for our accessions.
#Since there are multiple SRRs for each accession, it is critical that the SRR we chose in dp_sheet matches the rest of the inputs

# Load the pre-created dp_sheet
dp_sheet <- read.csv(file.path(getwd(), "data", "dp_sheet.csv"), row.names = NULL)

#Now just extracting the SRRs

# Extract just the Run column
srrs <- data.frame(Run = dp_sheet$Run) %>%
  arrange(desc(Run))

# Write to a .tsv file
write.table(srrs, file = file.path(getwd(), "nextflow", "fetchngs", "input","dp_srrs.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


#also pull out subsets for the test runs

# Short reads:
short_reads <- c("SRR7279519", "SRR7279511", "SRR7279567", "SRR7279572", "SRR7279630")

present_reads <- short_reads %in% srrs$Run
present_reads

write.table(short_reads, file = file.path(getwd(), "nextflow", "fetchngs", "input", "short_reads.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Long reads:
long_reads <- c("SRR10343551", "SRR11073828", "SRR10343557", "SRR11073834", "SRR11113509")

present_reads <- long_reads %in% srrs$Run
present_reads

write.table(long_reads, file = file.path(getwd(), "nextflow", "fetchngs", "input", "long_reads.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Combined:
combined_reads <- c(short_reads, long_reads)

write.table(combined_reads, file = file.path(getwd(), "nextflow", "fetchngs", "input", "combined_reads.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# # Load burst phenotype data
# diversity_panel <- read.csv('data/pollen_diversity_panel_phenotypes.csv')
# 
# # Save columns of interest in a new dataframe
# dp_burst <- diversity_panel[, c('accession', 
#                                 'adjusted_integral_increase', 
#                                 'burst_integral_34C_adjusted_mean')] %>%
#   rename(name_CW = accession)
# 
# # Join var_sheet with burst phenotype data
# var_sheet_phenos <- var_sheet %>%
#   left_join(dp_burst, by = "name_CW") %>%
#   drop_na()
# 
# # Load and process pistil-anther ratio data
# flower_phenotypes <- read.csv("data/flower_phenotypes.csv") %>%
#   group_by(accession_id) %>%
#   summarize(mean_anther_pistil_ratio = mean(anther_over_pistil)) %>%
#   rename(name_CW = accession_id)
# 
# var_sheet_phenos <- var_sheet_phenos %>%
#   left_join(flower_phenotypes, by = "name_CW")
# 
# # Load and process locule number data
# locule_num <- read.table(file = file.path(getwd(), "data", "varitome_locule_number.fam"),
#                          sep = '\t', header = FALSE) %>%
#   select(c(1, 6)) %>%
#   rename(name_original_razifard = V1, locule = V6)
# 
# locule_num <- locule_num %>%
#   left_join(var_sheet, by = c("name_original_razifard" = "Sample.Name")) %>%
#   select(accession_id = name_CW, locule) %>%
#   drop_na() %>%
#   rename(name_CW = accession_id)
# 
# var_sheet_phenos <- var_sheet_phenos %>%
#   left_join(locule_num, by = "name_CW")
# 
# # Load and join tube length data
# tube_length <- diversity_panel[, c('accession', 
#                                    'tube_length_26C_adjusted_mean', 
#                                    'tube_length_34C_adjusted_mean')] %>%
#   rename(name_CW = accession)
# 
# var_sheet_phenos <- var_sheet_phenos %>%
#   left_join(tube_length, by = "name_CW")
# 
# #Load remainder of the phenotypes:
# pheno_fill <- diversity_panel[, c('accession',
#                                   'burst_integral_26C_adjusted_mean',
#                                   'burst_2h_34C_adjusted_mean',
#                                   'burst_2h_26C_adjusted_mean')] %>%
#   rename(name_CW = accession)
# 
# var_sheet_phenos <- var_sheet_phenos %>%
#   left_join(pheno_fill, by = "name_CW")
# 
# # Drop all NAs
# var_sheet_phenos <- var_sheet_phenos %>%
#   drop_na()


#var_sheet_phenos now has everything necessary to build paths for fastqs, build samples.tsv, build pheno.tsv for multiple different phenos


# nextflow inputs (paths for hpc)------------------------------------------------------

# loading necessary libraries
library(dplyr)
library(stringr)
library(tidyr)

# Set working directory
setwd('/Users/cperkins/Desktop/palanivelu_lab/kmers_gwas/')

# dp_sheet has all necessary info for pulling all the SRR IDs that match exactly with the SRRs we have for our accessions.
# Since there are multiple SRRs for each accession, it is critical that the SRR we chose in dp_sheet matches the rest of the inputs

# Load the pre-created dp_sheet
dp_sheet <- read.csv(file.path(getwd(), "data", "dp_sheet.csv"), row.names = NULL)

# Path to the fastq files
base_path <- "/xdisk/rpalaniv/cjperkins1/kmers-gwas/genomes/dp/test/fastq"

# Step 1: Generate fq1 and fq2 conditionally
samples <- dp_sheet %>%
  select(Run, Sample.Name, name_CW, Experiment, read_type) %>%
  mutate(
    fq1 = if_else(
      read_type == "paired",
      file.path(base_path, paste0(Experiment, '_', Run, '_1.fastq.gz')),
      file.path(base_path, paste0(Experiment, '_', Run, '.fastq.gz'))
    ),
    fq2 = if_else(
      read_type == "paired",
      file.path(base_path, paste0(Experiment, '_', Run, '_2.fastq.gz')),
      NA_character_
    )
  )

# Step 2: Pivot to long format and assign paired_end based on file suffix
samples_long <- samples %>%
  select(name_CW, fq1, fq2, Run, read_type) %>%
  pivot_longer(cols = c(fq1, fq2), names_to = "fq_column", values_to = "fq") %>%
  filter(!is.na(fq)) %>%
  mutate(
    paired_end = case_when(
      read_type == "single" ~ 0,
      grepl("_1\\.fastq\\.gz$", fq) ~ 1,
      grepl("_2\\.fastq\\.gz$", fq) ~ 2,
      TRUE ~ NA_integer_
    )
  ) %>%
  rename(
    accession_id = name_CW,
    srr_id = Run
  ) %>%
  select(accession_id, fq, paired_end, srr_id)

# Save as .tsv
write.table(samples_long, file = file.path(getwd(), "nextflow", "kmc", "input", "dp_samplesheet.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Pull out the short read test

# Define your SRR list
short_srr_ids <- c("SRR7279519", "SRR7279511", "SRR7279567", "SRR7279572", "SRR7279630")

# Filter samples_long
short_test_samples <- samples_long %>%
  filter(srr_id %in% short_srr_ids)

# Save as .tsv
write.table(short_test_samples, file = file.path(getwd(), "nextflow", "kmc", "input", "short_read_test.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Pull out the long read test

# Define your SRR list
long_srr_ids <- c("SRR10343551", "SRR11073828", "SRR10343557", "SRR11073834", "SRR11113509")

# Filter samples_long
long_test_samples <- samples_long %>%
  filter(srr_id %in% long_srr_ids)

# Save as .tsv
write.table(long_test_samples, file = file.path(getwd(), "nextflow", "kmc", "input", "long_read_test.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
