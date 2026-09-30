## load packages and data
if(!require(pacman)) install.packages(pacman)
pacman::p_load(ape, linelist, epitrix, fitdistrplus,
               outbreaker2, tidyverse, rio,
               magrittr, remotes, epicontacts,
               visNetwork, lubridate, pairsnp,
               ggplot2, chron,glue, here)


library("remotes")

remotes::install_github("finlaycampbell/outbreaker2@timeline", force = TRUE)



setwd("/Users/akanteh/Documents/backup11:08:2022/phylodynamic_newAnalysis/outbreak_analysis/")

source("src/functions_wards.R")
source("src/functions.R")
source("src/onset_distribution.R")
source("src/functions_analyse.R")


scale <- 6
reporting.probability <- 0.5
number.import.chains <- 5
output.name <- "w1_i0.7_pi0.5"

## make unique file output name.
output.name <- paste0(output.name, "_",
                      format(Sys.time(), "%m%d%H%M%S"), "_",
                      sample(1:99999, 1))


# load data, clean dates and filter-out observations
mydata <- read.csv("data/metaData.csv", stringsAsFactors = FALSE )


mydata$DateOfCollection <- as.Date(mydata$DateOfCollection, format="%d/%m/%Y")
mydata$OnsetDate <- as.Date(mydata$OnsetDate, format="%d/%m/%Y")

# Convert to Date format if not already
mydata <- mydata %>%
  mutate(DateOfCollection = ymd(DateOfCollection),
         wave = case_when(
           year(DateOfCollection) == 2020 ~ "wave1",
           year(DateOfCollection) == 2021 ~ "wave2",
           TRUE ~ NA_character_
         ))

#### create wave 1, wave2, sequence only or JulyOnly dataset (explore)
mydata_wave1 <- mydata[mydata$wave %in% c("wave1"),]

mydata_wave2 <- mydata[mydata$wave %in% c("wave2"),]

sequnce_only <- mydata[1:20,]

julyOnly <- mydata %>%
  filter(year(DateOfCollection) == 2020, month(DateOfCollection) == 7)



#################################
#process onset dates
#################################
date_onset<- mydata_wave1$OnsetDate
names(date_onset) <- rep("local", length(date_onset))

#################################
#load dna data in DNAbin format,
#################################
#one sequence per individual
dna <- read.FASTA("analysis/oubreakerSamplesfirstwave.fasta")

# Assign meaningful names (e.g. sample IDs)
names(dna) <- mydata_wave1$lab_no[1:16] ### only first 16 samples have dna sequnces

# Save and re-import as sparse matrix, calculate pairwise SNP distances
write.FASTA(dna, paste0("data/", output.name, "error-removed_masked.ord.aln"))
dna.SNP <- import_fasta_sparse(paste0("data/", output.name, "error-removed_masked.ord.aln")) %>% snp_dist(.)

# Assign row and column names correctly
mylabels <- names(dna)
rownames(dna.SNP) <- colnames(dna.SNP) <- mylabels
rm(mylabels)
unlink(paste0("data/", output.name,  "error-removed_masked.ord.aln"))



#create meta object
meta <- mydata_wave1[,c("lab_no","OnsetDate")]
names(meta) <- c("id", "date_onset")
row.names(meta) <- NULL

# Incubation (infection -> onset) and Generation time (infection -> infection)
incub <- distcrete::distcrete("gamma", shape = 5.89, scale = 1.068555, interval = 1)
gen   <- distcrete::distcrete("gamma", shape = 4.83,  scale = 1.051851, interval = 1)


# Use directly (no extra scaling) – 1..30 days carries most mass
w_dens <- gen$d(1:30)
f_dens <- cbind(local = incub$d(1:30))


data <- outbreaker_data(
  dates = date_onset,
  dna = dna,
  ids = meta$id,
  w_dens = w_dens,
  f_dens = f_dens,
  genetic_model = "MRCA"
)


cfg <- outbreaker2::create_config(
  n_iter        = 1000000,
  sample_every  = 500
  
)

run_with_Only_symptomsOnset <- outbreaker(data = data, config = cfg)

run_with_with_symptomsOnset_DNA <- outbreaker(data = data, config = cfg)

smry_basic_symtomOnly <- summary(Pre_run_with_Only_symptomsOnset)
hist(smry_basic_symtomOnly$tree$support, col = "grey", border = "white", 
     main = "Consensus ancestry: support", xlim = c(0,1))

smry_basic_withDNA <- summary(Pre_run_with_with_DNA)
hist(smry_basic_withDNA$tree$support, col = "grey", border = "white", 
     main = "Consensus ancestry: support", xlim = c(0,1))




# All cases (same order as your outbreaker dates/ids)
ids_all   <- mydata_wave1$lab_no
dates_all <- mydata_wave1$OnsetDate

# The IDs that actually have sequences (names(dna) already set)
seq_ids <- names(dna)

# Ensure dna.SNP is square and named on seq_ids only
dna.SNP <- dna.SNP[seq_ids, seq_ids, drop = FALSE]

# Build full N x N (NA-filled) and drop the sequenced block in the right rows/cols
N <- length(ids_all)
D_full <- matrix(NA_real_, N, N, dimnames = list(ids_all, ids_all))
D_full[match(seq_ids, ids_all), match(seq_ids, ids_all)] <- as.matrix(dna.SNP)




data <- outbreaker2::outbreaker_data(
  dates  = date_onset,   # onset dates
  ids    = meta$id,
  dna    = dna,         # only the sequenced subset (length = number of seqs)
  w_dens = w_dens,
  f_dens = f_dens,
  D      = D_full       # N x N with NA for unsequenced pairs
  
)


cfg <- outbreaker2::create_config(
  n_iter        = 1000000,
  sample_every  = 500
  
)

pre_run_onset_DNA_matrix <- outbreaker(data = data, config = cfg)

smry_basic_DNA_matrix <- summary(pre_run_onset_DNA_matrix)
hist(smry_basic_DNA_matrix$tree$support, col = "grey", border = "white", 
     main = "Consensus ancestry: support", xlim = c(0,1))




# ----------------------------------------------
# 4) STAFF LOCATION CONTACT MATRIX (robust)
# ----------------------------------------------

loc.meta <- read.csv("data/metaData.csv", stringsAsFactors = FALSE)

# 1) Standardize location labels to match weight table
std <- function(x) tolower(trimws(x))
mydata_wave1$locations <- std(mydata_wave1$locations)


# 2) Make sure we are in the same order as your model IDs
#    (use data$ids if you’ve already built the outbreaker_data)
ord <- match(mydata_wave1$lab_no, mydata_wave1$lab_no)     # mydata order is the target
if (any(is.na(ord))) stop("Some mydata$lab_no not found in loc.meta$lab_no")
mydata_wave1 <- mydata_wave1[ord, , drop = FALSE]

# 3) Define the location-to-location base weights
weight_table <- matrix(c(
  1,  0, 0, 0,
  0,  1,   0, 0,
  0,  0, 1,   0,
  0, 0,0,1
), 4, byrow = TRUE,
dimnames = list(c("ward","lab","opd","other"),
                c("ward","lab","opd","other")))


# 4) Build the NxN matrix by indexing (no loops)
locs <- mydata_wave1$locations
if (!all(locs %in% rownames(weight_table))) stop("Location label not in weight_table")

contact_matrix <- weight_table[locs, locs, drop = FALSE]

# 5) Optional: downweight pairs far apart in time (e.g., > 14 days)
mydata_wave1$OnsetDate <- as.Date(mydata_wave1$OnsetDate, format="%d/%m/%Y")
if (any(is.na(mydata_wave1$OnsetDate))) warning("Some OnsetDate are NA; time weights will be NA there.")

# Exponential time-decay with scale tau (days)
tau <- 7  # adjust: 7 means ~e^-1 (0.37x) every 7 days apart
dt  <- abs(outer(mydata_wave1$OnsetDate, mydata_wave1$OnsetDate, "-"))


# 6) Clean up and align with outbreaker2 data object
#    Set identical row/col names to your IDs (mydata$lab_no or data$ids)
rownames(contact_matrix) <- mydata_wave1$lab_no
colnames(contact_matrix) <- mydata_wave1$lab_no

# outbreaker2 usually ignores the diagonal; set to 0 to avoid self-bias
diag(contact_matrix) <- 0

# 7) Reorder exactly to the order of your outbreaker data object (if needed)
if (exists("data")) {
  cm_ord <- match(data$ids, rownames(contact_matrix))
  contact_matrix <- contact_matrix[cm_ord, cm_ord, drop = FALSE]
  stopifnot(all(rownames(contact_matrix) == data$ids),
            all(colnames(contact_matrix) == data$ids))
}


# 9) Attach to outbreaker2 data
data$contact <- contact_matrix


cfg <- outbreaker2::create_config(
  n_iter        = 100000,
  sample_every  = 50
  
)

run_withContactMetrixWave1 <- outbreaker(data =data, config =cfg)



### wave 2 dataset

run3withContactMetrixWave2 <- outbreaker(data =data)

### wave 1 dataset
run3withContactMetrixSequenceOnly <- outbreaker(data =data)

### wave JulyOnly dataset
run3withContactMetrixJulyOnly <- outbreaker(data =data, config = cfg )




# Prior on mutation rate μ: still centered ~0.1 but wider than ultra-tight
prior_mu <- function(param) dbeta(param$mu,89.9,809.1, log = TRUE)

priors <- outbreaker2::custom_priors(mu = prior_mu)


cfg <- outbreaker2::create_config(
  n_iter        = 100000,
  sample_every  = 50,
  move_mu       = TRUE,  init_mu = 0.2, sd_mu = 0.05,
  move_pi       = TRUE,  prior_pi = c(5, 20),   # mean ≈ 0.20; relax if needed
  move_eps      = FALSE,
  move_tau      = FALSE,
  move_model    = FALSE,
  sd_t_inf      = 3
  
)


# Quick checks
stopifnot(length(data$dates) == length(data$ids))
stopifnot(all(rownames(data$D) == data$ids), all(colnames(data$D) == data$ids))
if (!is.null(data$contact)) stopifnot(all(dim(data$contact) == c(length(data$ids), length(data$ids))))

# Run
res_withPriorswave1JulyOnly <- outbreaker2::outbreaker(
  data        = data,
  priors      = priors,
  config      = cfg
)

res_withPriorswave1 <- outbreaker2::outbreaker(
  data        = data,
  priors      = priors,
  config      = cfg
)
