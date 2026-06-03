library(poppr)
library(vcfR)
library(adegenet)

metafile <- '../../../../hybrid_diversity/5_snp_analyses/6_hybrid_analyses/metadata.txt'
vcffile <- '../../../../hybrid_diversity/5_snp_analyses/2_atype_network/atype_variants_qualfiltered_maffiltered_renamed_justhybrids.recode.vcf'
permutations <- 999

cat("Loading VCF file...\n")
# Read the VCF - adjust path as needed
vcf <- read.vcfR(vcffile)
 
# Convert to genind object (required for AMOVA)
cat("Converting to genind object...\n")
gen <- vcfR2genind(vcf)
 
cat("Sample names in VCF:\n")
print(indNames(gen))

cat("\nLoading metadata...\n")
metadata <- read.table(metafile, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
print(head(metadata))

# Check that all VCF samples have metadata
missing_meta <- setdiff(indNames(gen), metadata$strain_ID)
if (length(missing_meta) > 0) {
  stop("Missing metadata for strains: ", paste(missing_meta, collapse = ", "))
}

# Reorder metadata to match VCF sample order
metadata <- metadata[match(indNames(gen), metadata$strain_ID), ]

print(head(metadata))

# Attach metadata as strata to the genind object
strata_df <- data.frame(
  strain_ID = metadata$strain_ID,
  cytotype = factor(metadata$cytotype),
  region = factor(metadata$region)
)
strata(gen) <- strata_df
 
cat("\nStrata summary:\n")
print(table(strata_df$cytotype, strata_df$region))

cat("\nPermutations:\n")
print(permutations)

# Set population structure and run AMOVA without within-individual variance
setPop(gen) <- ~region/cytotype
amova_result1 <- poppr.amova(gen, ~region/cytotype, method = "pegas", within = FALSE, nperm = permutations)
cat("\nAMOVA ~region/cytotype\n")
print(amova_result1)

# Test the other hierarchy
setPop(gen) <- ~cytotype/region  
amova_result2 <- poppr.amova(gen, ~cytotype/region, method = "pegas", within = FALSE, nperm = permutations)
cat("\nAMOVA ~cytotype/region\n")
print(amova_result2)

# Set population structure and run AMOVA without within-individual variance
amova_result3 <- poppr.amova(gen, ~region, method = "pegas", within = FALSE, nperm = permutations)
cat("\nAMOVA just region\n")
print(amova_result3)

# Test the other hierarchy
amova_result4 <- poppr.amova(gen, ~cytotype, method = "pegas", within = FALSE, nperm = permutations)
cat("\nAMOVA just cytotype\n")
print(amova_result4)



