#!/usr/bin/env python3
"""
maketsv.py
Extracts CHROM, POS, REF, A1_ALLELE (CC2192B4), and A2_ALLELE (RCC3703) from a VCF file.

Usage:
    python maketsv.py input.vcf output.tsv

The script auto-detects which column index corresponds to the CC2192B4 (A1) and RCC3703 (A2) samples
by reading the #CHROM header line.
"""

import sys

def get_allele(ref, alt_list, genotype_field):
    """
    Resolve a sample's genotype to the actual allele sequence.
    Handles GT fields like '0/0', '0/1', '1/1', '0|1', './.', etc.
    Returns the allele string (e.g. 'A', 'T', 'ATCG') or '.' if missing.
    """
    gt = genotype_field.split(":")[0]       # GT is always first sub-field
    gt = gt.replace("|", "/")               # treat phased same as unphased
    alleles = [ref] + alt_list              # index 0 = REF, 1+ = ALT alleles

    calls = gt.split("/")
    resolved = []
    for c in calls:
        if c == ".":
            return "."                      # missing genotype
        idx = int(c)
        if idx < len(alleles):
            resolved.append(alleles[idx])
        else:
            return "."                      # malformed index

    # Return unique allele (homozygous → one value; het → both joined by '/')
    unique = list(dict.fromkeys(resolved))  # preserve order, deduplicate
    return "/".join(unique)


def vcf_to_tsv(vcf_path, tsv_path,
               strain1_name="CC2192B4", strain2_name="RCC3703",
               chrom_filter="CM083880.1"):

    strain1_idx = None
    strain2_idx = None

    with open(vcf_path, "r") as vcf, open(tsv_path, "w") as tsv:
        # Write header
        tsv.write("\t".join(["CHROM", "POS", "REF", "A1_ALLELE", "A2_ALLELE"]) + "\n")

        for line in vcf:
            # Skip meta-information lines
            if line.startswith("##"):
                continue

            fields = line.rstrip("\n").split("\t")

            # Parse the column header to find sample indices
            if line.startswith("#CHROM"):
                sample_names = fields[9:]           # samples start at column 9
                try:
                    strain1_idx = sample_names.index(strain1_name) + 9
                except ValueError:
                    sys.exit(f"ERROR: Sample '{strain1_name}' not found in VCF header.\n"
                             f"Available samples: {sample_names}")
                try:
                    strain2_idx = sample_names.index(strain2_name) + 9
                except ValueError:
                    sys.exit(f"ERROR: Sample '{strain2_name}' not found in VCF header.\n"
                             f"Available samples: {sample_names}")
                continue

            # Data lines
            if strain1_idx is None or strain2_idx is None:
                sys.exit("ERROR: No #CHROM header line found before data lines.")

            chrom   = fields[0]

            # Skip records not on the target chromosome
            # if chrom != chrom_filter:
            #     continue

            pos     = fields[1]
            ref     = fields[3]
            alt_str = fields[4]
            alt_list = alt_list = [] if alt_str == "." else alt_str.split(",")

            a1 = get_allele(ref, alt_list, fields[strain1_idx])
            a2 = get_allele(ref, alt_list, fields[strain2_idx])

            tsv.write("\t".join([chrom, pos, ref, a1, a2]) + "\n")

    print(f"Done. Output written to: {tsv_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python vcf_to_tsv.py <input.vcf> <output.tsv>")
        sys.exit(1)

    vcf_to_tsv(sys.argv[1], sys.argv[2], chrom_filter="CM083880.1")