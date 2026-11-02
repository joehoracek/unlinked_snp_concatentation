from Bio import SeqIO
import csv
import glob
import os
import random


def get_all_snps(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        raise ValueError("No sequences found in the FASTA file.")

    seq_length = len(records[0].seq)
    for rec in records:
        if len(rec.seq) != seq_length:
            raise ValueError("Not all sequences are the same length. Is your FASTA aligned?")

    snp_columns = []
    for i in range(seq_length):
        column = [str(rec.seq[i]).upper() for rec in records]
        unique_bases = set(column)
        if "-" in unique_bases:
            unique_bases.discard("-")
        if len(unique_bases) > 1:
            snp_columns.append(column)

    # Convert to dict: specimen -> list of SNPs
    snp_dict = {rec.id: [] for rec in records}
    for col in snp_columns:
        for j, rec in enumerate(records):
            snp_dict[rec.id].append(col[j])

    return snp_dict, len(snp_columns)




def calculate_maf(column):
    counts = {}
    for base in column:
        base = base.upper()
        if base in ("-", "N"):
            continue
        counts[base] = counts.get(base, 0) + 1

    if len(counts) <= 1:
        return 0.0  # monomorphic

    total = sum(counts.values())
    freqs = [c / total for c in counts.values()]
    return min(freqs)  # minor allele frequency

def get_best_snp(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        return None, None

    seq_length = len(records[0].seq)
    for rec in records:
        if len(rec.seq) != seq_length:
            raise ValueError(f"Not all sequences are the same length in {fasta_file}.")

    # Get all polymorphic SNP columns
    snp_dict_all, num_snps = get_all_snps(fasta_file)
    if num_snps == 0:
        # No SNPs at all
        print(f"No SNPs found in {fasta_file}. Filling with '-' for all vouchers.")
        return {rec.id: "-" for rec in records}, 0

    # Identify informative SNPs (not all N/gap)
    informative_positions = []
    for i in range(num_snps):
        column = [snp_dict_all[rec.id][i] for rec in records]
        if set(column) - set(["N", "-"]):
            informative_positions.append(i)

    if informative_positions:
        # Calculate MAF for each informative SNP and pick the highest
        best_maf = -1
        best_pos = None
        for i in informative_positions:
            column = [snp_dict_all[rec.id][i] for rec in records]
            maf = calculate_maf(column)
            if maf > best_maf:
                best_maf = maf
                best_pos = i
        best_column = {rec.id: snp_dict_all[rec.id][best_pos] for rec in records}
    else:
        # All SNPs are only N/gaps -> pick one SNP column at random
        chosen_pos = random.randint(0, num_snps - 1)
        best_column = {rec.id: snp_dict_all[rec.id][chosen_pos] for rec in records}
        best_maf = 0.0
        print(f"All SNPs in {fasta_file} are only N/gaps. Random SNP selected at position {chosen_pos}.")

    return best_column, best_maf


if __name__ == "__main__":
    fasta_folder = "bol-loci"
    voucher_file = "bol-voucher.csv"  # CSV with column "voucher"
    output_csv = "all_loci_snps.csv"
    output_fasta = "all_loci_snps.fasta"

    # Load voucher list
    vouchers = []
    with open(voucher_file, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            vouchers.append(row["voucher"])

    fasta_files = sorted(glob.glob(os.path.join(fasta_folder, "*.fasta")))

    if not fasta_files:
        print("No FASTA files found in the folder.")
    else:
        all_snps = {v: [] for v in vouchers}
        loci_names = []

        for fasta_file in fasta_files:
            locus_name = os.path.basename(fasta_file).replace(".fasta", "")
            loci_names.append(locus_name)

            snp_dict, maf = get_best_snp(fasta_file)
            if snp_dict is None:
                print(f"No SNPs found in {fasta_file}. Filling with '-' for all vouchers.")
                for v in vouchers:
                    all_snps[v].append("-")
                continue

            # Fill in SNPs, using '-' for missing specimens
            for v in vouchers:
                base = snp_dict.get(v, "-")
                all_snps[v].append(base)

            print(f"Processed {fasta_file} (MAF={maf})")

        # Write combined CSV
        with open(output_csv, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Specimen"] + loci_names)
            for v in vouchers:
                writer.writerow([v] + all_snps[v])

        print(f"All SNPs saved to {output_csv}")

        # Write combined FASTA
        with open(output_fasta, "w") as fasta_out:
            for v in vouchers:
                seq = "".join(all_snps[v])
                fasta_out.write(f">{v}\n{seq}\n")

        print(f"FASTA alignment saved to {output_fasta}")
