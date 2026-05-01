#!/usr/bin/env python3

import os
import subprocess
import argparse

def run_command(command):
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

def main(temp_gbz, read_file, tmp_dir, threads, output_fasta):
    os.makedirs(tmp_dir, exist_ok=True)

    # Step 1: Create distance index
    run_command(f"vg index -j {tmp_dir}/temp.dist {temp_gbz}")

    # Step 2: Generate r-index
    run_command(f"vg gbwt -p --num-threads {threads} -r {tmp_dir}/temp.ri -Z {temp_gbz}")

    # Step 3: Generate haplotype info
    run_command(f"vg haplotypes -v 2 -t {threads} -H {tmp_dir}/temp.hapl {temp_gbz}")

    # Step 4: KMC k-mer counting
    run_command(f"kmc -k29 -m128 -okff -t{threads} -hp {read_file} {tmp_dir}/sample {tmp_dir}")

    # Step 5: Diploid haplotype sampling
    run_command(f"vg haplotypes --diploid-sampling -v 2 -t {threads} --num-haplotypes 2 "
                f"-i {tmp_dir}/temp.hapl -k {tmp_dir}/sample.kff -g {tmp_dir}/sample.gbz {temp_gbz}")

    # Step 6: Extract paths to FASTA
    run_command(f"vg paths -x {tmp_dir}/sample.gbz -F -S recombination > {output_fasta}")

    # Step 7: Reverse complement
    run_command(f"seqtk seq -r {output_fasta} > {tmp_dir}/temp_rc.fa && mv {tmp_dir}/temp_rc.fa {output_fasta}")

    print(f"Done. Output FASTA: {output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="VG haplotype sampling pipeline.")
    parser.add_argument("-g", "--gbz", type=str, required=True, help="Input GBZ file")
    parser.add_argument("-r", "--reads", type=str, required=True, help="Input read file (interleaved FASTQ)")
    parser.add_argument("-t", "--threads", type=int, default=128, help="Number of threads (default: 128)")
    parser.add_argument("-d", "--tmp-dir", type=str, required=True, help="Temporary directory")
    parser.add_argument("-o", "--output", type=str, default="sample.fa", help="Output FASTA file")

    args = parser.parse_args()
    main(args.gbz, args.reads, args.tmp_dir, args.threads, args.output)
