#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <reference.fa> <hap1.fa> <hap2.fa> <output.vcf>"
    exit 1
fi

# Assign command-line arguments to variables
REFERENCE=$1
HAP1=$2
HAP2=$3
OUTPUT_VCF=$4

# Check if the input files exist
if [[ ! -f "$REFERENCE" || ! -f "$HAP1" || ! -f "$HAP2" ]]; then
    echo "Error: Make sure reference.fa, hap1.fa, and hap2.fa are provided and exist."
    exit 1
fi

# Step 1: Index the reference sequence
echo "Indexing reference sequence..."
samtools faidx $REFERENCE

# Step 2: Align haplotypes to the reference using minimap2
echo "Aligning haplotypes to reference with minimap2..."
minimap2 -ax asm5 $REFERENCE $HAP1 | samtools view -bS - > hap1.bam
minimap2 -ax asm5 $REFERENCE $HAP2 | samtools view -bS - > hap2.bam

# Step 3: Sort and index the BAM files
echo "Sorting and indexing BAM files..."
samtools sort hap1.bam -o hap1_sorted.bam
samtools sort hap2.bam -o hap2_sorted.bam
samtools index hap1_sorted.bam
samtools index hap2_sorted.bam

# Step 4: Merge the BAM files to simulate a diploid organism
echo "Merging BAM files to simulate diploid sample..."
samtools merge merged.bam hap1_sorted.bam hap2_sorted.bam
samtools index merged.bam

# Step 5: Generate pileup and call variants with bcftools
echo "Calling variants with bcftools..."
bcftools mpileup -f $REFERENCE merged.bam -Ou -o merged.bcf
bcftools call -vmO v -o raw_output.vcf merged.bcf

# Step 6: Filter the VCF file for quality (optional step)
echo "Filtering VCF file..."
# bcftools filter -e 'QUAL<20 || DP<10' raw_output.vcf -o $OUTPUT_VCF
cp raw_output.vcf $OUTPUT_VCF

# Clean up intermediate files
echo "Cleaning up intermediate files..."
rm hap1.bam hap2.bam hap1_sorted.bam hap2_sorted.bam merged.bam merged.bcf raw_output.vcf

# Final message
echo "Genotyped VCF file generated: $OUTPUT_VCF"

