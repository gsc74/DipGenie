
# Reproduce Results

## 1. Prerequisites

Before starting the process, ensure the following tools and datasets are available:

- **Required Tools**:
    - DipGenie
    - seqkit
    - bcftools (with plugins)
    - gfa2gbwt
    - agc
    - whatshap
    - truvari

Ensure all dependencies are installed or available via conda.

---

## 2. Prerequisites

### Download the datasets and install DipGenie
```bash
git clone https://github.com/gsc74/DipGenie
cd DipGenie
make

# Test diploid inference
./DipGenie -t32 -p2 -g test/MHC_4.gfa.gz -r test/HG002.mhc.2x.fq.gz -o HG002

# Download the datasets
wget https://zenodo.org/api/records/17685087/files-archive -O data.zip
unzip data.zip
unzip Graphs.zip -d Graphs
unzip Truth.zip
unzip hprc_haps.zip
unzip Reads.zip -d Reads
```

---

## 3. Running DipGenie

```bash
# run for different set of samples
sh run_DipGenie_batch.sh samples_batch_1.txt
sh run_DipGenie_batch.sh samples_batch_2.txt
sh run_DipGenie_batch.sh samples_batch_3.txt
sh run_DipGenie_batch.sh samples_batch_4.txt
sh run_DipGenie_batch.sh samples_batch_5.txt
```

---

## 4. Running VG

```bash
# run for different set of samples
sh run_VG_batch.sh samples_batch_1.txt
sh run_VG_batch.sh samples_batch_2.txt
sh run_VG_batch.sh samples_batch_3.txt
sh run_VG_batch.sh samples_batch_4.txt
sh run_VG_batch.sh samples_batch_5.txt
```

---

## 5. Generate truth and test VCFs
```bash
sh gen_VCFs_VG.sh
sh gen_VCFs_DipGenie.sh

# Optionally generate the ground truth VCFs (Truth VCFs are already provided in datasets)
sh gen_VCFs_truth.sh
```

---

## 6. Evalaute results

```bash
sh get_F1_DipGenie.sh
sh get_SER_DipGenie.sh
sh get_F1_VG.sh
sh get_SER_VG.sh
```

## 6. Print results

```bash
sh print_F1_DipGenie.sh
sh print_len_DipGenie.sh
sh print_SER_DipGenie.sh
sh print_SVs_DipGenie.sh
sh print_F1_VG.sh
sh print_len_VG.sh
sh print_SER_VG.sh
sh print_SVs_VG.sh
```

---
