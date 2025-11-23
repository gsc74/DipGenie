
# Reproduce Results

## 1. Prerequisites

Before starting the process, ensure the following tools and datasets are available:

- **Miniforge**: [Miniforge](https://github.com/conda-forge/miniforge)
- **Required Tools**:
    - DipGenie
    - seqkit
    - bcftools (with plugins)
    - gfa2gbwt

Ensure all dependencies are installed or available via conda.

---

## 2. Prerequisites

### Download the datasets and install DipGenie
```bash
git clone https://github.com/gsc74/DipGenie
cd DipGenie

# Set installation root
export DipGenie_HOME="$(pwd)"

# Install dependencies (Miniforge is required)
"$DipGenie_HOME"/Installdeps

# Environment setup
export PATH="$DipGenie_HOME/extra/bin:$PATH"
export LD_LIBRARY_PATH="$DipGenie_HOME/extra/lib:$LD_LIBRARY_PATH"

# Build
make

# Test diploid inference
./DipGenie -t32 -p2 -g test/MHC_4.gfa.gz -r test/HG002.mhc.2x.fq.gz -o HG002

# Download the datasets from Zenodo
wget zenodo_link -O data.zip && unzip data.zip && Unpack_commands
```

---

## 3. Running DipGenie

For ILP and IQP, execute the following scripts:

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

## 5. Evalaute results with whatshap for SER

```bash

```

---

## 7. Postprocessing

After the runs are completed, extract relevant metrics:

```bash
# ILP 
python3 postprocessing_3_MILP.py 

# ILP (no relaxation)
python3 postprocessing_3.py 

# IQP 
python3 postprocessing_2_MIQP.py 

# IQP (no relaxation)
python3 postprocessing_2.py

# Progressive imputation
python3 postprocessing_9.py # 13 haps
python3 postprocessing_10.py # 25 haps
python3 postprocessing_11.py # 49 haps
python3 postprocessing_12.py # 7 haps
python3 postprocessing_13.py # 3 haps

# VG
python3 postprocessing_VG.py

# PanGenie
python3 postprocessing_PG.py
```
---