<div align="center">
  <!-- <img src="test/logo/logo_phi.png" alt="PHI Logo" width="200"> -->
</div>

## <div align="center"><span style="color:red;"><b>DipGenie</b></span> (<span style="color:red;"><b>D</b></span>iploid <span style="color:red;"><b>G</b></span>enome <span style="color:red;"><b>i</b></span>nference)</div>


## <a name="started"></a>Getting Started

### Prerequisites

Before running **DipGenie**, make sure **Miniforge** is installed: [Miniforge Installation Guide](https://github.com/conda-forge/miniforge). This is required for installing dependencies such as **VG**, **BayesOpt**, and **samtools**.


## <a name="get_phi"></a>Get DipGenie

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

# Test runs
# For haploid
./DipGenie -t32 -p1 -g test/MHC_4.gfa.gz -r test/CHM13_reads.fq.gz -o CHM13
# For diploid
./DipGenie -t32 -p2 -g test/MHC_4.gfa.gz -r test/HG002.mhc.2x.fq.gz -o HG002
```
