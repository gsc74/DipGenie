<div align="center">
  <!-- <img src="test/logo/logo_phi.png" alt="PHI Logo" width="200"> -->
</div>

## <div align="center"><span style="color:red;"><b>DipGenie</b></span> (<span style="color:red;"><b>D</b></span>iploid <span style="color:red;"><b>G</b></span>enome <span style="color:red;"><b>i</b></span>nference)</div>


## <a name="started"></a>Getting Started

### Prerequisites

Before using DipGenie, please ensure that **Miniforge** is installed: [Miniforge Installation Guide](https://github.com/conda-forge/miniforge). This package installer is used for installing a few dependencies such as VG, BayesOpt and samtools.

## <a name="get_phi"></a>Get DipGenie

```bash
git clone https://github.com/gsc74/DipGenie
cd DipGenie
# Install dependencies (Miniforge is required)
./Installdeps
export PATH="$(pwd)/extra/bin:$PATH"
export LD_LIBRARY_PATH="$(pwd)/extra/lib:$LD_LIBRARY_PATH"
make

# test runs 
./DipGenie -t32 -p1 -g test/MHC_4.gfa.gz -r test/CHM13_reads.fq.gz -o CHM13 # For haploid
./DipGenie -t32 -p2 -g test/MHC_4.gfa.gz -r test/HG002.mhc.2x.fq.gz -o HG002 # For diploid
```
