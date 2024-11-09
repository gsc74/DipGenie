<div align="center">
  <img src="test/logo/logo_phi.png" alt="PHI Logo" width="200">
</div>

## <div align="center"><span style="color:red;"><b>PHI</b></span> (<span style="color:red;"><b>P</b></span>angenome-based <span style="color:red;"><b>H</b></span>aplotype <span style="color:red;"><b>I</b></span>nference)</div>


## <a name="started"></a>Getting Started

### Prerequisites

Before using PHI, please ensure that **Miniforge** is installed: [Miniforge Installation Guide](https://github.com/conda-forge/miniforge). This package installer is used for installing a few dependencies such as VG and samtools. To run PHI, you also need a Gurobi license. You can get a free academic license [here](https://www.gurobi.com/academia/academic-program-and-licenses/). You should download and save `gurobi.lic` file in your home directory.

## <a name="get_phi"></a>Get PHI

```bash
git clone https://github.com/gsc74/PHI2
cd PHI2
# Install dependencies (Miniforge is required)
./Installdeps
export PATH="$(pwd)/extra/bin:$PATH"
export LD_LIBRARY_PATH="$(pwd)/extra/lib:$LD_LIBRARY_PATH"
make

# test run 
./PHI -t32 -g test/MHC_4.gfa.gz -r test/CHM13_reads.fq.gz -o CHM13
```

#### Adding Binary and Library Paths to `.bashrc`
To ensure that the `extra/bin` and `extra/lib` directories are automatically loaded for every terminal session, you can export them to your `~/.bashrc`. This will make sure the required binaries and libraries for `PHI2` are available.

```bash
# Add extra/bin and extra/lib to .bashrc
echo 'export PATH="$(pwd)/extra/bin:$PATH"' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="$(pwd)/extra/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc
source ~/.bashrc
```
