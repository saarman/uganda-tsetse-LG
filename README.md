# Landscape Genetics of Tsetse in Uganda (uganda-tsetse-LG)

**Brief overview:**  
This repository supports population-genetics and landscape-genetics analyses of *Glossina fuscipes fuscipes* (tsetse fly) in Uganda.  

Core objectives include mapping genetic clusters, modeling connectivity via Random Forest, and conducting leave-one-point-out cross-validation. The repository contains all code and small metadata files for our landscape‐genetics analyses of *Glossina fuscipes fuscipes* in Uganda.  Large raw genomic and geospatial inputs, as well as final outputs, are stored externally (see below).

## Folder structure

- **`scripts/`**  (see `scripts/README.md` for table of contents)   
  – R, Python, and Bash scripts in numbered order (e.g., `01_prepare_data.R`, `02_run_PCA.R`).  
  
- **`docs/`**  (Long-form documentation; see `docs/README.md` for details)  
      1. Connecting to RStudio on CHPC  
      2. Submitting Slurm (interactive & batch) jobs  
      3. Linking GitHub with RStudio (creating a PAT, cloning, committing, pushing)  

- **`metadata/`** (e.g., sample metadata, lookup tables)      
  – Small reference files (e.g., sample metadata, lookup tables).      
  – Large raw genomic and geospatial inputs, as well as final outputs, are stored externally (see below).
          
- **`data/`** (ignored by GitHub; see `data/README.md` for external locations)

- **`results/`** (ignored by GitHub; see `results/README.md` for external locations)

- **`figures/`**  (ignored by GitHub; see `figures/README.md` for list of figures)
  High-resolution output figures from published/popgen/ML analyses.


## Getting started

1. Clone the repo:
    ```bash
    git clone https://github.com/your-username/uganda-tsetse-LG.git
    ```

2. Install required R packages:
    ```bash
    Rscript scripts/00_install_packages.R
    ```
3. Follow the pipeline:
    - Run `01_prepare_data.R` to load and format data.
    - Then run subsequent scripts in numeric order (e.g., `02_run_PCA.R`, `03_build_RF_model.R`, etc.).

4. For HPC setup (RStudio + Slurm), see:
   `docs/HPC_instructions.md`

5. For GitHub ⇄ RStudio integration, see:
   `docs/Git_RStudio_setup.md`


## Citation & License

If you use these scripts or refer to our data, please cite:  
> Saarman NP et al. (2025) “Landscape genetics of *Glossina fuscipes fuscipes* in Uganda,” In prepartion for publication in *Insects: Spatial Population Genetics in Insects*, [URL:https://www.mdpi.com/journal/insects/special_issues/C9FL984799](https://www.mdpi.com/journal/insects/special_issues/C9FL984799).

Licensed under MIT. See `LICENSE` for full terms.
