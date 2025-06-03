# Landscape Genetics of Tsetse in Uganda (uganda-tsetse-LG)

**Brief overview:**  
This repository supports population-genetics and landscape-genetics analyses of *Glossina fuscipes fuscipes* (tsetse fly) in Uganda.  

Core objectives include mapping genetic clusters, modeling connectivity via Random Forest, and conducting leave-one-point-out cross-validation. The repository contains all code and small metadata files for our landscape‐genetics analyses of *Glossina fuscipes fuscipes* in Uganda.  Large raw genomic and geospatial inputs, as well as final outputs, are stored externally (see below).

## Folder structure

- **`scripts/`**  
  – R, Python, and Bash scripts in numbered order (e.g., `01_prepare_data.R`, `02_run_PCA.R`).

- **`docs/`**  
  Long-form documentation:
  - `HPC_instructions.md` (RStudio + Slurm setup)
  - `Git_RStudio_setup.md` (GitHub⇄RStudio details)
  - `GoogleDrive_links.md` (external data locations)
  - `TOC_scripts.md` (annotated list of all major analysis scripts)

- **`metadata/`**  
  – Small reference files (e.g., sample metadata, lookup tables).
  - Other information/instructions:
  1. `HPC_instructions.md` – How to start RStudio/Slurm on CHPC.  
  2. `Git_RStudio_setup.md` – How to link GitHub and RStudio (personal access token, etc.).  
  3. `GoogleDrive_links.md` – Where to find TIFs, external data sources.  
  4. `scripts_TOC.md` – Annotated list of each major script and its purpose.

- **`data/`** (ignored by GitHub; see `data/README.md` for external locations)

- **`results/`** (ignored by GitHub; see `results/README.md` for external locations)

- **`figures/`**  (ignored by GitHub)
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
