# HPC, RStudio & GitHub Setup

This document explains how to:  
1. Connect to RStudio Server on CHPC (HPC cluster)  
2. Submit jobs via Slurm (interactive and batch)  
3. Link your GitHub repository with RStudio for version control  

---

## 1. RStudio Server on CHPC

1. Navigate to the CHPC OnDemand portal:  
   https://ondemand.chpc.utah.edu/pun/sys/dashboard/batch_connect/sys/rstudio_server_app/session_contexts/new

2. Configure a new RStudio session:  
   - **Cluster:** notchpeak  
   - **R version:** R 4.4.0 (Geospatial packages)  
   - **Number of cores:** 4 (up to 32 available)  
   - **Wall time (hours):** 100 (max 336 per allocation; choose at least 24)  
   - **Account:** saarman-np  
   - **Partition:** saarman-shared-np (allows multiple simultaneous jobs)  
   - **Memory per job:** 100G (cluster limit: 1000G total; avoid exceeding half)

3. Start the session. Once running, click “Connect to RStudio” to open.

4. To re-enter a running session, visit:  
   https://ondemand.chpc.utah.edu/pun/sys/dashboard/batch_connect/sessions

---

## 2. Slurm Interactive Session

To allocate resources for an interactive shell job, run the following in a terminal (each backslash represents a line continuation):

    salloc --time=1:00:00 \
           --ntasks=1 \
           --mem=100G \
           --account=saarman-np \
           --partition=saarman-shared-np

Once the node is allocated, launch an interactive R session by typing:

    R

---

## 3. Submitting Batch Jobs with Slurm

Create an SBATCH script (for example, save as scripts/run_analysis.sbatch) with contents similar to the following (each line beginning with #SBATCH is an SBATCH directive):

    #!/bin/bash
    #SBATCH --job-name=uganda_analysis
    #SBATCH --account=saarman-np
    #SBATCH --partition=saarman-shared-np
    #SBATCH --time=24:00:00
    #SBATCH --ntasks=1
    #SBATCH --mem=100G
    #SBATCH --output=logs/uganda_%j.out

    module load R/4.0.3

    Rscript scripts/01_prepare_data.R

Then submit the job with:

    sbatch scripts/run_analysis.sbatch

---

## 4. GitHub ⇄ RStudio Integration

### 4.1 Create a Personal Access Token (PAT) on GitHub

1. In your browser, go to:  
   https://github.com/settings/tokens

2. Click “Generate new token” (give it a descriptive name like “CHPC RStudio”) and select at least these scopes:  
   - repo (for access to public/private repositories)  
   - workflow (if you plan to use GitHub Actions)  
   - any additional scopes you require

3. Copy the generated token and store it securely (you will need it in the next steps).

### 4.2 Enable Git in RStudio Server

1. In RStudio (within your CHPC session), open a new R console.  
2. Install the usethis package if it’s not already installed by running:  

    install.packages("usethis")

3. Run the following command to cache your PAT:  

    usethis::browse_github_token()

   - A browser window will open; paste your PAT when prompted.  
   - After saving, RStudio will be able to authenticate Git operations automatically.

### 4.3 Clone the Repository in RStudio

1. In RStudio, go to File → New Project → Version Control → Git.  
2. For the repository URL, enter:  
   https://github.com/your-username/uganda-tsetse-LG.git  
3. Choose a local directory on the CHPC file system where you want to store the project.

RStudio will clone the repository and open it as a project.

### 4.4 Commit, Pull & Push

- In the Git pane (usually in the top right of RStudio), you can:  
  - Stage changed files by checking the box next to each.  
  - Commit your changes with a concise commit message.  
  - Pull to fetch and merge upstream changes from GitHub.  
  - Push to send your local commits to the remote GitHub repository.

---

_This document assumes you have an active CHPC account with the “saarman-np” project allocation and valid GitHub permissions to read/write the uganda-tsetse-LG repository._

