# codespaces_RStudio

Reproducible teaching labs for RNA-seq and R programming built around GitHub Codespaces. This repo now ships a full RStudio Server IDE, a lightweight desktop viewer, and a mamba-managed toolkit so students always land in the same environment regardless of their laptop.

## What you get

- **RStudio Server 2023.x** delivered through the Rocker devcontainer feature so the familiar IDE is available on port 8787 inside every Codespace.
- **XFCE desktop + noVNC bridge** (desktop-lite feature) so you can interact with graphical tools straight from VS Code or the browser, including launching RStudio in a desktop session if preferred.
- **Conda/mamba toolchain** (Miniforge feature) for installing CLI bioinformatics utilities such as FastQC, MultiQC, Salmon, STAR, samtools, seqtk, and bedtools. The default env spec lives at `env/mamba-environment.yml`.
- **Pre-installed R/Bioconductor stack** that covers every package used across the included practicals (`DESeq2`, `airway`, `IsoformSwitchAnalyzeR`, `clusterProfiler`, etc.).
- **Ready-to-teach content**: notebooks, figures, and data mirrored from the original Gitpod-based course.

## Start in GitHub Codespaces

1. Fork or import this repo (rename it `codespaces_RStudio` if you want to mirror the upstream plan).
2. Click **Code -> Codespaces -> Create codespace on main**.
3. Codespaces automatically runs the `.devcontainer/devcontainer.json` stack. The first build can take ~8-10 minutes because R, Bioconductor, and the rnaseq-tools mamba env are installed.
4. When the ports tab lists **8787**, open it in the browser and log into RStudio Server with username `rstudio` and password `rstudio`.
5. If you prefer a desktop workflow, open **port 6080** (password `rstudio`). That launches the XFCE session provided by desktop-lite; from there you can start a browser and navigate to `http://localhost:8787` or launch any GUI apps you add later.

## Local VS Code / Dev Container usage

```bash
# Requires Docker + VS Code Dev Containers extension
git clone https://github.com/FairTeach/codespaces_RStudio.git
cd codespaces_RStudio
code .
# When prompted, "Reopen in Container"
```
All devcontainer features behave the same locally because GitHub Codespaces and VS Code Dev Containers share the identical runtime model.

## Managing CLI dependencies with mamba

The Miniforge feature installs mamba into `/opt/conda`. The bootstrap script automatically creates/updates the `rnaseq-tools` environment defined in `env/mamba-environment.yml`. Use it like this once the container is up:

```bash
mamba activate rnaseq-tools
multiqc --help
```

Feel free to add more tools (e.g., `gffcompare`, `bwa`, `hisat2`) to that YAML file; rebuild or run `mamba env update -f env/mamba-environment.yml` to propagate changes.

## Pre-installed R packages

The post-create hook runs `.devcontainer/scripts/install-r-packages.R`, ensuring the following are ready on first launch:

- airway, SummarizedExperiment, DESeq2, vsn, apeglm, AnnotationDbi, org.Hs.eg.db
- clusterProfiler, enrichplot, DOSE, IsoformSwitchAnalyzeR
- tidyverse, ggplot2, ggnewscale, ggbeeswarm, ggrepel, pheatmap, RColorBrewer, dbplyr, tidyr, readr, GGally

Add more in that script if your lessons need extra libraries.

## Repository layout

- `docs/` - grading material and supporting markdown.
- `env/` - shared config such as `mamba-environment.yml` and R data sets.
- `figures/` - screenshots used in the exercises.
- `rocker/` - legacy Dockerfiles kept for reference.
- `*.Rmd` - the teaching notebooks students will open inside RStudio.

## Troubleshooting tips

- Rebuild the container (`Codespaces / Rebuild container`) whenever you change `.devcontainer/**`, `env/mamba-environment.yml`, or the bootstrap scripts.
- If RStudio Server does not start, check the **Ports** tab to confirm 8787 is forwarded. The service runs under `systemctl --user` managed by the `rstudio-server` feature, so restarts are persistent across Codespace suspends.
- Desktop-lite sessions can be reset from the Command Palette (`Codespaces: Open in Browser`) if the VNC view freezes.

Happy teaching!
