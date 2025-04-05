# Mass Spectrometry imaging and spatially-resolved Transcriptomics Integrator(MSTI)
## Overview

With the advancement of technology, spatial omics is gradually transforming our understanding of the complexity of life. It provides a new perspective, enabling us to explore biological processes in the spatial dimension. Among the commonly used spatial omics approaches are spatially-resolved Transcriptomics and Mass spectrometry imaging; however, the joint analysis of spatially-resolved Transcriptomics and Mass spectrometry imaging remains challenging. we introduce the MSTI (Mass Spectrometry imaging and spatially-resolved Transcriptomics Integrator) framework, which utilizes a KEGG-guided graph to correlate features across the spatially-resolved Transcriptomics and Mass spectrometry imaging datasets. MSTI employs a Variational Graph Auto-Encoders (VGAE) for dimensionality reduction and a LouvainBayesian clustering strategy to achieve accurate co-clustering of multi-omics data , facilitating the integration and correlation of molecular features across these layers[1]. In the MSTI framework, we accomplished cross-sample and cross-patient integration of spatial Transcriptomics and Mass spectrometry imaging data.

[1]	KIPF T N, WELLING M. Variational Graph Auto-Encoders[J/OL] 2016, arXiv:1611.07308[https://ui.adsabs.harvard.edu/abs/2016arXiv161107308K.


## Quick Start Guide for MSTI for SRT data and MSI data

### 🚀 Quick Start Guide for LouvainBayesian Analysis for SRT data

We provide a demo R script named **`ST_analysis.R`** and a demo spatial transcriptomics (SRT) dataset named **`ST_outs`** for LouvainBayesian analysis.

Please follow the steps below to run the demo:

1. **Place Files in the Same Directory**  
   Make sure the **`ST_analysis.R`** file and the **`ST_outs`** folder are located in the **same directory**.

2. **Open the Script**  
   Open **`ST_analysis.R`** using RStudio or your preferred R environment.

3. **Edit Parameters**  
   Modify any relevant parameters in the script as needed (e.g., input path, resolution settings, etc.).

4. **Run the Code**  
   Select all the code (`Ctrl + A`), then click the **"Run"** button to execute the analysis.

### 🚀 Quick Start Guide for LouvainBayesian Analysis for MSI data

We provide a demo R script named **`SM_analysis.R`** and a demo spatial transcriptomics (SRT) dataset named **`SM_outs`** for LouvainBayesian analysis.

Please follow the steps below to run the demo:

1. **Place Files in the Same Directory**  
   Make sure the **`SM_analysis.R`** file and the **`SM_outs`** folder are located in the **same directory**.

2. **Open the Script**  
   Open **`SM_analysis.R`** using RStudio or your preferred R environment.

3. **Edit Parameters**  
   Modify any relevant parameters in the script as needed (e.g., input path, resolution settings, etc.).

4. **Run the Code**  
   Select all the code (`Ctrl + A`), then click the **"Run"** button to execute the analysis.

### 🧪 SRT&MSI Integrator

We provide demo folders for the **ST&SM Integrator**, which can be used to perform integrated analysis of spatial transcriptomics and spatial metabolomics data.  
Currently, two datasets are available: **DCIS** and **IBC**.

Please follow the steps below to run either demo:

#### 📁 Folder Structure
- `MSTI/ST&SM Integrator/workflow/DCIS/`
- `MSTI/ST&SM Integrator/workflow/IBC/`

Each folder contains:
- A Jupyter Notebook (`.ipynb`)
- A model file (`.dill`)
- A dataset file (`.h5ad`)

---

#### ▶️ Steps to Run the Demo

1. **Place Files in the Same Directory**  
   Ensure the following three files are located in the **same folder** (either `DCIS` or `IBC`):
   - `<dataset>.ipynb`
   - `<dataset>.dill`
   - `<dataset>.h5ad`

2. **Open the Notebook**  
   Launch the notebook (e.g., `DCIS.ipynb` or `IBC.ipynb`) using **Anaconda**, via **Jupyter Notebook** or **JupyterLab**.

3. **Edit Parameters**  
   Modify any relevant parameters in the notebook as needed (e.g., input paths, resolution settings, etc.).

4. **Run the Code**  
   Select all cells (`Ctrl + A`), then click the **"Run"** button (or use **`Cell > Run All`**) to execute the analysis.
---

If you encounter any issues, please refer to the [Troubleshooting](#troubleshooting) section or open an issue in the repository.

## Requirements
All benchmark tests were performed on a personal computer with 10th Gen Intel® Core™ i5-1035G1- Core Processor, 8GB memory, and installed with Windows10 operation system , R-4.3.1 and Python v. 3.11.11 including packages CARD (v1.1), clustree (v0.5.0), confuns (v1.0.3), cowplot (v1.1.1), data.table (v1.14.8), distances (v0.1.11), dplyr (v1.1.4), future (v1.33.0), ggplot2 (v3.5.1), ggpubr (v0.6.0), ggsci (v3.0.0), GSEABase (v1.62.0), gson (v0.1.0), GSVA (v1.48.3), mistyR (v1.10.0), monocle (v2.26.0), msigdbr (v7.5.1), pals (v1.8), patchwork (v1.2.0), RColorBrewer (v1.1-3), recipes (v1.0.8), Seurat (v4.4.0), spacexr (v2.2.1), SPATA2 (v2.0.4), STdeconvolve (v1.3.1), stringr (v1.5.0), tibble (v3.2.1), tidyr (v1.3.0), tidyverse (v2.0.0), anndata (v0.9.2), anyio (v4.2.0), argon2-cffi (v21.3.0), argon2-cffi-bindings (v21.2.0), asttokens (v2.0.5), async-lru (v2.0.4), attrs (v24.2.0), babel (v2.11.0), backcall (v0.2.0), beautifulsoup4 (v4.12.3), blas (v1.0), bleach (v4.1.0), brotli-python (v1.0.9), ca-certificates (v2024.11.26), certifi (v2024.8.30), cffi (v1.17.1), charset-normalizer (v3.3.2), colorama (v0.4.6), comm (v0.2.1), contourpy (v1.1.1), cycler (v0.12.1), debugpy (v1.6.7), decorator (v5.1.1), defusedxml (v0.7.1), dill (v0.3.9), et_xmlfile (v1.1.0), exceptiongroup (v1.2.0), executing (v0.8.3), filelock (v3.16.1), fonttools (v4.54.1), fsspec (v2024.10.0), get-annotations (v0.1.2), h11 (v0.14.0), h5py (v3.11.0), httpcore (v1.0.2), httpx (v0.27.0), icc_rt (v2022.1.0), icu (v73.1), idna (v3.7), igraph (v0.11.8), importlib-metadata (v8.5.0), importlib-resources (v6.4.5), importlib_metadata (v7.0.1), importlib_resources (v6.4.0), intel-openmp (v2021.4.0), ipykernel (v6.29.5), ipython (v8.12.2), ipywidgets (v8.1.2), jedi (v0.19.1), jinja2 (v3.1.4), joblib (v1.4.2), jpeg (v9e), json5 (v0.9.6), jsonschema (v4.23.0), jsonschema-specifications (v2023.7.1), jupyter (v1.0.0), jupyter-lsp (v2.2.0), jupyter_client (v8.6.0), jupyter_console (v6.6.3), jupyter_core (v5.7.2), jupyter_events (v0.10.0), jupyter_server (v2.14.1), jupyter_server_terminals (v0.4.4), jupyterlab (v4.2.5), jupyterlab_pygments (v0.1.2), jupyterlab_server (v2.27.3), jupyterlab_widgets (v3.0.10), kiwisolver (v1.4.7), krb5 (v1.20.1), leidenalg (v0.10.2), libclang (v14.0.6), libclang13 (v14.0.6), libffi (v3.4.4), libpng (v1.6.39), libpq (v12.20), libsodium (v1.0.18), llvmlite (v0.41.1), louvain (v0.8.2), lz4-c (v1.9.4), markupsafe (v2.1.5), matplotlib (v3.7.5), matplotlib-inline (v0.1.6), mistune (v2.0.4), mkl (v2021.4.0), mkl-service (v2.4.0), mkl_fft (v1.3.1), mkl_random (v1.2.2), mpmath (v1.3.0), natsort (v8.4.0), nbclient (v0.8.0), nbconvert (v7.16.4), nbformat (v5.10.4), nest-asyncio (v1.6.0), networkx (v3.1), notebook (v7.2.2), notebook-shim (v0.2.3), numba (v0.58.1), numpy (v1.24.4), numpy-base (v1.24.3), openpyxl (v3.1.5), openssl (v3.0.15), overrides (v7.4.0), packaging (v24.2), pandas (v2.0.3), pandocfilters (v1.5.0), parse (v1.20.2), parso (v0.8.3), patsy (v1.0.1), pickleshare (v0.7.5), pillow (v10.4.0), pip (v24.2), pkgutil-resolve-name (v1.3.10), platformdirs (v3.10.0), ply (v3.11), prometheus_client (v0.14.1), prompt-toolkit (v3.0.43), protobuf (v5.28.3), psutil (v5.9.0), pure_eval (v0.2.2), pycparser (v2.21), pygments (v2.15.1), pynndescent (v0.5.13), pynvml (v11.5.3), pyparsing (v3.1.4), pyqt (v5.15.10), pyqt5-sip (v12.13.0), pysocks (v1.7.1), python-dateutil (v2.9.0), python-fastjsonschema (v2.16.2), python-json-logger (v2.0.7), pytorch-ignite (v0.5.1), pytz (v2024.2), pywin32 (v305), pywinpty (v2.0.10), pyyaml (v6.0.2), pyzmq (v25.1.2), qt-main (v5.15.2), qtconsole (v5.6.0), qtpy (v2.4.1), referencing (v0.30.2), requests (v2.32.3), rfc3339-validator (v0.1.4), rfc3986-validator (v0.1.1), rpds-py (v0.10.6), scanpy (v1.9.8), scglue (v0.3.2), scikit-learn (v1.3.2), scipy (v1.10.1), seaborn (v0.13.2), session-info (v1.0.0), setuptools (v75.1.0).

* Operating System: Windows
* Compilation Environment
~~~
Packages need in R code.
devtools::load_all(".\\LouvainBayesian")
devtools::load_all(".\\monocle")
library(Seurat)
library(SPATA2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(msigdbr)
library(cowplot)
library(clustree)
library(confuns)
library(ggsci)
library(GSVA)
library(gson)
library(GSEABase)
library(data.table)
library(spacexr)
library(tidyverse)
library(stringr)
library(patchwork)
library(STdeconvolve)
library(pals)
library(mistyR)
library(CARD)
library(distances)
library(future)
library(recipes)
~~~
~~~
Requirements need in python code.
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import scanpy as sc
from itertools import chain
import itertools
import networkx as nx
import scglue
import seaborn as sns
from matplotlib import rcParams
import sys
import os
import MSTI_Guidancegraph as MG
import MSTI_preprocess as MP
~~~
* Hardware Requirements: no special requires.

* Estimated Time for Installing Required Packages
~~~
The R code for installing packages takes approximately 10 minutes, depending on individual network conditions.
The python code for installing requirements takes approximately 10 minutes, depending on individual network conditions.
~~~
