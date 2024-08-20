# WEC as a multi-port

## Goals

 - Consistent and correct conventions
 - More explanation, details, and examples than previously published papers
 - (Maybe?) publish this repository with the paper

## Instructions

 - See the PDF: https://github.com/ryancoe/wec_as_multiport/blob/build/wec_as_multiport.pdf
 - Edit the paper: edit the LaTex source (note that you don't need to commit the PDF to the git repo, it is compiled by GitHub Automation)
 	- use "`XX`" to make a comment
 	- line break after each sentence
 	- use "`\,`" for a small space between a value and unit (e.g., "`5\,m`")
 	- to refer to a figure, use “`\figurename~\ref{fig:fig_label}`”
 	- to refer to a table, use “`Table~\ref{tab:tab_label}`”
 	- use the subequations environment when possible if you have multiple equations you're presenting together (see, e.g., https://tex.stackexchange.com/questions/38996/referencing-main-subequation)
 	- plots should not generally have titles as this information is captured in the figure caption
 	- use PDF vector figures when possible
 - Edit diagrams: edit the PDF files using [IPE](https://ipe.otfried.org/)
 - Update plots: 
 	1. Install package (assuming you have [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html))
 		1. Create an environment: `conda create -n wam pip "python=3.11"`
 		2. Activate environment: `conda activate wam`
 		3. Install package in editable mode (from within the root directory of this repository): `pip install -e .`
 	2. Run/edit 
 		- plotting script: [`wec_as_multiport.ipynb`](wec_as_multiport.ipynb)
 		- source: [`core.py`](wec_as_multiport/core.py)
 	3. Output figures to `gfx` directory for inclusion in paper

## Structure

 - `gfx`: figures (preferably vector PDFs) for the paper
 - `paper`: LaTeX source (including BibTex file)
 - `wec_as_multiport`: Python package
