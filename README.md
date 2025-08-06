# WEC as a multi-port
Modeling wave energy converters (WECs) using a multi-port network framework

## Goals
 - Consistent and correct conventions
 - More explanation, details, and examples than previously published papers
 - Publish this repository with the paper

## Structure
`wec_as_multiport`
├── `tests`: testing using [`pytest`](https://docs.pytest.org/en/stable/)
├── `papers`: LaTeX source and plotting scripts
└── `wec_as_multiport`: Python package

## LaTeX
This repo uses actions to compile PDFs from the LaTex source and deploy those PDFs to its [build branch](/../build)

 - Edit the paper: edit the LaTex source (do not commit the PDF to the git repo, it is compiled by GitHub Automation)
 	- use "`XX`" to make a comment
 	- line break after each sentence
 	- use "`\,`" for a small space between a value and unit (e.g., "`5\,m`")
 	- to refer to a figure, use “`\figurename~\ref{fig:fig_label}`”
 	- to refer to a table, use “`Table~\ref{tab:tab_label}`”
 	- use the [subequations environment](https://tex.stackexchange.com/questions/38996/referencing-main-subequation) when possible if you have multiple equations you're presenting together
 	- plots should not generally have titles as this information is captured in the figure caption
 	- use PDF vector figures when possible
 - Edit diagrams: edit the PDF files using [IPE](https://ipe.otfried.org/)
 - Update plots:
 	1. Edit python notebooks
 	2. Output vector PDFs to `gfx` directory

## Python
 1. Install `wec_as_multiport` package (assuming you have [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html))
	1. Create environment: `conda create -n wam pip "python=3.11"`
	2. Activate environment: `conda activate wam`
	3. Install package in editable mode (from within the root directory of this repository): `pip install -e .`
 2. Run/edit 
	- source: [`core.py`](wec_as_multiport/core.py)
	- plotting scripts: see `.ipynb` files
 3. Output figures to `gfx` directory for inclusion in paper
