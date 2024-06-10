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
 - Update plots: edit `wec_as_multiport.ipynb`, output figures as PDFs to `gfx` directory

## Structure

 - `gfx`: figures (preferably vector PDFs) for the paper
 - `paper`: LaTeX source (including BibTex file)
 - `src`: Python and MATLAB code
