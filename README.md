# diffloop

[![Travis-CI Build Status](https://travis-ci.org/aryeelab/diffloop.svg?branch=master)](https://travis-ci.org/aryeelab/diffloop)

An R package for identifying differential features from chromatin capture confirmation experiments. 

## Quick Start
Raw FASTQ read files can be preprocessed with the `mango`
(https://github.com/dphansti/mango) software, an `R` package that relies 
on `bowtie`, `bedtools`, and `MACS2` to align reads,
call anchor peaks, and summarize PETs per sample in the ChIA-PET experiment.

## Sample analysis with diffloop
Check out the [vignette](http://rpubs.com/caleblareau/diffloop_vignette)
for a working example on quality control, visualization, association, and annotation.

## Interactive visualization
Check out [DNAlandscapeR](dnalandscaper.aryeelab.org) for an R/Shiny environment to
visualize topology data interactively. 

## Contact
[Caleb Lareau](https://caleblareau.github.io) is the primary maintainer/developer of diffloop. 