# diffloop

[![Travis-CI Build Status](https://travis-ci.org/aryeelab/diffloop.svg?branch=master)](https://travis-ci.org/aryeelab/diffloop)

An R package for identifying differential features from chromatin capture confirmation experiments. 

## Quick Start
Raw FASTQ read files can be preprocessed with the `mango`
(https://github.com/dphansti/mango) software, an `R` package that relies 
on `bowtie`, `bedtools`, and `MACS2` to align reads,
call anchor peaks, and summarize PETs per sample in the ChIA-PET experiment.
<br> <br>
We recommend setting the `reportallpairs` flag in Mango to `TRUE` to maintain the
integrity of differential loop calling (non-zero interactions may be discarded
otherwise). The `mangoCorrection()` function can then be applied to a `loops`
object in the R environment by aggregating loop data across samples after
the data are imported using the `loopsMake.mango` function. Other preprocessing
software may be appropriate and can more generally be imported using the `loopsMake`
function, which reads in `.bedpe` files. 

## Example analyses with diffloop
Check out the [vignette](https://dl.dropboxusercontent.com/u/210183/diffloop_2016-04-13.html)
for a working example on quality control, visualization, association, and annotation.

## Contact
[Caleb Lareau](https://caleblareau.github.io) is the primary maintainer/developer of diffloop. 