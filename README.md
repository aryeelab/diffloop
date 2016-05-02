# diffloop

[![Travis-CI Build Status](https://travis-ci.org/aryeelab/diffloop.svg?branch=master)](https://travis-ci.org/aryeelab/diffloop)

An R package for identifying differential features from ChIA-PET expierments. 

## Preprocessing
Raw FASTQ read files can be preprocessed with the `dnaloop`
(https://github.com/aryeelab/dnaloop) software, a `PyPi` package that relies 
on `samtools`, `bedtools`, `cutadapt`, and `MACS2` to align reads,
call anchor peaks, and summarize PETs per sample in the ChIA-PET experiment.
<br> <br>
To use the `loopdataMake` function from a different preprocessing step,
have files `X.loop_counts.bedpe`,`Y.loop_counts.bedpe`, 
`Z.loop_counts.bedpe` in `bed_dir` for `samples = (X,Y,Z)` where the 
first lines should resemble:
```
1 10002272 10004045 10 120968807 120969483 . 1
1 10002272 10004045 10 99551498  99552470  . 1
1 10002272 10004045 1  10002272  10004045  . 17
```

where the first three columns specify the position (chr:start:stop) of
the first anchor, the second three columns specify the position
(chr:start:stop) of the secnd anchor, the 7th column is "." and the 8th
column is the number of paired-end reads support that particular PET.

## Loading data
We read in data from the preprocessing pipeline using the `loopdataMake`
function. The example in the [vignette](https://dl.dropboxusercontent.com/u/210183/diffloop_2016-04-13.html) uses sample data included in the `diffloopdata`
package. 

## Analysis
Check out the [vignette](https://dl.dropboxusercontent.com/u/210183/diffloop_2016-04-13.html) for a working example on quality control, visualization, 
association, and annotation
