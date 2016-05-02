# diffloop
[![Travis-CI Build Status](https://travis-ci.org/aryeelab/diffloop.svg?branch=master)](https://travis-ci.org/aryeelab/diffloop)
An R package for identifying differential features from ChIA-PET expierments. 

## Preprocessing
Raw FASTQ read files can be preprocessed with the `dnaloop`
(https://github.com/aryeelab/dnaloop) software, a `PyPi` package that relies 
on `samtools`, `bedtools`, `cutadapt`, and `MACS2` to align reads,
call anchor peaks, and summarize PETs per sample in the ChIA-PET experiment.
<br> <br>
To use the `looptestMake` function from a different preprocessing step,
have files `X.loop_counts.bedpe`,`Y.loop_counts.bedpe`, 
`Z.loop_counts.bedpe` in `bed_dir` for `samples = (X,Y,Z)` where the 
first lines should resemble: <br> <br>
1 10002272 10004045 10 120968807 120969483 . 1 <br>
1 10002272 10004045 10 99551498 99552470 . 1 <br>
1 10002272 10004045 1 10002272 10004045 . 17 <br>
<br> <br>
where the first three columns specify the position (chr:start:stop) of
the first anchor, the second three columns specify the position
(chr:start:stop) of the secnd anchor, the 7th column is "." and the 8th
column is the number of paired-end reads support that particular PET.

## Loading data
We read in data from the preprocessing pipeline using the `loopdataMake`
function. The example below uses sample data included in the `diffloopdata`
package. It contains interactions involving chromosome 1 from a set of six
Cohesin ChIA-PET samples.$^{1,2}$ You would normally set `bed_dir` to your 
real data directory.
