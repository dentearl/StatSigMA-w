# StatSigMA-w
(c) 2007-2012 The Authors, see LICENSE for details.

## Authors
Xiaoyu Chen, Amol Prakash, Martin Tompa. All of University of Washington.

## Adaptor
[Dent Earl](https://github.com/dentearl/) of UC Santa Cruz. Dent Earl adapted the StatSigMA-w code to allow parameters to be changed on the command line and to add input file verification steps.

## I. StatSigMAw v2 

Given a multiple sequence alignment and a phylogeny relating the aligned sequences, StatSigMAw assesses every site of the alignment and identifies suspiciously aligned regions in the alignment. 

This package has two major programs: SigMA and combine. We also provide a supplementary program, getSpeciesName. The function and usage of these programs are described in doc/execute.txt.

## II. How to install

1. Download or clone the repo.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## III. Command lines and outputs

Detailed help on function, usage, command line and output for each program are described in doc/execute.txt.

Assumptions on input alignments:
* The reference species must be present in every block.
* There may be no duplicate species in any block. A species may be present only once per block.
* The blocks must be ordered such that there is a monotonic increase in position relative to the reference sequence.

## IV. More about parameter setting

There is more information in doc/parameters.txt about advanced parameter setting.


## V. Running time

The running time of StatSigMAw is determined by the number of aligned species and the length of aligned sequences. Using the default tree provided in our code, we give two examples below to illustrate the running time.

For ENCODE region ENm003 (0.5 Mbp), the running time of SigMA on all 53 branches is 15-17 hours (17-19 minutes per branch on average), depending on the branch multiplier; the running time of combine is 18 minutes. 

For ENCODE region ENm002 (1 Mbp), the running time of SigMA on all 53 branches is 23-24 hours (26-27 minutes per branch on average), depending on the branch multiplier; the running time of combine is 35 minutes.

The computer on which we ran these two examples has a 2.33GHz CPU and a 4G memeory. 

## VI. Conditions of Use

Xiaoyu Chen, Amol Prakash, and Martin Tompa give permission for Dent Earl and his institution to use the StatSigmMA-w software developed at the University of Washington for research purposes, on the following conditions:

1. [Dent Earl](https://github.com/dentearl/) and his institutional colleagues may modify the StatSigmMA-w software and distribute the resulting modified software for research purposes, provided (a) that it is distributed together with these Conditions of Use, (b) that Martin Tompa receive a copy of the finalized modified software, and (c) that Xiaoyu Chen, Amol Prakash, and Martin Tompa are credited with the authorship of the software.

2. The StatSigmMA-w software will be used by you and/or your institution solely for noncommercial purposes, except with express permission from the authors.

3. Any risk associated with using the StatSigmMA-w software at your institution is with you and your institution.

4. StatSigmMA-w will be cited in any publication(s) reporting on data obtained from it as:
  * Amol Prakash and Martin Tompa, "Measuring the Accuracy of Genome-Size Multiple Alignments". Genome Biology, vol. 8, issue 6, June 2007, R124.
  * Xiaoyu Chen and Martin Tompa, "Comparative assessment of methods for aligning multiple genome sequences". Nature Biotechnology, vol. 28, no. 6, June 2010, 567-572.

## References
* Amol Prakash and Martin Tompa, "Measuring the Accuracy of Genome-Size Multiple Alignments". Genome Biology, vol. 8, issue 6, June 2007, R124. doi:10.1186/gb-2007-8-6-r124 http://genomebiology.com/2007/8/6/r124
* Xiaoyu Chen and Martin Tompa, "Comparative assessment of methods for aligning multiple genome sequences". Nature Biotechnology, vol. 28, no. 6, June 2010, 567-572. doi:10.1038/nbt.1637 http://www.nature.com/nbt/journal/v28/n6/abs/nbt.1637.html
