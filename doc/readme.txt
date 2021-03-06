I. StatSigMAw v2 

    Given a multiple sequence alignment and a phylogeny relating the aligned sequences, StatSigMAw assesses every site of the alignment and identifies suspiciously aligned regions in the alignment. 

    This package has two major programs: SigMA and combine. We also provide a supplementary program, getSpeciesName. The function and usage of these programs are described in execute.txt.


II. How to install

    Information about installing the programs is listed in install.txt.


III. Command lines and outputs

    Detailed help on function, usage, command line and output for each program are described in execute.txt.


IV. More about parameter setting

    There is more information in parameters.txt about advanced parameter setting.


V. Running time

    The running time of StatSigMAw is determined by the number of aligned species and the length of aligned sequences. Using the default tree provided in our code, we give two examples below to illustrate the running time.

    For ENCODE region ENm003 (0.5 Mbp), the running time of SigMA on all 53 branches is 15-17 hours (17-19 minutes per branch on average), depending on the branch multiplier; the running time of combine is 18 minutes. 

    For ENCODE region ENm002 (1 Mbp), the running time of SigMA on all 53 branches is 23-24 hours (26-27 minutes per branch on average), depending on the branch multiplier; the running time of combine is 35 minutes.

    The computer on which we ran these two examples has a 2.33GHz CPU and a 4G memeory. 


VI. Conditions of Use

   Xiaoyu Chen, Amol Prakash, and Martin Tompa give permission for Dent Earl and his institution to use the StatSigmMA-w software developed at the University of Washington for research purposes, on the following conditions:

   1. Dent Earl and his institutional colleagues may modify the StatSigmMA-w software and distribute the resulting modified software for research purposes, provided (a) that it is distributed together with these Conditions of Use, (b) that Martin Tompa receive a copy of the finalized modified software, and (c) that Xiaoyu Chen, Amol Prakash, and Martin Tompa are credited with the authorship of the software.

   2. The StatSigmMA-w software will be used by you and/or your institution solely for noncommercial purposes, except with express permission from the authors.

   3. Any risk associated with using the StatSigmMA-w software at your institution is with you and your institution.

   4. StatSigmMA-w will be cited in any publication(s) reporting on data obtained from it as:
      Amol Prakash and Martin Tompa, "Measuring the Accuracy of Genome-Size Multiple Alignments". Genome Biology, vol. 8, issue 6, June 2007, R124.
      Xiaoyu Chen and Martin Tompa, "Comparative assessment of methods for aligning multiple genome sequences". Nature Biotechnology, vol. 28, no. 6, June 2010, 567-572.
