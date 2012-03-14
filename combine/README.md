# combine and getSpeciesName
(c) 2007-2012 The Authors, see LICENSE for details.

## Authors
Xiaoyu Chen, Amol Prakash, [Martin Tompa](http://www.cs.washington.edu/homes/tompa/). All of University of Washington.

## Adaptor
[Dent Earl](https://github.com/dentearl/) of UC Santa Cruz. Dent Earl adapted the StatSigMA-w code to allow parameters to be changed on the command line and to add input file verification steps.

## Combine
### Function
Given the directory where the SigMA output files are located and the prefix of SigMA output files, the alignment segments identified by SigMA for all branches in the phylogenetic tree and for all three branch multipliers are assessed. The program finally identifies suspiciously aligned regions with respect to each branch.

NOTICE: This program assumes that SigMA has been run on the combinations of ALL possible branches of the phylogenetic tree and THREE branch multipliers (0.01, 1 and 100).
### Usage
    combine <maf file> <SigMA output dir> <SigMA output prefix> [options]

* <maf file>: the multiple sequence alignment in maf format.
* <SigMA output dir>: the directory where the SigMA output files are located.
* <SigMA output prefix>: the prefix of SigMA output file names.
    -h, --help                prints a help message and exits
    --phylogeny=PHYLOGENY     phylogeny in newick format default=(((((((((((((hg:0.006690,chimp:0.007571):0.024272,
                              (colobus_monkey:0.015404,(baboon:0.008258,macaque:0.028617):0.008519):0.022120)
                              :0.023960,(dusky_titi:0.025662,(owl_monkey:0.012151,marmoset:0.029549):0.008236)
                              :0.027158):0.066101,(mouse_lemur:0.059024,galago:0.121375):0.032386):0.017073,
                              ((rat:0.081728,mouse:0.077017):0.229273,rabbit:0.206767):0.023340):0.023026,
                              (((cow:0.159182,dog:0.147731):0.004946,rfbat:0.138877):0.010150,(hedgehog:0.193396,
                              shrew:0.261724):0.054246):0.024354):0.028505,armadillo:0.149862):0.015994,
                              (elephant:0.104891,tenrec:0.259797):0.040371):0.218400,monodelphis:0.371073)
                              :0.065268,platypus:0.468116):0.123856,chicken:0.454691):0.123297,xenopus:0.782453)
                              :0.156067,((tetraodon:0.199381,fugu:0.239894):0.492961,zebrafish:0.782561):0.156067)
    --refSpecies=REFSPECIES   reference species in phylogeny default=hg
    --minSegSize=MINSEGSIZE   specifies the minimum size of suspicious or good regions.
                              default=50
    --pThreshBad=PTHRESHBAD   specifies the p-value threshold for suspicious regions. Every 
                              position in a suspicious region has a p-value greater than the 
                              specified value. default=0.1
    --pThreshGood=PTHRESHGOOD specifies the p-value threshold for good regions. Every position 
                              in a good region has a p-value less than the specified value.
                              default=1e-10

### Output
The program generates the output to stdout. The output consists of three parts:
a) Regions that are well aligned with respect to ALL branches. The format of this part of output is:
- my_good_region (simply a flag) 
- start position of this region (in the alignment coordinate, i.e. which column in the alignment is this position)
- length of this region
- chromosome
- start position 
- end position

The last three items show the genomic coordinates of the reference species for this region.

b) Regions that are suspiciously aligned w.r.t. a branch. The format of this part of output is:
- my_bad_region (simply a flag) 
- start position of this region (in the alignment coordinate, i.e. which column in the alignment is this position)
- branch index (the region is suspiciously aligned with respect to this particular branch)
- length of this region
- chromosome
- start position 
- end position

Again, the last three items show the genomic coordinates of the reference species for this region.

Note: For regions that are suspiciously aligned with respect to a branch incident on a leaf (i.e., a species), you can read them as regions where that species is suspiciously aligned to the other species.

c) Statistics for suspiciously-aligned regions. The format of this part of output is:
- my_count (simply a flag)
- branch index
- number of aligned bases w.r.t. the branch
- number of suspiciously aligned bases w.r.t. the branch
- % of suspiciously aligned bases w.r.t. the branch
- number of suspiciously aligned regions w.r.t. the branch

NOTICE: Please refer to parameters.txt for parameter settings of well- or suspiciously-aligned regions.

## getSpeciesName
### Function
The program generates a mapping between branch indices and species names. For each of the branches incident on leaves, it tells which species that branch is incident on.

### Usage
    getSpeciesName [options]
    
    -h,--help               prints this message and exits
    --phylogeny=PHYLOGENY   phylogeny in newick format default=(((((((((((((hg:0.006690,chimp:0.007571):0.024272,
                            (colobus_monkey:0.015404,(baboon:0.008258,macaque:0.028617):0.008519):0.022120)
                            :0.023960,(dusky_titi:0.025662,(owl_monkey:0.012151,marmoset:0.029549):0.008236)
                            :0.027158):0.066101,(mouse_lemur:0.059024,galago:0.121375):0.032386):0.017073,
                            ((rat:0.081728,mouse:0.077017):0.229273,rabbit:0.206767):0.023340):0.023026,
                            (((cow:0.159182,dog:0.147731):0.004946,rfbat:0.138877):0.010150,(hedgehog:0.193396,
                            shrew:0.261724):0.054246):0.024354):0.028505,armadillo:0.149862):0.015994,
                            (elephant:0.104891,tenrec:0.259797):0.040371):0.218400,monodelphis:0.371073)
                            :0.065268,platypus:0.468116):0.123856,chicken:0.454691):0.123297,xenopus:0.782453)
                            :0.156067,((tetraodon:0.199381,fugu:0.239894):0.492961,zebrafish:0.782561):0.156067)
    --printAll              print every branch in the tree

### Output
The program generates the output to stdout. The output lists the mapping between branch indices and species names. For example, given the default tree provided in our code, branch #6 correponds to platypus, which means the branch incident on platypus has index 6.
