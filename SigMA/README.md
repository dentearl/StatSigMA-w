# SigMA
(c) 2007-2012 The Authors, see LICENSE for details.

## Authors
Xiaoyu Chen, Amol Prakash, [Martin Tompa](http://www.cs.washington.edu/homes/tompa/). All of University of Washington.

## Adaptor
[Dent Earl](https://github.com/dentearl/) of UC Santa Cruz. Dent Earl adapted the StatSigMA-w code to allow parameters to be changed on the command line and to add input file verification steps.

## SigMA
### Function
Given a branch of the phylogenetic tree and a branch multiplier, SigMA identifies alignment segments w.r.t. the given branch and assigns a p-value to each segment. 

In order to perform the downstream procedure of StatSigMAw, SigMA should be run on the combinations of ALL possible branches of the phylogenetic tree and THREE branch multipliers (0.01, 1 and 100).


### Usage
    SigMA <path/maf file> <branch index> <branch multiplier> <path/multi_segment_pvalue.txt> [options]

    maf file: the multiple sequence alignment in maf format.
    branch index: the index of one branch of the phylogenetic tree corresponding to the alignment. 
    branch multiplier: the multiplier of branch length of the phylogenetic tree.
    pvalue.txt: path to the (supplied) pvalue distribution file.

    Options
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
    --maxBlockSize=MBS        maximum total size of contiguous alignment blocks allowed.
                              If the length of one alignment block is larger than value assigned
                              maxBlockSize will be assigned to the actual alignment block size. default=100000
    --maxSegments=MS          max number of high-scoring segments. default=100000
    
    These are used to estimate Karlin-Altschul parameters:
    --totalNumberTuples=TNT   max number of tuples to be read per block (after randomization). default=50000
    --totalIterateParam=TIP   number of times that Karlin-Altschul parameter is computed. default=100

NOTICE: We need to run SigMA on the combinations of ALL possible branches of the tree and THREE branch multipliers (0.01, 1 and 100). For example, the default tree provided in our code has 53 branches, and the branch index ranges between 0-52.  We will then run SigMA 53X3 times.

### Output
The program generates the output to stdout, which includes information about identifed alignment segments, such as p-values. 

   The stdout should be REDIRECTED to a file whose PREFIX has to be "<branch index>_<branch multiplier>.out".

### Example
    SigMA  myAlign.maf  6  0.01  >  myOutDir/myPrefix.6_0.01.out 
    
We show an example of running SigMA on branch #6 and multipler 0.01. Notice the output file has a prefix "myPrefix".
