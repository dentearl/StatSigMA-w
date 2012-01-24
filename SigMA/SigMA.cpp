/*VI. Conditions of Use
 *
 * Xiaoyu Chen, Amol Prakash, and Martin Tompa give permission for Dent Earl and his institution to use the StatSigmMA-w software developed at the University of Washington for research purposes, on the following conditions:
 *
 * 1. Dent Earl and his institutional colleagues may modify the StatSigmMA-w software and distribute the resulting modified software for research purposes, provided (a) that it is distributed  together with these Conditions of Use, (b) that Martin Tompa receive a copy of the finalized modified software, and (c) that Xiaoyu Chen, Amol Prakash, and Martin Tompa are credited with the authorship of the software.
 *  2. The StatSigmMA-w software will be used by you and/or your institution solely for noncommercial purposes, except with express permission from the authors.
 * 3. Any risk associated with using the StatSigmMA-w software at your institution is with you and your institution.
 * 4. StatSigmMA-w will be cited in any publication(s) reporting on data obtained from it as:
 *    Amol Prakash and Martin Tompa, "Measuring the Accuracy of Genome-Size Multiple Alignments". Genome Biology, vol. 8, issue 6, June 2007, R124.
 *    Xiaoyu Chen and Martin Tompa, "Comparative assessment of methods for aligning multiple genome sequences". Nature Biotechnology, vol. 28, no. 6, June 2010, 567-572.
 */

#include <getopt.h>
#include <iostream>
#include <math.h>
#include "alpha.h"
#include "block.h"
#include "global.h"
#include "tree.h"

Tree* root;
Tree* list_of_leaves[MAX_NUM_SPECIES];

#include "pvalue.h"

#define length(x) (sizeof(x) / sizeof(x[0]))

// Helper Functions
Tree* read_tree(char newick_format[], int& sp_num, int& pos);
void generate_and_read_newick_tree();
void read_multi_pv(char* filename);
void init_background();
void loadDefaultParameters(void);
void parseArgs(int argc, char **argv);
void usage(void);

// Main program

int main(int ARGC, char* ARGV[]) {
    loadDefaultParameters();
    if(ARGC != 4){
        fprintf(stderr, "Error: wrong number of arguments.\n");
        usage();
    }
    parseArgs(ARGC, ARGV);
    
    BRANCH_MULTIPLIER = atof(ARGV[3]);
    char* pos = ARGV[0];
    while (strstr(pos, "SigMA") != NULL)
        pos = strstr(pos, "SigMA") + 1;
    pos--;
    char binary_path[1000];
    strncpy(binary_path, ARGV[0], pos - ARGV[0]);
    binary_path[int(pos - ARGV[0])] = '\0';

    generate_and_read_newick_tree();

#ifdef DEBUG
    cout << "# Generated Tree " << root->tree_name << " " << "#species = " << species_num << endl;
#endif

    // Read the pvalue distribution of the various scores and number of segments
    char filename[1000];
    strcpy(filename, ARGV[4]);
    // sprintf(filename, "%smulti_segment_pvalue.txt", binary_path);
    read_multi_pv(filename);
#ifdef DEBUG
    cout << "# Read multipv " << filename << endl;
#endif

#ifdef PROTEIN
    sprintf(filename, "%s/pmb.dat", binary_path);
    read_aa_matrix(filename);
#ifdef DEBUG
    cout << "# Read AA matrix " << filename << endl;
#endif
#endif

    init_background();  

    blocks_type blocks;
    /*
      blocks.read_from_stream(ARGV[2]);
      #ifdef DEBUG
      cout << "# Read #blocks=" << blocks.chain.size() << " ; last blocklength= " << blocks.chain[blocks.chain.size()-1]->length << endl;
      #endif
    */

    double pvalue[2 * MAX_NUM_SPECIES];
    Tree* branch_ptr[2 * MAX_NUM_SPECIES];
    int branch=0;
    compute_pvalue(root, blocks, pvalue, branch, branch_ptr, atoi(ARGV[2]), ARGV[1]);

    double max_pvalue = -1;
    Tree* weakest_branch = NULL;
    for (int i = 0; i < branch; i++) {
        if (pvalue[i] > max_pvalue) {
            max_pvalue = pvalue[i];
            weakest_branch = branch_ptr[i];
        }
    }

    cout << "p-value= " << max_pvalue << "\tbranch= " << weakest_branch->tree_name 
         << "\ttree= " << root->tree_name << endl;

    return 0;
}

/*--------------------- Helper Functions ------------------ */

Tree* read_tree(char newick_format[], int& sp_num, int& pos) {
    while (isblank(newick_format[pos]))
        pos++;
    if (newick_format[pos] == '(') {
        pos++;
        Tree* left = read_tree(newick_format, sp_num, pos);
        while (newick_format[pos] != ':')
            pos++;
        pos++;
        while (isblank(newick_format[pos]))
            pos++;
        int start_pos = pos;
        while (newick_format[pos] != ',')
            pos++;
        char branch_length[100];
        strncpy(branch_length, newick_format + start_pos, pos-start_pos);
        branch_length[pos - start_pos] = '\0';
        float left_branch = atof(branch_length);
        if (left_branch < 0.001) left_branch = 0.001;
        left_branch *= BRANCH_MULTIPLIER;
        if (left_branch > MAX_BRANCH_LENGTH) 
            left_branch = MAX_BRANCH_LENGTH;
    
        pos++;
        Tree* right = read_tree(newick_format, sp_num, pos);
        while (newick_format[pos] != ':')
            pos++;
        pos++;
        while (isblank(newick_format[pos]))
            pos++;
        start_pos = pos;
        while (newick_format[pos] != ')')
            pos++;
        strncpy(branch_length, newick_format + start_pos, pos - start_pos);
        branch_length[pos - start_pos] = '\0';
        float right_branch = atof(branch_length);
        if (right_branch < 0.001) right_branch = 0.001;
        right_branch *= BRANCH_MULTIPLIER;
        if (right_branch > MAX_BRANCH_LENGTH) 
            right_branch = MAX_BRANCH_LENGTH; 
    
        pos++;
        return (new Tree(left, right, left_branch, right_branch));
    }else {
        int start_pos = pos;
        while (newick_format[pos] != ':')
            pos++;
        strncpy(species_name[sp_num], newick_format + start_pos, pos - start_pos);
        species_name[sp_num][pos - start_pos] = '\0';    
        Tree* node = new Tree(sp_num, species_name[sp_num]);
        list_of_leaves[sp_num] = node;
        sp_num++;
        return node;
    }
}

void generate_and_read_newick_tree() {

    int pos = 0;
    species_num = 0;
    root = read_tree(globalOptions.PHYLOGENY, species_num, pos);
}

void read_multi_pv(char* filename) {
    multi_pv = new double*[globalOptions.MAX_SEGMENTS + 1];
    for (int i = 0; i < globalOptions.MAX_SEGMENTS + 1; i++)
        multi_pv[i] = new double[int(MAX_TOTAL_SCORE - MIN_TOTAL_SCORE + 1)];

    for (int i = 0; i < (globalOptions.MAX_SEGMENTS + 1); i++)
        for (int j = 0; j < (MAX_TOTAL_SCORE - MIN_TOTAL_SCORE + 1); j++)
            multi_pv[i][j]=-1;
  
    for (int j = MIN_TOTAL_SCORE; j <= MAX_TOTAL_SCORE; j++)
        multi_pv[1][j - MIN_TOTAL_SCORE] = 1 - exp(-1 * exp(-1 * j));
    ifstream ifs(filename);
    int a, b;
    double c;
    if (!ifs) {
        cout << "multi_segment_pvalue.txt not found ... quitting\n";
        exit(1);
    }
    while (ifs) {
        ifs >> a >> b >> c;
        if ((a <= globalOptions.MAX_SEGMENTS) && (b <= MAX_TOTAL_SCORE) && (b >= MIN_TOTAL_SCORE)) {
            //      if (b<0.5) c=1; 
            multi_pv[a][b - MIN_TOTAL_SCORE] = c;
        }
    }
    ifs.close();
}

void init_background() {
    back[0] = 0;
#ifdef PROTEIN
    for (int i = 1; i < ALPHABET_SIZE; i++)
        back[i] = aa_back[i - 1] * (1 - GAP_BACKGRD);
#else
    for (int i = 1; i < ALPHABET_SIZE; i++)
        back[i] = (1 - GAP_BACKGRD - back[0]) / float(ALPHABET_SIZE - 1);
#endif
    back[GAP] = GAP_BACKGRD;
}

void usage(void)
{
    fprintf(stderr, "Usage: SigMA <maf file> <branch index> <branch multiplier> <pvale dist file>[options]\n\n");
    fprintf(stderr, "  <maf file>: the multiple sequence alignment in maf format.\n");
    fprintf(stderr, "  <branch index>: the index of one branch of the phylogenetic tree\n"
            "                  corresponding to the alignment.\n");
    fprintf(stderr, "  <branch multiplier>: the multiplier of branch length of the phylogenetic tree.\n");
    fprintf(stderr, "  <pvale dist file>: pvalue distribution of the various scores and number of segments.\n");

    fprintf(stderr, "\nDESCRIPTION\n");
    fprintf(stderr,"Given a branch of the phylogenetic tree and a branch multiplier,\n"
            "SigMA identifies alignment segments w.r.t. the given branch and\n"
            "assigns a p-value to each segment.\n\n"
            "In order to perform the downstream procedure of StatSigMAw, SigMA\n"
            "should be run on the combinations of ALL possible branches of the\n"
            "phylogenetic tree and THREE branch multipliers (0.01, 1 and 100).\n\n");
    fprintf(stderr, "We need to run SigMA on the combinations of ALL possible branches of the\n"
            "tree and THREE branch multipliers (0.01, 1 and 100). For example,\n"
            "the default tree provided in our code has 53 branches, and the \n"
            "branch index ranges between 0-52.  We will then run SigMA (53 * 3) = 159 times.\n");

    fprintf(stderr, "\nOUTPUT\n");
    fprintf(stderr, "The program generates the output to stdout, which includes information\n"
            "about identifed alignment segments, such as p-values.\n"
            "The stdout should be REDIRECTED to a file whose SUFFIX has to be\n"
            "\"<branch index>_<branch multiplier>.out\".\n");

    fprintf(stderr, "\nEXAMPLE\n");
    fprintf(stderr, "$ SigMA  myAlign.maf  6  0.01  >  myOutDir/myPrefix.6_0.01.out\n"
            "We show an example of running SigMA on branch #6 and multipler 0.01. \n"
            "Notice the output file has a suffix \"6_0.01.out\".\n");

    fprintf(stderr, "\nOPTIONS\n");
    fprintf(stderr, "   --phylogeny       input phylogeny in newick format default = %s\n", 
            globalOptions.PHYLOGENY);
    fprintf(stderr, "   --refSpecies      reference species default = %s\n", globalOptions.REF_SPECIES);
    fprintf(stderr, "   --maxBlockSize    maximum total size of contiguous alignment blocks allowed.\n" 
            "                     If the length of one alignment block is larger than value\n"
            "                     assigned, maxBlockSize will be assigned to the actual\n"
            "                     alignment block size."
            " default = %d\n", 
            globalOptions.MAX_BLOCK_SIZE);
    fprintf(stderr, "   --maxSegments     max number of high-scoring segments default = %d\n", 
            globalOptions.MAX_SEGMENTS);
    fprintf(stderr, "\nThese are used to estimate Karlin-Altschul parameters:\n\n");
    fprintf(stderr, "   --totalNumTuples  max number of tuples to be read per block (after \n"
            "                     randomization). default = %d\n", globalOptions.TOTAL_NUM_TUPLES);
    fprintf(stderr, "   --totalIterateParam  number of times that Karlin-Altschul parameter \n"
            "                     is computed. default = %d\n", globalOptions.TOTAL_ITERATE_PARAM);
    exit(1);
}

void loadDefaultParameters(void)
{
    strcpy(globalOptions.PHYLOGENY, "(((((((((((((hg:0.006690,chimp:0.007571):0.024272,(colobus_monkey:0.015404,(baboon:0.008258,macaque:0.028617):0.008519):0.022120):0.023960,(dusky_titi:0.025662,(owl_monkey:0.012151,marmoset:0.029549):0.008236):0.027158):0.066101,(mouse_lemur:0.059024,galago:0.121375):0.032386):0.017073,((rat:0.081728,mouse:0.077017):0.229273,rabbit:0.206767):0.023340):0.023026,(((cow:0.159182,dog:0.147731):0.004946,rfbat:0.138877):0.010150,(hedgehog:0.193396,shrew:0.261724):0.054246):0.024354):0.028505,armadillo:0.149862):0.015994,(elephant:0.104891,tenrec:0.259797):0.040371):0.218400,monodelphis:0.371073):0.065268,platypus:0.468116):0.123856,chicken:0.454691):0.123297,xenopus:0.782453):0.156067,((tetraodon:0.199381,fugu:0.239894):0.492961,zebrafish:0.782561):0.156067)");
    strcpy(globalOptions.REF_SPECIES, "hg");
}

void parseArgs(int argc, char **argv)
{
    static const char *optString = "vh?";
    static const struct option longOpts[] = {
        {"phylogeny", required_argument, NULL, 0},
        {"refSpecies", required_argument, NULL, 0},
        {"maxBlockSize", required_argument, NULL, 0},
        {"maxSegments", required_argument, NULL, 0},
        {"totalNumTuples", required_argument, NULL, 0},
        {"totalIterateParam", required_argument, NULL, 0},
        
        {"verbose", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
 
       {NULL, no_argument, NULL, 0}
    };
    int longIndex;
    int opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    while(opt != -1){
        switch(opt) {
        case 0:
            if(strcmp("phylogeny", longOpts[longIndex].name) == 0){
                if (!sscanf(optarg, "%4095s", globalOptions.PHYLOGENY)){
                    fprintf(stderr, "Unable to read --phylogeny\n");
                    exit(1);
                    break;
                }
            }else if(strcmp("refSpecies", longOpts[longIndex].name) == 0){
                if (!sscanf(optarg, "%31s", globalOptions.REF_SPECIES)){
                    fprintf(stderr, "Unable to read --refSpecies\n");
                    exit(1);
                    break;
                }
            }else if(strcmp("maxBlockSize", longOpts[longIndex].name) == 0){
                globalOptions.MAX_BLOCK_SIZE = atoi(optarg);
                break;
            }else if(strcmp("maxSegments", longOpts[longIndex].name) == 0){
                globalOptions.MAX_SEGMENTS = atoi(optarg);
                break;
            }else if(strcmp("totalNumTuples", longOpts[longIndex].name) == 0){
                globalOptions.TOTAL_NUM_TUPLES = atoi(optarg);
                break;
            }else if(strcmp("totalIterateParam", longOpts[longIndex].name) == 0){
                globalOptions.TOTAL_ITERATE_PARAM = atoi(optarg);
                break;
            }
        case 'v':
            globalOptions.verbose++;
            break;
        case '?':
        case 'h':
            usage();
            break;
        default:
            break;
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
        // globalOptions.inputFiles = argv + optind;
        // globalOptions.numInputFiles = argc - optind;
    }
}
