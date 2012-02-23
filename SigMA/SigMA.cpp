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
#include "common.h"
#include "global.h"
#include "tree.h"

Tree* root;
Tree* list_of_leaves[kMaxNumberBranches];

#include "pvalue.h"

// Helper Functions
Tree* read_tree(char newick_format[], int& sp_num, int& pos);
void generate_and_read_newick_tree();
void read_multi_pv(char* filename);
void init_background();
void verifyNewickBounds(int pos);
void loadDefaultParameters(void);
int parseArgs(int argc, char **argv);
void usage(void);
void reportOptions(void);

// Main program
int main(int argc, char* argv[]) {
    loadDefaultParameters();
    int n = parseArgs(argc, argv);
    if (argc - n != 4) {
        fprintf(stderr, "Error: wrong number of arguments.\n");
        usage();
    }
    for (int i = 0; i < argc; ++i)
        debug("argv[%d]: %s\n", i, argv[i]);
    strncpy(g_options.mafFile, argv[n], kMaxFilePath);
    g_options.branchIndex = atoi(argv[n + 1]);
    g_options.branchMultiplier = atof(argv[n + 2]);
    // Read the pvalue distribution of the various scores and number of segments
    strncpy(g_options.multipvFile, argv[n + 3], kMaxFilePath);
    reportOptions();
    char* pos = argv[0];
    while (strstr(pos, "SigMA") != NULL)
        pos = strstr(pos, "SigMA") + 1;
    pos--;
    char binary_path[1000];
    strncpy(binary_path, argv[0], pos - argv[0]);
    binary_path[int(pos - argv[0])] = '\0';
    
    read_multi_pv(g_options.multipvFile);
    debug("Read multipv %s\n", g_options.multipvFile);
    debug("Phylogeny: %s\n", g_options.phylogeny);
    generate_and_read_newick_tree();
    debug("Generated Tree: %s\n", root->tree_name);
    debug("Number of species = %d\n", g_NUM_SPECIES);
#ifdef PROTEIN
    sprintf(filename, "%s/pmb.dat", binary_path);
    read_aa_matrix(filename);
#ifdef DEBUG
    cout << "# Read AA matrix " << filename << endl;
#endif // DEBUG
#endif // PROTEIN

    init_background();  
    blocks_type blocks;
    double pvalue[2 * kMaxNumberBranches];
    Tree* branch_ptr[2 * kMaxNumberBranches];
    int branch = 0;
    compute_pvalue(root, blocks, pvalue, branch, branch_ptr, g_options.branchIndex, g_options.mafFile);
    double max_pvalue = -1.0;
    Tree* weakest_branch = NULL;
    for (int i = 0; i < branch; i++) {
        if (pvalue[i] > max_pvalue) {
            max_pvalue = pvalue[i];
            weakest_branch = branch_ptr[i];
        }
    }
    verbose("Number of blocks skipped due to missing reference: %d\n", g_NUM_BLOCKS_SKIPPED);
    cout << "p-value= " << max_pvalue << "\tbranch= " << weakest_branch->tree_name 
         << "\ttree = " << root->tree_name << endl;
    return EXIT_SUCCESS;
}

/*--------------------- Helper Functions ------------------ */
void verifyNewickBounds(int pos) {
    if (kMaxLengthNewick <= pos) {
        fprintf(stderr, "Error, poorly formed newick input from --phylogeny.\n");
        exit(EXIT_FAILURE);
    }
}
Tree* read_tree(char newick_format[], int& sp_num, int& pos) {
    verifyNewickBounds(pos);
    while (isblank(newick_format[pos])) {
        pos++;
        verifyNewickBounds(pos);
    }
    if (newick_format[pos] == '(') {
        pos++;
        Tree* left = read_tree(newick_format, sp_num, pos);
        while (newick_format[pos] != ':') {
            pos++;
            verifyNewickBounds(pos);
        }
        pos++;
        while (isblank(newick_format[pos])) {
            pos++;
            verifyNewickBounds(pos);
        }
        int start_pos = pos;
        while (newick_format[pos] != ',') {
            pos++;
            verifyNewickBounds(pos);
        }
        char branch_length[100];
        strncpy(branch_length, newick_format + start_pos, pos-start_pos);
        branch_length[pos - start_pos] = '\0';
        float left_branch = atof(branch_length);
        if (left_branch < 0.001) left_branch = 0.001;
        left_branch *= g_options.branchMultiplier;
        if (left_branch > kMaxBranchLength) 
            left_branch = kMaxBranchLength;
        pos++;
        Tree* right = read_tree(newick_format, sp_num, pos);
        while (newick_format[pos] != ':') {
            pos++;
            verifyNewickBounds(pos);
        }
        pos++;
        verifyNewickBounds(pos);
        while (isblank(newick_format[pos])) {
            pos++;
            verifyNewickBounds(pos);
        }
        start_pos = pos;
        while (newick_format[pos] != ')') {
            pos++;
            verifyNewickBounds(pos);
        }
        strncpy(branch_length, newick_format + start_pos, pos - start_pos);
        branch_length[pos - start_pos] = '\0';
        float right_branch = atof(branch_length);
        if (right_branch < 0.001) right_branch = 0.001;
        right_branch *= g_options.branchMultiplier;
        if (right_branch > kMaxBranchLength) 
            right_branch = kMaxBranchLength; 
        pos++;
        verifyNewickBounds(pos);
        return (new Tree(left, right, left_branch, right_branch));
    } else {
        int start_pos = pos;
        while (newick_format[pos] != ':') {
            pos++;
            verifyNewickBounds(pos);
        }
        strncpy(g_SPECIES_NAMES[sp_num], newick_format + start_pos, pos - start_pos);
        g_SPECIES_NAMES[sp_num][pos - start_pos] = '\0';    
        Tree* node = new Tree(sp_num, g_SPECIES_NAMES[sp_num]);
        list_of_leaves[sp_num] = node;
        sp_num++;
        return node;
    }
}
void generate_and_read_newick_tree() {
    int pos = 0;
    g_NUM_SPECIES = 0;
    root = read_tree(g_options.phylogeny, g_NUM_SPECIES, pos);
}
void read_multi_pv(char* filename) {
    g_MULTI_PV = new double*[g_options.maxSegments + 1];
    for (int i = 0; i < g_options.maxSegments + 1; i++)
        g_MULTI_PV[i] = new double[int(kMaxTotalScore - kMinTotalScore + 1)];
    for (int i = 0; i < (g_options.maxSegments + 1); i++)
        for (int j = 0; j < (kMaxTotalScore - kMinTotalScore + 1); j++)
            g_MULTI_PV[i][j]=-1;
    for (int j = kMinTotalScore; j <= kMaxTotalScore; j++)
        g_MULTI_PV[1][j - kMinTotalScore] = 1 - exp(-1 * exp(-1 * j));
    ifstream ifs(filename);
    int a, b;
    double c;
    if (!ifs) {
        cerr << "Unable to open pvalue file `" << filename << "'" << endl;
        exit(EXIT_FAILURE);
    }
    while (ifs) {
        ifs >> a >> b >> c;
        if ((a <= g_options.maxSegments) && (b <= kMaxTotalScore) && (b >= kMinTotalScore)) {
            g_MULTI_PV[a][b - kMinTotalScore] = c;
        }
    }
    ifs.close();
}
void init_background() {
    g_back[0] = 0;
#ifdef PROTEIN
    for (int i = 1; i < kAlphabetSize; i++)
        g_back[i] = aa_back[i - 1] * (1 - kGapBackground);
#else
    for (int i = 1; i < kAlphabetSize; i++)
        g_back[i] = (1. - kGapBackground - g_back[0]) / float(kAlphabetSize - 1.);
#endif // PROTEIN
    g_back[kGap] = kGapBackground;
}
void usage(void) {
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
            "branch index ranges between 0-52.  We will then run \n"
            "StatSigMAw (53 * 3) = 159 times.\n");
    
    fprintf(stderr, "\nASSUMPTIONS\n");
    fprintf(stderr, "  * The reference species must be present in every block.\n"
            "  * There may be no duplicate species in any block. A species may be \n"
            "    present only once per block.\n"
            "  * The blocks must be ordered such that there is a monotonic increase \n"
            "    in position relative to the reference sequence.\n");

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
    fprintf(stderr, "   --phylogeny       input phylogeny in newick format default=%s\n", 
            g_options.phylogeny);
    fprintf(stderr, "   --refSpecies      reference species. Must be present in an alignment\n"
            "                     block in order for that block to be analyzed. default=%s\n", 
            g_options.refSpecies);
    fprintf(stderr, "   --maxBlockSize    maximum total size of contiguous alignment blocks allowed.\n" 
            "                     If the length of one alignment block is larger than value\n"
            "                     assigned, maxBlockSize will be assigned to the actual\n"
            "                     alignment block size."
            " default=%d\n", 
            g_options.maxBlockSize);
    fprintf(stderr, "   --maxSegments     max number of high-scoring segments default=%d\n", 
            g_options.maxSegments);
    fprintf(stderr, "\nThese are used to estimate Karlin-Altschul parameters:\n\n");
    fprintf(stderr, "   --totalNumberTuples  max number of tuples to be read per block (after \n"
            "                     randomization). default=%d\n", g_options.totalNumberTuples);
    fprintf(stderr, "   --totalIterateParam  number of times that Karlin-Altschul parameter \n"
            "                     is computed. default=%d\n", g_options.totalIterateParam);
    exit(EXIT_FAILURE);
}
void loadDefaultParameters(void) {
    strcpy(g_options.phylogeny, "(((((((((((((hg:0.006690,chimp:0.007571):0.024272,(colobus_monkey:0.015404,(baboon:0.008258,macaque:0.028617):0.008519):0.022120):0.023960,(dusky_titi:0.025662,(owl_monkey:0.012151,marmoset:0.029549):0.008236):0.027158):0.066101,(mouse_lemur:0.059024,galago:0.121375):0.032386):0.017073,((rat:0.081728,mouse:0.077017):0.229273,rabbit:0.206767):0.023340):0.023026,(((cow:0.159182,dog:0.147731):0.004946,rfbat:0.138877):0.010150,(hedgehog:0.193396,shrew:0.261724):0.054246):0.024354):0.028505,armadillo:0.149862):0.015994,(elephant:0.104891,tenrec:0.259797):0.040371):0.218400,monodelphis:0.371073):0.065268,platypus:0.468116):0.123856,chicken:0.454691):0.123297,xenopus:0.782453):0.156067,((tetraodon:0.199381,fugu:0.239894):0.492961,zebrafish:0.782561):0.156067)");
    strcpy(g_options.refSpecies, "hg");
}
int parseArgs(int argc, char **argv) {
    static const char *optString = "dvh?";
    static const struct option longOpts[] = {
        {"phylogeny", required_argument, NULL, 0},
        {"refSpecies", required_argument, NULL, 0},
        {"maxBlockSize", required_argument, NULL, 0},
        {"maxSegments", required_argument, NULL, 0},
        {"totalNumberTuples", required_argument, NULL, 0},
        {"totalIterateParam", required_argument, NULL, 0},
        {"rseed", required_argument, NULL, 0},
        {"debug", no_argument, NULL, 'd'},
        {"verbose", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
        {NULL, no_argument, NULL, 0}
    };
    int longIndex, setSeed = 0;
    int opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    while (opt != -1) {
        switch(opt) {
        case 'v':
            g_options.verbose++;
            break;
        case 'd':
            g_options.debug = true;
            break;
        case 0:
            if (strcmp("phylogeny", longOpts[longIndex].name) == 0) {
                if (!sscanf(optarg, "%4095s", g_options.phylogeny)) {
                    fprintf(stderr, "Unable to read --phylogeny\n");
                    exit(EXIT_FAILURE);
                    break;
                }
                if (strlen(g_options.phylogeny) == (unsigned) kMaxLengthNewick - 1) {
                    fprintf(stderr, "Unable to read --phylogeny, too large!\n");
                    exit(EXIT_FAILURE);
                    break;
                }
                break;
            } else if (strcmp("refSpecies", longOpts[longIndex].name) == 0) {
                if (!sscanf(optarg, "%31s", g_options.refSpecies)) {
                    fprintf(stderr, "Unable to read --refSpecies\n");
                    exit(EXIT_FAILURE);
                    break;
                }
                if (strlen(g_options.refSpecies) == (unsigned) kMaxLengthRefSpecies - 1) {
                    fprintf(stderr, "Unable to read --refSpecies, too large!\n");
                    exit(EXIT_FAILURE);
                    break;
                }
                break;
            } else if (strcmp("maxBlockSize", longOpts[longIndex].name) == 0) {
                g_options.maxBlockSize = atoi(optarg);
                break;
            } else if (strcmp("maxSegments", longOpts[longIndex].name) == 0) {
                g_options.maxSegments = atoi(optarg);
                break;
            } else if (strcmp("totalNumberTuples", longOpts[longIndex].name) == 0) {
                g_options.totalNumberTuples = atoi(optarg);
                break;
            } else if (strcmp("totalIterateParam", longOpts[longIndex].name) == 0) {
                g_options.totalIterateParam = atoi(optarg);
                break;
            } else if (strcmp("rseed", longOpts[longIndex].name) == 0) {
                int tmp = atoi(optarg);
                if (tmp < 0)
                    tmp = -tmp;
                g_options.rseed = (unsigned) tmp;
                break;
            }
        case '?':
        case 'h':
            usage();
            break;
        default:
            break;
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    if (!setSeed) {
        g_options.rseed = ((time(NULL) & 0xFFFF) | (getpid() << 16)); // likely unique
    }
    srand(g_options.rseed);
    return optind;
}
void message(char const *type, char const *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "%s: ", type);
    vfprintf(stderr, fmt, args);
    va_end(args);
}
void verbose(char const *fmt, ...) {
    char str[kMaxLengthMessage];
    va_list args;
    va_start(args, fmt);
    if (g_options.verbose) {
        int n = vsprintf(str, fmt, args);
        if (n >= kMaxLengthMessage) {
            fprintf(stderr, "Error, failure in debug(), (n = %d) > "
                    "(kMaxLengthMessage = %d)\n", n, kMaxLengthMessage);
        }
        message("Verbose", str, args);
    }
    va_end(args);
}
void debug(char const *fmt, ...) {
    char str[kMaxLengthMessage];
    va_list args;
    va_start(args, fmt);
    if (g_options.debug) {
        int n = vsprintf(str, fmt, args);
        if (n >= kMaxLengthMessage) {
            fprintf(stderr, "Error, failure in debug(), (n = %d) > "
                    "(kMaxLengthMessage = %d)\n", n, kMaxLengthMessage);
            exit(EXIT_FAILURE);
        }
        message("Debug", str, args);
    }
    va_end(args);
}
void reportOptions(void) {
    debug("Options:\n");
    debug("  mafFile=%s\n", g_options.mafFile);
    debug("  branchIndex=%d\n", g_options.branchIndex);
    debug("  branchMultiplier=%f\n", g_options.branchMultiplier);
    debug("  multipvFile=%s\n", g_options.multipvFile);
    debug("  --maxBlockSize=%d\n", g_options.maxBlockSize);
    debug("  --maxSegments=%d\n", g_options.maxSegments);
    debug("  --totalNumberTuples=%d\n", g_options.totalNumberTuples);
    debug("  --totalIterateParam=%d\n", g_options.totalIterateParam);
    debug("  --phylogeny=%s\n", g_options.phylogeny);
    debug("  --refSpecies=%s\n", g_options.refSpecies);
    debug("  --rseed=%u\n", g_options.rseed);
    debug("  --verbose=%d\n", g_options.verbose);
    debug("  --debug=%d\n", g_options.debug);
}
