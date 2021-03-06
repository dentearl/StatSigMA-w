/*VI. Conditions of Use
 *
 * Xiaoyu Chen, Amol Prakash, and Martin Tompa give permission for
 * Dent Earl and his institution to use the StatSigmMA - w software
 * developed at the University of Washington for research purposes, on
 * the following conditions:
 *
 * 1. Dent Earl and his institutional colleagues may modify the
 *  StatSigmMA - w software and distribute the resulting modified
 *  software for research purposes, provided (a) that it is
 *  distributed together with these Conditions of Use, (b) that Martin
 *  Tompa receive a copy of the finalized modified software, and (c)
 *  that Xiaoyu Chen, Amol Prakash, and Martin Tompa are credited with
 *  the authorship of the software.  
 * 2. The StatSigmMA - w software
 *  will be used by you and/or your institution solely for
 *  noncommercial purposes, except with express permission from the
 *  authors.  
 * 3. Any risk associated with using the StatSigmMA - w
 *  software at your institution is with you and your institution.
 * 4. StatSigmMA - w will be cited in any publication(s) reporting on
 *  data obtained from it as: 
 * * Amol Prakash and Martin Tompa,
 *   "Measuring the Accuracy of Genome - Size Multiple
 *   Alignments". Genome Biology, vol. 8, issue 6, June 2007, R124.
 * * Xiaoyu Chen and Martin Tompa, "Comparative assessment of methods
 *   for aligning multiple genome sequences". Nature Biotechnology,
 *   vol. 28, no. 6, June 2010, 567 - 572.  
 */
#include <getopt.h>
#include "global_combine.h"
#include "tree_combine.h"

void loadDefaultParameters(void);
int parseArgs(int argc, char **argv);
void usage(void);

int main(int argc, char* argv[]) {
    loadDefaultParameters();
    parseArgs(argc, argv);
    cout << "#branch number\tspecies" << endl;
    read_newick_tree();
    if (g_options.printAll)
        return 0;
    int br;
    for (int i = 0; i < g_NUM_SPECIES; i++) {
        br = name2branch(species_name[i]);
        if (br > -1)
            cout << setw(2) << br << "\t" << species_name[i] << endl;
    }
    return 0;
}

void loadDefaultParameters(void) {
    strncpy(g_options.phylogeny, "(((((((((((((hg:0.006690,chimp:0.007571):0.024272,(colobus_monkey:0.015404,(baboon:0.008258,macaque:0.028617):0.008519):0.022120):0.023960,(dusky_titi:0.025662,(owl_monkey:0.012151,marmoset:0.029549):0.008236):0.027158):0.066101,(mouse_lemur:0.059024,galago:0.121375):0.032386):0.017073,((rat:0.081728,mouse:0.077017):0.229273,rabbit:0.206767):0.023340):0.023026,(((cow:0.159182,dog:0.147731):0.004946,rfbat:0.138877):0.010150,(hedgehog:0.193396,shrew:0.261724):0.054246):0.024354):0.028505,armadillo:0.149862):0.015994,(elephant:0.104891,tenrec:0.259797):0.040371):0.218400,monodelphis:0.371073):0.065268,platypus:0.468116):0.123856,chicken:0.454691):0.123297,xenopus:0.782453):0.156067,((tetraodon:0.199381,fugu:0.239894):0.492961,zebrafish:0.782561):0.156067)", kMaxLengthNewick);
}

void usage(void) {
    fprintf(stderr, "Usage: getSpeciesName [options]\n\n");
    fprintf(stderr, "\nDESCRIPTION\n");
    fprintf(stderr,"getSpeciesName is a utility for StatSigMAw,\n"
            "The program generates a mapping between branch indices and species names. \n"
            "For each of the branches incident on leaves, it tells which species that \n"
            "branch is incident on.\n"
            );
    fprintf(stderr, "\nOUTPUT\n");
    fprintf(stderr, "The program generates the output to stdout. The output lists \n"
            "the mapping between branch indices and species names. For example, \n"
            "given the default tree provided in our code, branch #6 correponds to \n"
            "platypus, which means the branch incident on platypus has index 6.\n"
            );
    fprintf(stderr, "\nEXAMPLE\n");
    fprintf(stderr, "$ getSpeciesName\n"
            "This will print out the branch number and names of all species in the tree.\n");

    fprintf(stderr, "\nOPTIONS\n");
    fprintf(stderr, "   --phylogeny       input phylogeny in newick format. default=%s\n", 
            g_options.phylogeny);
    fprintf(stderr, "   --printAll        print every branch in the tree.\n");
    exit(EXIT_FAILURE);
}

int parseArgs(int argc, char **argv) {
    static const char *optString = "vh?";
    static const struct option longOpts[] = {
        {"phylogeny", required_argument, NULL, 0},
        {"printAll", no_argument, NULL, 0},
        {"verbose", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
        {NULL, no_argument, NULL, 0}
    };
    int longIndex;
    int opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    while(opt != -1) {
        switch(opt) {
        case 0:
            if(strcmp("phylogeny", longOpts[longIndex].name) == 0) {
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
            } else if (strcmp("printAll", longOpts[longIndex].name) == 0) {
                g_options.printAll = 1;
                break;
            }
        case 'v':
            g_options.verbose++;
            break;
        case '?':
        case 'h':
            usage();
            break;
        default:
            break;
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    return optind;
}
