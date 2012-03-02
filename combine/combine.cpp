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
#include "common.h"

void loadDefaultParameters(void);
int parseArgs(int argc, char **argv);
void usage(void);

float get_PV(int i) {
    return g_CUR_PV[i];
}
void set_PV(int i, float pv) { 
    g_CUR_PV[i] = pv;
}
bool get_G_IsPresent(int i, int sp) { 
    return (g_isPresent[sp][i]); 
}
bool is_NoneG_IsPresent(int i) { 
    bool res = true;
    for (int s = 0; s < branch_num; s++)
        if (g_isPresent[s][i]) {
            res = false;
            break;
        }
    return res; 
}
void set_G_IsPresent(int i, int sp, bool cur_g_isPresent) { 
    g_isPresent[sp][i] = cur_g_isPresent;
}
bool get_Gap(int i, int sp) { 
    return (g_isGaps[sp][i]); 
}
void set_Gap(int i, int sp, bool cur_g_isGaps) {
    g_isGaps[sp][i] = cur_g_isGaps;
}
// Initialize g_isWorstSet for all species at position i.
void init_WorstSet(int i) {
    for (int s = 0; s < (branch_num + 1); s++)
        g_isWorstSet[s][i] = false;
}
void add_WorstBranch(int i, int sp) {
    g_isWorstSet[sp][i] = true;
    if (sp < branch_num)
        for (int s = 0; s < branch_num; s++)
            // check if sp is on the path from s to reference
            if ((get_G_IsPresent(i, s)) && (is_OnPathToReference(s, sp)))
                g_isWorstSet[s][i] = true;
}
bool isin_WorstSet(int i, int sp) {
    return (g_isWorstSet[sp][i]);
}
void read_Present() {
    // Initialize
    for (int i = 0; i < g_ALIGN_LEN; i++)
        for (int s = 0; s < branch_num; s++)
            set_G_IsPresent(i, s, false);
    int cur_pos = 0;
    int cur_len = 0;
    bool *cur_g_isPresent = new bool[branch_num];
    for (int s = 0; s < branch_num; s++)
        cur_g_isPresent[s] = false; 
    ifstream ifs(g_options.mafFile);
    if (!ifs.is_open()) {
        cerr << "Error, unable to open maf file " << g_options.mafFile << endl;
        exit(EXIT_FAILURE);
    }
    debug("Opening maf %s\n", g_options.mafFile);
    char *buffer = new char[g_MAX_BLOCK_SIZE];
    // boolb firstBlock = true;
    // int start_pos = 0;
    while (ifs) {
        buffer[0]='\0';
        ifs.getline(buffer, g_MAX_BLOCK_SIZE);
        if ((strncmp(buffer, "a score=", 8) == 0) || (!ifs)) {
            insert_inner_branches(cur_g_isPresent);
            for (int i = 0; i < cur_len; i++)
                for (int s = 0; s < branch_num; s++)
                    set_G_IsPresent(cur_pos + i, s, cur_g_isPresent[s]);
            // Reset
            cur_pos += cur_len;
            cur_len = 0;
            for (int s = 0; s < branch_num; s++)
                cur_g_isPresent[s] = false;
        } else {
            char* a;
            if (buffer[0] == 's') {
                // Checking for reference
                bool is_reference = false;
                if (strncmp(buffer + 2, g_options.refSpecies, strlen(g_options.refSpecies)) == 0)
                    is_reference = true;
                a = strtok(buffer, " "); // s
                a = strtok(NULL, " "); // name
                int species_line = name2branch(a);
                a = strtok(NULL, " "); // start
                a = strtok(NULL, " "); // length
                a = strtok(NULL, " "); // strand
                a = strtok(NULL, " "); // source length
                a = strtok(NULL, " "); // alignment field
                if (is_reference)
                    cur_len = strlen(a);
                if (a[0] != 'N')
                    cur_g_isPresent[species_line] = true;
            }
        }
    }
    delete [] cur_g_isPresent;
    delete [] buffer;
    ifs.close();
}
void read_Gaps() {
    for (int i = 0; i < g_ALIGN_LEN; i++)
        for (int s = 0; s < branch_num; s++)
            set_Gap(i, s, false);
  
    int cur_pos = 0;
    int cur_len = 0;
    char *buffer = new char[g_MAX_BLOCK_SIZE];
    bool **cur_g_isGaps;
    cur_g_isGaps = new bool*[branch_num];
    for (int s = 0; s < branch_num; s++)
        cur_g_isGaps[s] = new bool[g_MAX_BLOCK_SIZE];
    for (int s = 0; s < branch_num; s++)
        for (int i = 0; i < g_MAX_BLOCK_SIZE; i++)
            cur_g_isGaps[s][i] = true;
    ifstream ifs(g_options.mafFile);
    while (ifs) {
        buffer[0]='\0';
        ifs.getline(buffer, g_MAX_BLOCK_SIZE);
        if ((strncmp(buffer, "a score=", 8) == 0) || (!ifs)) {      
            bool *tmp_g_isGaps = new bool[branch_num];
            for (int i = 0; i < cur_len; i++) {
                for (int s = 0; s < branch_num; s++)
                    tmp_g_isGaps[s] = cur_g_isGaps[s][i];
                check_gapped_branches(tmp_g_isGaps);
                for (int s = 0; s < branch_num; s++)
                    set_Gap(cur_pos + i, s, tmp_g_isGaps[s]);
            }
            delete [] tmp_g_isGaps;
            // Reset
            cur_pos += cur_len;
            for (int s = 0; s < branch_num; s++)
                for (int i = 0; i < g_MAX_BLOCK_SIZE; i++)
                    cur_g_isGaps[s][i]=true;
            cur_len = 0;
        } else {
            char* a;
            if (buffer[0] == 's') {
                // Checking for reference species
                bool is_reference = false;
                if (strncmp(buffer + 2, g_options.refSpecies, strlen(g_options.refSpecies)) == 0)
                    is_reference = true;
                a = strtok(buffer, " ");
                a = strtok(NULL, " ");
                int species_line = name2branch(a);
                a = strtok(NULL, " ");
                a = strtok(NULL, " ");
                a = strtok(NULL, " ");
                a = strtok(NULL, " ");
                a = strtok(NULL, " ");
                if (is_reference)
                    cur_len = strlen(a);
                // a includes g_isGaps in the reference species here.
                for (unsigned i = 0; i < strlen(a); i++)
                    if ((a[i]!='-') && (a[i] != 'N'))
                        cur_g_isGaps[species_line][i] = false;
            }
        }
    }
    ifs.close();
    for (int s = 0; s < branch_num; s++)
        delete [] cur_g_isGaps[s];
    delete [] cur_g_isGaps;
    delete [] buffer;
}
bool read_segment(ifstream &ifs, int& start, int& end, int& left_nearest, int& right_nearest, double& pv) {
    char text[10000];
    int segment, block, block_size;
    double score;
    if (!ifs) return false;
    char errorMsg[50] = "# number of sets of parameter computed= 0";
    text[0] = text[2] = '\0';
    ifs.getline(text, 10000);
    while ((ifs) && (text[2] != 'S')) {
        // if no parameter estimated, abort this file.
        if (strcmp(text, errorMsg) == 0) {
            return false;
        }
        text[0]=text[2] = '\0';
        ifs.getline(text, 10000);
    }
    if (text[2]=='S') {
        char* a = strtok(text + 2, " \t");
        a = strtok(NULL, " \t");
        segment = atoi(a);
        a = strtok(NULL, " \t");
        block = atoi(a);
        a = strtok(NULL, " \t");
        block_size = atoi(a);
        a = strtok(NULL, " \t");
        score = atof(a);
        a = strtok(NULL, " \t");
        a = strtok(NULL, " \t");
        a = strtok(NULL, " \t");
        pv = atof(a);
        a = strtok(NULL, " \t");
        left_nearest = atoi(a);
        a = strtok(NULL, " \t");
        right_nearest = atoi(a);
        a = strtok(NULL, " \t");
        a = strtok(NULL, " \t");
        a = strtok(NULL, " \t");
        a = strtok(NULL, " \t");
        start = atoi(a);
        a = strtok(NULL, " \t");
        end = atoi(a);
        if (end > g_ALIGN_LEN) end = g_ALIGN_LEN;
        return true;
    }
    return false;
}
// Identifying good regions  (removed gap_check, because now g_isGaps are counted after outputting good regions)
void output_good_regions() {
    int start = -1;
    for (int i = 0; i < g_ALIGN_LEN; i++) {
        if (get_PV(i) < g_options.pvThreshGood) {
            if (start == -1)
                start = i;
        }
        else if (start != -1) {
            int currLen = i - start;
            if (currLen >= g_options.minSegSize) {
                // the first and the last reference coordinates in current region                
                int c = start;
                while ((refSpeciesCoord[c] == -1) && (c < g_ALIGN_LEN - 1))
                    c++;
                int firstHgCoor = refSpeciesCoord[c];
	
                c = start + currLen - 1;
                while ((refSpeciesCoord[c] == -1) && (c > 0))
                    c--;
                int lastHgCoor = refSpeciesCoord[c];
	
                cout << "my_good_region " << " " << start << " " << currLen <<  " " 
                     << g_refSpeciesChr << " " << firstHgCoor + 1 << " " << lastHgCoor + 1 << endl;
            }
            start = -1;
        }
    }
}
void read_output_pv(const char* branch_multiplier) {
    for (int sp = 0; sp < branch_num; sp++) {
        int start, end, left_nearest, right_nearest;
        double pv;
        char filename[1000];
        sprintf(filename, "%s/%s.%d_%s.out", 
                g_options.SigMAwOutDir, 
                g_options.SigMAwOutPrefix, 
                sp, branch_multiplier);
        //cout<<filename<<endl;
        ifstream ifs(filename);
        while (read_segment(ifs, start, end, left_nearest, right_nearest, pv))
            for (int i = start; i<=end; i++)
                // get the max PValue for all branches
                if ((get_G_IsPresent(i, sp)) && (get_PV(i) < pv))
                    set_PV(i, pv);
        ifs.close();
    }
}
void read_output_worstbranch(const char* branch_multiplier) {
    // get the tree ID for current tree
    char treeID;
    if (strcmp(branch_multiplier, "100") == 0) 
        treeID = '2';
    if (strcmp(branch_multiplier, "1") == 0) 
        treeID = '1';
    if (strcmp(branch_multiplier, "0.01") == 0) 
        treeID = '0';
    for (int sp = 0; sp < branch_num; sp++) {
        int start, end, left_nearest, right_nearest;
        double pv;
        char filename[1000];
        sprintf(filename, "%s/%s.%d_%s.out", 
                g_options.SigMAwOutDir, 
                g_options.SigMAwOutPrefix, 
                sp, branch_multiplier);
        ifstream ifs(filename);
        while (read_segment(ifs, start, end, left_nearest, right_nearest, pv))
            for (int i = start; i<=end; i++) {
                // when the branch is g_isPresent at the position, 
                // the pvalue of the branch equals to min(max pvalue), 
                // and current tree is the reponsible tree for such pvalue.
                if (get_G_IsPresent(i, sp) && (get_PV(i) == pv) && (g_TREE_PV[i] == treeID))
                    add_WorstBranch(i, sp);
            }
        ifs.close();
    }
}
void read_pv_dist(char* branch_multiplier, int sp) {
    double pv_dist[10];
    for (int i = 0; i < 10; i++)
        pv_dist[i]=0;
    int total = 0;
    int start, end, left_nearest, right_nearest;
    double pv;
    char filename[1000];
    sprintf(filename, "%s/%s.%d_%s.out", 
            g_options.SigMAwOutDir, 
            g_options.SigMAwOutPrefix, 
            sp, branch_multiplier);
    ifstream ifs(filename);
    while (read_segment(ifs, start, end, left_nearest, right_nearest, pv))
        if (end - start >= 0) {
            int total_start_end = 0;
            for (int j = start; j < end; j++)
                if (get_G_IsPresent(j, sp))
                    total_start_end++;
            if (pv < 1e-9) 
                pv = 1e-9;
            pv = -1 * log10(pv);
            pv_dist[int(pv)] += total_start_end;
            total += total_start_end;
        }
    ifs.close();
    for (int i = 0; i < 10; i++)
        cout << "my_pv " << sp << " " << i << " " << pv_dist[i]/float(total) 
             << " " << pv_dist[i] << endl; 
}
void scan_alignFile() {
    g_ALIGN_LEN = 0;
    int cur_len = 0;
    char *buffer = new char[g_MAX_BLOCK_SIZE];
    bool firstBlock = true;
    int currStart = 0;
    int currNtNum = 0;
    ifstream ifs(g_options.mafFile);
    while (ifs) {
        buffer[0]='\0';
        ifs.getline(buffer, g_MAX_BLOCK_SIZE);
        if ((strncmp(buffer, "a score=", 8)==0) || (!ifs)) {
            g_ALIGN_LEN += cur_len;
            cur_len = 0;
        }
        else {
            char* a;
            if (buffer[0]=='s') {
                // Checking for reference
                bool is_reference = false;
                if (strncmp(buffer + 2, g_options.refSpecies, strlen(g_options.refSpecies)) == 0)
                    is_reference = true;
                if (is_reference) {
                    a = strtok(buffer, " ");
                    a = strtok(NULL, " ");
                    a = strtok(NULL, " ");
                    currStart = atoi(a);
                    if (firstBlock) {
                        g_CHR_START = atoi(a);
                        firstBlock = false;
                    }
                    a = strtok(NULL, " ");
                    currNtNum = atoi(a);
                    a = strtok(NULL, " ");
                    a = strtok(NULL, " ");
                    a = strtok(NULL, " ");
                    cur_len = strlen(a);
                }
            }
        }
    }
    ifs.close();
    delete [] buffer;
}
void map_referenceCoordinates() {
    int cumStart = 0;
    int cur_len = 0;
    char *buffer = new char[g_MAX_BLOCK_SIZE];
    bool firstBlock = true;
    char *orient;
    int currHgPos = 0;
    int chrLen = 0;
    ifstream ifs(g_options.mafFile);
    while (ifs) {
        buffer[0] = '\0';
        ifs.getline(buffer, g_MAX_BLOCK_SIZE);
        if ((strncmp(buffer, "a score=", 8) == 0) || (!ifs)) {
            cumStart += cur_len;
            cur_len = 0;
        }
        else {
            char* a;
            if (buffer[0]=='s') {
                // Checking for reference
                bool is_reference = false;
                if (strncmp(buffer + 2, g_options.refSpecies, strlen(g_options.refSpecies)) == 0)
                    is_reference = true;
                if (is_reference) {
                    a = strtok(buffer, " ");
                    a = strtok(NULL, " ");
                    // get the chr info
                    if (firstBlock) {
                        char* temp = strstr(a, "chr");
                        strcpy(g_refSpeciesChr, temp);
                        firstBlock = false;
                    }
                    a = strtok(NULL, " ");
                    currHgPos = atoi(a);
                    a = strtok(NULL, " ");
                    orient = strtok(NULL, " ");
	  
                    a = strtok(NULL, " ");
                    chrLen = atoi(a);
                    a = strtok(NULL, " ");
                    cur_len = strlen(a);
                    if (strcmp(orient, "-") == 0) {
                        cout<<"Warning: there is negative strand in reference!"<<endl;
                        currHgPos = chrLen - currHgPos - 1;
                    }
                    for (int i = 0; i < cur_len; i++) {
                        // 'N' is actually genomic positions
                        if (a[i]=='-')
                            refSpeciesCoord[cumStart + i] = -1;
                        else {
                            refSpeciesCoord[cumStart + i] = currHgPos;
                            currHgPos++;
                        }
                    }
                }
            }
        }
    }
    ifs.close();
    delete [] buffer;
}
void usage(void) {
    fprintf(stderr, "Usage: combine <maf file> <SigMA output dir> <SigMA output prefix> [options]\n\n");
    fprintf(stderr, "  <maf file>: the multiple sequence alignment in maf format.\n");
    fprintf(stderr, "  <SigMA output dir>: the directory where the SigMA output files are located.\n");
    fprintf(stderr, "  <SigMA output prefix>: the prefix of SigMA output file names.\n");
    fprintf(stderr, "\nDESCRIPTION\n");
    fprintf(stderr, "Given the directory where the SigMA output files are located \n"
            "and the prefix of SigMA output files, the alignment segments \n"
            "identified by SigMA for all branches in the phylogenetic tree \n"
            "and for all three branch multipliers are assessed. The program \n"
            "finally identifies suspiciously aligned regions with respect to \n"
            "each branch.\n\n"
            "NOTICE: This program assumes that SigMA has been run on the \n"
            "combinations of ALL possible branches of the phylogenetic tree \n"
            "and THREE branch multipliers (0.01, 1 and 100).\n"
            );
    fprintf(stderr, "\nOUTPUT\n"
            "The program generates the output to stdout. The output consists \n"
            "of three parts:\n"
            "  a) Regions that are well aligned with respect to ALL branches.\n"
            "     The format of this part of output is:\n"
            "    - my_good_region (simply a flag) \n"
            "    - start position of this region (in the alignment coordinate,\n"
            "      i.e. which column in the alignment is this position)\n"
            "    - length of this region\n"
            "    - chromosome\n"
            "    - start position\n"
            "    - end position\n\n"
            "The last three items show the genomic coordinates of the reference \n"
            "species for this region.\n\n"
            "  b) Regions that are suspiciously aligned w.r.t. a branch. The \n"
            "     format of this part of output is:\n"
            "    - my_bad_region (simply a flag) \n"
            "    - start position of this region (in the alignment coordinate, \n"
            "      i.e. which column in the alignment is this position)\n"
            "    - branch index (the region is suspiciously aligned with respect \n"
            "      to this particular branch)\n"
            "    - length of this region\n"
            "    - chromosome\n"
            "    - start position\n"
            "    - end position\n\n"
            "Again, the last three items show the genomic coordinates of the \n"
            "reference species for this region.\n\n"
            "Note: For regions that are suspiciously aligned with respect to a \n"
            "branch incident on a leaf (i.e., a species), you can read them as \n"
            "regions where that species is suspiciously aligned to the other \n"
            "species.\n\n"
            "  c) Statistics for suspiciously-aligned regions. The format of \n"
            "     this part of output is:\n"
            "    - my_count (simply a flag)\n"
            "    - branch index\n"
            "    - number of aligned bases w.r.t. the branch\n"
            "    - number of suspiciously aligned bases w.r.t. the branch\n"
            "    - %% of suspiciously aligned bases w.r.t. the branch\n"
            "    - number of suspiciously aligned regions w.r.t. the branch\n\n"
            "NOTICE: Please refer to parameters.txt for parameter settings of \n"
            "well- or suspiciously-aligned regions.\n"
            );
    fprintf(stderr, "\nEXAMPLE\n"
            "   combine  myAlign.maf  myOutDir  myPrefix  >  myStatSigMAw.res\n"
            );
    exit(EXIT_FAILURE);
}
void loadDefaultParameters(void) {
    strcpy(g_options.phylogeny, "(((((((((((((hg:0.006690,chimp:0.007571):0.024272,(colobus_monkey:0.015404,(baboon:0.008258,macaque:0.028617):0.008519):0.022120):0.023960,(dusky_titi:0.025662,(owl_monkey:0.012151,marmoset:0.029549):0.008236):0.027158):0.066101,(mouse_lemur:0.059024,galago:0.121375):0.032386):0.017073,((rat:0.081728,mouse:0.077017):0.229273,rabbit:0.206767):0.023340):0.023026,(((cow:0.159182,dog:0.147731):0.004946,rfbat:0.138877):0.010150,(hedgehog:0.193396,shrew:0.261724):0.054246):0.024354):0.028505,armadillo:0.149862):0.015994,(elephant:0.104891,tenrec:0.259797):0.040371):0.218400,monodelphis:0.371073):0.065268,platypus:0.468116):0.123856,chicken:0.454691):0.123297,xenopus:0.782453):0.156067,((tetraodon:0.199381,fugu:0.239894):0.492961,zebrafish:0.782561):0.156067)");
    strcpy(g_options.refSpecies, "hg");
}
int parseArgs(int argc, char **argv) {
    static const char *optString = "vh?";
    static const struct option longOpts[] = {
        {"phylogeny", required_argument, NULL, 0},
        {"refSpecies", required_argument, NULL, 0},
        {"minSegSize", required_argument, NULL, 0},
        {"pThreshBad", required_argument, NULL, 0},
        {"pThreshGood", required_argument, NULL, 0},
        
        {"verbose", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
 
       {NULL, no_argument, NULL, 0}
    };
    int longIndex;
    int opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    while(opt != -1) {
        switch(opt) {
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
            } else if (strcmp("minSegSize", longOpts[longIndex].name) == 0) {
                g_options.minSegSize = atoi(optarg);
                break;
            } else if (strcmp("pThreshBad", longOpts[longIndex].name) == 0) {
                g_options.pvThreshBad = atof(optarg);
                break;
            } else if (strcmp("pThreshGood", longOpts[longIndex].name) == 0) {
                g_options.pvThreshGood = atof(optarg);
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
        // g_options.inputFiles = argv + optind;
        // g_options.numInputFiles = argc - optind;
    }
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
    char str[kMaxMessageLength];
    va_list args;
    va_start(args, fmt);
    if (g_options.verbose) {
        int n = vsprintf(str, fmt, args);
        if (n >= kMaxMessageLength) {
            fprintf(stderr, "Error, failure in debug(), (n = %d) > "
                    "(kMaxMessageLength = %d)\n", n, kMaxMessageLength);
        }
        message("Verbose", str, args);
    }
    va_end(args);
}
void debug(char const *fmt, ...) {
    char str[kMaxMessageLength];
    va_list args;
    va_start(args, fmt);
    if (g_options.debug) {
        int n = vsprintf(str, fmt, args);
        if (n >= kMaxMessageLength) {
            fprintf(stderr, "Error, failure in debug(), (n = %d) > "
                    "(kMaxMessageLength = %d)\n", n, kMaxMessageLength);
            exit(EXIT_FAILURE);
        }
        message("Debug", str, args);
    }
    va_end(args);
}
int main(int argc, char* argv[]) {
    loadDefaultParameters();
    int n = parseArgs(argc, argv);
    if (argc - n != 4) {
        fprintf(stderr, "Error: wrong number of arguments.\n");
        usage();
    }
    
    strcpy(g_options.mafFile, argv[n]);
    strcpy(g_options.SigMAwOutDir, argv[n + 1]);
    strcpy(g_options.SigMAwOutPrefix, argv[n + 2]);
  
    scan_alignFile();
    cout << "alignmnt length: "<< g_ALIGN_LEN << endl;
    g_CUR_PV = new float[g_ALIGN_LEN];
    g_TEMP_PV = new float[g_ALIGN_LEN];
    // record the tree ID that has min(max PValue)
    g_TREE_PV = new char[g_ALIGN_LEN];
    read_newick_tree();
    assert(branch_num * g_ALIGN_LEN < 1e9);
    cout << "read tree" << endl;  
    // memory allocate
    g_isPresent = new bool*[branch_num];
    for (int s = 0; s < branch_num; s++) 
        g_isPresent[s] = new bool[g_ALIGN_LEN];
    // read in g_isPresent[] for each branch - position, 
    // indicating if the branch is "g_isPresent" at the  position
    read_Present();
    cout << "read Present" << endl;
    // read in curr_pv[], the min(max pvalue over branches) over trees.
    // Firstly, identify the max pvalue among different branches. 
    // Then, identify the min pvalue among the three trees.
    // Using g_TEMP_PV as the structure to store temporary pvalue (as the 3 trees are read)
  
    // initialize
    for (int i = 0; i < g_ALIGN_LEN; i++) {
        set_PV(i, 0);
        g_TEMP_PV[i] = 1;
        g_TREE_PV[i] = 'n';
    }
    // for tree "100", ID is 2.
    read_output_pv("100");
    for (int i = 0; i < g_ALIGN_LEN; i++) {
        if ((g_TEMP_PV[i] > g_CUR_PV[i]) && (!is_NoneG_IsPresent(i))) {
            g_TEMP_PV[i] = g_CUR_PV[i];
            g_TREE_PV[i] = '2';
        }
        set_PV(i, 0);
    }
    // for tree "0.01", ID is 0
    read_output_pv("0.01");
    for (int i = 0; i < g_ALIGN_LEN; i++) {
        if ((g_TEMP_PV[i] > g_CUR_PV[i]) && (!is_NoneG_IsPresent(i))) {
            g_TEMP_PV[i] = g_CUR_PV[i];
            g_TREE_PV[i] = '0';
        }
        set_PV(i, 0);
    }
    // for tree "1", ID is 1.
    read_output_pv("1");
    for (int i = 0; i < g_ALIGN_LEN; i++) {
        if (g_TEMP_PV[i] < g_CUR_PV[i])
            g_CUR_PV[i] = g_TEMP_PV[i];
        else
            g_TREE_PV[i] = '1';
    }
    delete [] g_TEMP_PV;
    cout << "identified pvalues" << endl;
    // Identify the worst set that are branches having pvalue == curr_pv at each position.
    // memory allocate
    g_isWorstSet = new bool*[branch_num + 1];
    for (int s = 0; s < (branch_num + 1); s++)
        g_isWorstSet[s] = new bool[g_ALIGN_LEN];
    for (int i = 0; i < g_ALIGN_LEN; i++)
        init_WorstSet(i);
    read_output_worstbranch("100");
    read_output_worstbranch("0.01");
    read_output_worstbranch("1");
    cout << "identified g_isWorstSet" << endl;
    // Bonferonni Correction for  the curr_pv (due to three trees different)
    for (int i = 0; i < g_ALIGN_LEN; i++) {
        float f = get_PV(i);
        if (f > 0.33)
            set_PV(i, 1);
        else
            set_PV(i, f * 3);
    }
  
    int **count = new int*[branch_num + 1];
    for (int i = 0; i < 3; ++i)
        count[i] = new int[3];
    for (int sp = 0; sp < branch_num + 1; sp++) {
        count[sp][0] = count[sp][1] = count[sp][2] = 0;
        for (int i = 0; i < g_ALIGN_LEN; i++)
            if (sp < branch_num) {
                if (get_G_IsPresent(i, sp))
                    count[sp][0]++;
            }
            else {
                if (!is_NoneG_IsPresent(i)) 
                    count[sp][0]++;
            }
    }
    
    // memory allocate
    for (int s = 0; s < branch_num; s++)
        delete [] g_isPresent[s];
    delete [] g_isPresent;
    g_isGaps = new bool*[branch_num];
    for (int s = 0; s < branch_num; s++)
        g_isGaps[s] = new bool[g_ALIGN_LEN];
    // Reading gaps
    read_Gaps();
    cout << "read gaps" << endl;
    // map alignment coordinates to reference coordinates
    refSpeciesCoord = new int[g_ALIGN_LEN];
    map_referenceCoordinates();
    cout << "got reference coordinates map" << endl;
    // Identigying good regions
    output_good_regions();
    // Identifying bad regions
    int start;
    for (int sp = 0; sp < branch_num; sp++) {
        start = -1;
        for (int i = 0; i < g_ALIGN_LEN; i++) {
            // when current position is a bad one
            if ((get_PV(i) >= g_options.pvThreshBad) && (isin_WorstSet(i, sp))) {
                // starting position
                if (start == -1)
                    start = i;
            }
            // when current position is a good one, output the region before it if the region is long enough.
            else if (start != -1) {
                int currLen = i - start;
                if (currLen >= g_options.minSegSize) {
                    // count the number of gaps in the region.
                    int gap_count = 0;
                    for (int j = start; j < i; j++)
                        if (get_Gap(j, sp))
                            gap_count++;
                    // the first and the last reference coordinates in current region
                    int c = start;
                    while ((refSpeciesCoord[c] == -1) && (c < g_ALIGN_LEN - 1))
                        c++;
                    int firstHgCoor = refSpeciesCoord[c];
	  
                    c = start + currLen - 1;
                    while ((refSpeciesCoord[c] == -1) && (c > 0))
                        c--;
                    int lastHgCoor = refSpeciesCoord[c];
                    // output
                    if (gap_count < 0.5 * currLen) {
                        // ucsc browser coordinate start with 1, refSpeciesCoord + 1
                        cout << "my_bad_region " << start << " " << sp << " " << currLen 
                             << " " << g_refSpeciesChr << " " << firstHgCoor + 1 << " " 
                             << lastHgCoor + 1 << endl; 
                        count[sp][1] += currLen;
                        count[sp][2]++;
                        for (int my_pos = start; my_pos < i; my_pos++) {
                            // using 'branch_num' in the 'g_isWorstSet' to indicate this position is a bad one.
                            add_WorstBranch(my_pos, branch_num);
                        }
                    }
                    start = -1;
                }
                // reset the starting position, if the region before current good position is not long enough.
                else
                    start = -1;
            }
        }
    }
    // Compute stats for bad positions due to any branch.
    start = -1;
    for (int i = 0; i < g_ALIGN_LEN; i++) {
        if (isin_WorstSet(i, branch_num)) {
            if (start == -1)
                start = i;
        }
        else if (start != -1) {
            count[branch_num][1] += i - start;
            count[branch_num][2]++;
            start = -1;
        }
    }
    // Output summary
    for (int sp = 0; sp < branch_num + 1; sp++)
        printf("my_count %d %d %d(%4.2f%%) %d\n", sp, count[sp][0], count[sp][1], 
               count[sp][1] * 100.0 / float(count[sp][0]), count[sp][2]);
    cout << "printed count" << endl; 
    
    for (int i = 0; i < branch_num + 1; ++i)
        delete [] count[i];
    delete [] count;
    exit(EXIT_SUCCESS);
}
