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
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <assert.h>
#include <string.h>

using namespace std;
using namespace __gnu_cxx;

#define MAX_LENGTH_NEWICK 4096
#define MAX_LENGTH_REF_SPECIES 32

struct globalOptions_t{
    // parameters for identifying suspicious or good regions:
    //
    // the min size of suspicious or good regions
    int MIN_SEG_SIZE;
    // the pvalue threshold for suspicious regions
    float PVTHRESH_BAD;
    // the pvalue threshold for good regions
    float PVTHRESH_GOOD;
    char *mafFile;
    char *SigMAwOutDir;
    char *SigMAwOutPrefix;
    char *PHYLOGENY;
    char *REF_SPECIES;
    int verbose;
};
struct globalOptions_t globalOptions = {50, 0.1, 1e-10, (char *) malloc(100),
                                        (char*) malloc(100), (char*) malloc(20),
                                        (char*) malloc(MAX_LENGTH_NEWICK), 
                                        (char*) malloc(MAX_LENGTH_REF_SPECIES), 0};

// parameters for identifying suspicious or good regions:
//
// the min size of suspicious or good regions
int MIN_SEG_SIZE = 50;
// the pvalue threshold for suspicious regions
float PVTHRESH_BAD = 0.1;
// the pvalue threshold for good regions
float PVTHRESH_GOOD = 1e-10;


// ****************** No change below *********************

int chr_len;
int CHR_START;
int ALIGN_LEN;
int MAX_BLOCK_SIZE = 1000000;

// the .maf file
char mafFile[100] = "";
// the directory that contains output files generated by StatSigMA
char SigMAwOutDir[100] = "";
// the prefix of StatSigMA output files
char SigMAwOutPrefix[20] = "";

float* cur_pv;
float* temp_pv;
char* tree_pv;

// human chr
char hgChr[128] = "";
// human coordinates for each alignment column
int* hgCoor;

// arrays of presence and gaps at each position for each branch.
bool** present;
bool** gaps;
// the worst set of species
bool** worst_set;
