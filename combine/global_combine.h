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
#ifndef GLOBAL_COMBINE_H_
#define GLOBAL_COMBINE_H_

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <assert.h>
#include <string.h>

using namespace std;
using namespace __gnu_cxx;

const int kMaxLengthNewick = 4096;
const int kMaxLengthRefSpecies = 32;
const int kMaxMessageLength = 1024;

struct globalOptions_t{
    // parameters for identifying suspicious or good regions:
    int minSegSize; // the min size of suspicious or good regions
    float pvThreshBad; // the pvalue threshold for suspicious regions
    float pvThreshGood; // the pvalue threshold for good regions
    char *mafFile;
    char *SigMAwOutDir;
    char *SigMAwOutPrefix;
    char *phylogeny;
    char *refSpecies;
    int verbose;
    int printAll;
    bool debug;
};
struct globalOptions_t g_options = {50, 0.1, 1e-10, new char[100],
                                    new char[100], new char[20],
                                    new char[kMaxLengthNewick], 
                                    new char[kMaxLengthRefSpecies], 0, 0, false};

int g_CHR_START;
int g_ALIGN_LEN;
int g_MAX_BLOCK_SIZE = 1000000;

float *g_CUR_PV;
float *g_TEMP_PV;
char *g_TREE_PV;

char g_refSpeciesChr[128] = ""; // reference species chromosome name
int *refSpeciesCoord; // reference species coordinates for each alignment column

// arrays of presence and gaps at each position for each branch.
bool **g_isPresent;
bool **g_isGaps;
bool **g_isWorstSet; // the worst set of species

#endif // GLOBAL_COMBINE_H_
