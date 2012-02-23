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
// #define DEBUG

#ifndef GLOBAL_H_
#define GLOBAL_H_

const int kMaxFilePath = 1024;
const int kMaxLengthNewick = 4096;
const int kMaxLengthRefSpecies = 32;
const int kMaxLengthMessage = 1024;
const int kMinLineLength = 1024;

struct g_options_t{
    char *mafFile;
    int branchIndex;
    char *multipvFile;
    double branchMultiplier;
    // maximum total size of contiguous alignment blocks allowed.
    // If the length of one alignment block is larger than value assigned below,
    // MAX_BLOCK_SIZE will be assigned to the actual alignment block size.
    int maxBlockSize;
    // max number of high-scoring segments
    int maxSegments;
    // These are used to estimate Karlin-Altschul parameters:
    // max number of tuples to be read per block (after randomization)
    int totalNumberTuples;
    // number of times that Karlin-Altschul parameter is computed 
    int totalIterateParam;
    char *phylogeny;
    char *refSpecies;
    int verbose;
    bool debug;
    unsigned rseed; // random seed.
};
g_options_t g_options = {new char[kMaxFilePath],
                         0, new char[kMaxFilePath], 1.0,
                         100000, 100000, 50000, 100,
                         new char[kMaxLengthNewick],
                         new char[kMaxLengthRefSpecies], 0, false, 0};

int g_CHR_LEN; // chromosome length
int g_CHR_START; // The start chromosome position of alignment
int g_ALIGN_LEN; // the length of alignment

// Parameters for the divide & conquer algo, 
// which estimates a pvalue for each segment:
//

float g_PV_THRESH_DVDCQ = 1e-5; // pvalue threshold
const int kMaxNumberContexts = 10000; // the max number of contexts 

// Estimating parameters
const int kScoreBins = 1000;
const int kMaxTotalScore = 100;
const int kMinTotalScore = -100;
const float kSubstitutionPerSite = 1;
const float kGapBackground = 0.05;
const float kMaxBranchLength = 10000;
// Divide by length; this value is added to each score unit
const int kMaxNumberBranches = 1000; // largest number of branches allowed in a tree

int g_NUM_SPECIES = kMaxNumberBranches;
char g_SPECIES_NAMES[kMaxNumberBranches][32];
double** g_MULTI_PV;
int g_NUM_BLOCKS_SKIPPED = 0; // count of blocks due to missing reference

#endif // GLOBAL_H_
