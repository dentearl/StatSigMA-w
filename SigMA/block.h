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
#ifndef BLOCK_H_
#define BLOCK_H_

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include "alpha.h"
#include "common.h"
#include "global.h"

using namespace std;
using namespace __gnu_cxx;

int name2speciesInt(char *name) {
    // take a char* and return the index in g_SPECIES_NAMES that maches, if any
    // return -1 if no match.
    for (int i = 0; i < g_NUM_SPECIES; ++i) {
        if (strncmp(name, g_SPECIES_NAMES[i], strlen(g_SPECIES_NAMES[i])) == 0) {
            // note that g_SPECIES_NAMES is populated from the newick, so using the length
            // is important as the maf species name field will likely contain chromosome
            // extensions
            return i;
        }
    }
    return -1;
}
int sum(bool *a) {
    // take an array of bools and return the sum of their values.
    int s = 0;
    for (unsigned i = 0; i < sizeof(a) / sizeof(a[0]); ++i) {
        s += a[i];
    }
    return s;
}
class single_block {
 public:
    int **block;
    int length;
    int start;
    int align_start;
    int reference;
    list<double> block_score;
    single_block() {
        length = 0;
        start = 0;
        align_start = 0;
        reference = -1;
        block = new int*[g_options.maxBlockSize];
        for (int i = 0; i < g_options.maxBlockSize; i++) {
            block[i] = new int[g_NUM_SPECIES];
            for (int j = 0; j < g_NUM_SPECIES; j++)
                block[i][j] = kRandom;
        }
    }
    ~single_block() {
        for (int i = 0; i < g_options.maxBlockSize; i++)
            delete [] block[i];
        delete [] block;
    }
};
class blocks_type {
 public:
    vector<single_block*> chain;
    int totalLength;
    bool containsReference;
    unsigned lineNumber; // line number within the maf
    unsigned currentPosition; // current position on the reference species chromosome
    blocks_type() {
        totalLength = 0;
        containsReference = false;
        lineNumber = 0;
    }
    bool read_from_stream(char* filename);
    bool read_single_from_maf(ifstream& ifs);
};
bool blocks_type::read_single_from_maf(ifstream& ifs) {
    /* Reads one block out of the maf file pointed to by ifs and
     * stores it internally to the object
     */
    // printf("read_single_from_maf()\n");
    single_block *seq1 = NULL;
    single_block *seq2 = NULL;
    if (chain.size() == 0) {
        seq1 = new single_block;
        seq2 = new single_block;
        chain.push_back(seq1);
        chain.push_back(seq2);
    } else {
        seq1 = chain[0];
        seq2 = chain[1];
        for (int i = 0; i < g_options.maxBlockSize; i++) {
            for (int j = 0; j < g_NUM_SPECIES; j++) {
                seq1->block[i][j] = seq2->block[i][j];
            }
        }
        seq1->align_start += seq1->length;
        seq1->length = seq2->length;
        seq1->block_score = seq2->block_score;
        seq1->start = seq2->start;
        seq1->reference = seq2->reference;
    }
    totalLength = g_CHR_LEN; // Every time make block_length = g_CHR_LEN (used to compute statistics)
    size_t bufferLen = g_options.maxBlockSize + kMinLineLength;
    char *buffer = new char[bufferLen];
    buffer[0] = '\0';
    bool isSuccess = false;
    int pos = -1;
    do {
        for (int i = 0; i < g_options.maxBlockSize; i++) {
            for (int j = 0; j < g_NUM_SPECIES; j++) {
                seq2->block[i][j] = kRandom;
            }
        }
        seq2->reference = -1;
        seq2->length = 0;
        seq2->start = 0;
        // advance to the '^a' line
        while (ifs && (buffer[0] != 'a')) {
            buffer[0] = '\0';
            ifs.getline(buffer, bufferLen);
            ++lineNumber;
            if (strlen(buffer) == bufferLen - 1) {
                cerr << "This is a problem, buffer is full " << bufferLen - 1 << endl;
                exit(EXIT_FAILURE);
            }
            if (ifs.fail() && !ifs.eof()) {
                cerr << "Error, failbit has been set while reading line number " << lineNumber << endl;
                exit(EXIT_FAILURE);
            }
        }
        char* str1 = strtok(buffer, "=");
        if (str1) {
            str1 = strtok(NULL, "=");
            seq2->block_score.clear();
            seq2->block_score.push_back(atof(str1));
        }
        buffer[0] = '\0';
        ifs.getline(buffer, bufferLen);
        ++lineNumber;
        if (ifs.fail() && !ifs.eof()) {
            cerr << "Error, failbit has been set while reading line number " << lineNumber << endl;
            exit(EXIT_FAILURE);
        }
        containsReference = false;
        bool *speciesPresent = new bool[g_NUM_SPECIES]();
        while ((ifs) && (buffer[0] != 'a')) {
            if (buffer[0] == 's') {
                // sequence line business
                int sp = name2speciesInt(buffer + 2); // ptr offset
                if (sp == -1) {
                    buffer[0] = '\0'; // throw this line away, it does not contain a sequence from the newick
                    continue;
                }
                if (!speciesPresent[sp]) {
                    speciesPresent[sp] = true;
                    // cout << "read " << g_SPECIES_NAMES[sp] << endl;
                } else {
                    cerr << "Error, a block in the maf file (near line number " << lineNumber
                         << ") contains a duplicate species: " << g_SPECIES_NAMES[sp]
                         << ", which is species index: " << sp << endl;
                    exit(EXIT_FAILURE);
                }
                if (0 <= sp) {
                    char* str = strtok(buffer, " "); // s
                    str = strtok(NULL, " "); // name
                    str = strtok(NULL, " "); // start
                    int start = atoi(str);
                    if (strncmp(buffer + 2, g_options.refSpecies, strlen(g_options.refSpecies)) == 0) {
                        seq2->start = start;
                        seq2->reference = sp;
                        containsReference = true;
                        isSuccess = true;
                        if (seq2->start < seq1->start + seq1->length) {
                            cerr << "Error, detected a block where the reference either "
                                "contains a duplicated position, or is out of order, near line "
                                "number " << lineNumber << ". Prev: " << seq1->start + seq1->length <<
                                ", vs current: " << seq2->start << endl;
                            exit(EXIT_FAILURE);
                        }
                    }
                    str = strtok(NULL, " "); // length
                    str = strtok(NULL, " "); // strand
                    str = strtok(NULL, " "); // source length
                    str = strtok(NULL, " "); // alignment field
                    int len = strlen(str);
                    for (int m = 0; m < len; m++)
                        seq2->block[m][sp] = alpha2int(toupper(str[m]));
                    seq2->length = len;
                }
            }
            buffer[0] = '\0';
            // record the position of pointer
            pos = ifs.tellg();
            ifs.getline(buffer, bufferLen);
            ++lineNumber;
            if (ifs.fail() && !ifs.eof()) {
                cerr << "Error, failbit has been set while reading line number " << lineNumber << endl;
                exit(EXIT_FAILURE);
            }
        }
        if (!containsReference && sum(speciesPresent)) {
            cerr << "Error, a block in the maf (preceding line number " << lineNumber
                 << ") fails to contain the reference species." << endl;
            exit(EXIT_FAILURE);
        }
        int chars = 0;
        for (int i = 0; i < seq1->length; i++) {
            if ((seq1->reference != -1) && (seq1->block[i][seq1->reference] != kGap)) {
                chars++;
            }
        }
        debug("Reading sc=%d\t"
              "old_non_gap_char=%d\t"
              "new_seq_start=%d\t"
              "new_seq_length=%d\t"
              "old_ref=%d\t"
              "new_ref=%d\n", 
              *(seq2->block_score.begin()), 
              chars, 
              seq2->start,
              seq2->length, 
              seq1->reference, 
              seq2->reference);
        if (((chars == 0) || (chars + seq1->start == seq2->start))
            && (seq1->length + seq2->length < g_options.maxBlockSize)
            && (0 < seq2->start)) {
            for (int i = 0; i < seq2->length; i++) {
                for (int j = 0; j < g_NUM_SPECIES; j++) {
                    seq1->block[i + seq1->length][j] = seq2->block[i][j];
                }
            }
            seq1->length += seq2->length;
            seq1->block_score.push_back(*(seq2->block_score.begin()));
            seq1->reference = seq2->reference;
            if (chars == 0)
                seq1->start = seq2->start;
            seq2->length = 0;
            debug("Joined to the previous segment with total length = %d\n", seq1->length);
        }
    } while (ifs && (seq2->length == 0));
    if (ifs && (pos != -1)) {
        // reset the pointer back one line (starting with 'a')
        // get ready for reading the next SigMA block
        ifs.seekg(pos);
        --lineNumber;
    }
    return isSuccess;
}

bool blocks_type::read_from_stream(char* filename) {
    ifstream ifs(filename);
    size_t bufferLen = g_options.maxBlockSize + kMinLineLength;
    char *buffer = new char[bufferLen];
    bool isSuccess = false;
    if (strstr(filename, ".maf")) {
        debug("Reading .maf file %s\n", filename);
        while (ifs) {
            buffer[0] = '\0';
            while (ifs && (ifs.getline(buffer, bufferLen))
                   && (buffer[0] != 's'))
                continue;
            single_block* seq = new single_block;
            while ((ifs) && (buffer[0] == 's')) {
                int sp = name2speciesInt(buffer + 2); // ptr offset
                isSuccess = true;
                if (0 <= sp) {
                    char* str = strtok(buffer, " "); // s
                    str = strtok(NULL, " "); // name
                    str = strtok(NULL, " "); // start
                    // int start = atoi(str);
                    str = strtok(NULL, " "); // length
                    str = strtok(NULL, " "); // strand
                    str = strtok(NULL, " "); // source length
                    str = strtok(NULL, " "); // sequence
                    int len = strlen(str);
                    for (int m = 0; m < len; m++)
                        seq->block[m][sp] = alpha2int(toupper(str[m]));
                    seq->length = len;
                }
                buffer[0] = '\0';
                ifs.getline(buffer, bufferLen);
                if (ifs.fail() && !ifs.eof()) {
                    cerr << "Error, failbit has been set" << endl;
                    exit(EXIT_FAILURE);
                }
            }
            if (0 < seq->length) {
                chain.push_back(seq);
                totalLength += seq->length;
            }
        }
        return isSuccess;
    } else if (strstr(filename, ".fa")) {
        debug("Reading .fasta file %s\n", filename);
        int sp = -1;
        single_block* seq = new single_block;
        while (ifs) {
            ifs.getline(buffer, bufferLen);
            if (ifs.fail() && !ifs.eof()) {
                cerr << "Error, failbit has been set" << endl;
                exit(EXIT_FAILURE);
            }
            if (buffer[0] == '>') {
                char* a = strtok(buffer + 1, " ");
                while(!a)
                    a = strtok(NULL, " ");
                sp = name2speciesInt(a);
                seq->length = 0;
            }
            else if (sp>=0) {
                isSuccess = true;
                for (unsigned i = 0; i < strlen(buffer); i++)
                    seq->block[seq->length + i][sp] = alpha2int(toupper(buffer[i]));
                seq->length+=strlen(buffer);
            }
            buffer[0] = '\0';
        }
        if (seq->length > 0) {
            chain.push_back(seq);
            totalLength += seq->length;
        }
        return isSuccess;
    } else if (strstr(filename, ".msf")) {
        debug("Reading .msf file %s\n", filename);
        int cur_length = 0;
        int firstspecies=-1;
        single_block* seq = new single_block;
        while (ifs) {
            int sp = -1;
            ifs.getline(buffer, bufferLen);
            if (ifs.fail() && !ifs.eof()) {
                cerr << "Error, failbit has been set" << endl;
                exit(EXIT_FAILURE);
            }
            char* name = strtok(buffer, " ");
            if (name)
                sp = name2speciesInt(name);
            if (sp >= 0) {
                isSuccess = true;
                if (firstspecies == -1)
                    firstspecies = sp;
                else if (firstspecies == sp)
                    seq->length += cur_length;
                int pos = 0;
                int startpos = strlen(name)+2;
                for (unsigned i = 0; i < strlen(buffer + startpos); i++)
                    if ((buffer[i + startpos] == '.') || (isalpha(buffer[i + startpos]))) {
                        seq->block[seq->length + pos][sp] = alpha2int(toupper(buffer[i + startpos]));
                        pos++;
                    }
                cur_length = pos;
            }
        }
        seq->length += cur_length;
        if (seq->length > 0) {
            chain.push_back(seq);
            totalLength += seq->length;
        }
        return isSuccess;
    } else {
        cerr << "Alignment file format not identified" << endl;
        exit(EXIT_FAILURE);
    }
}

#endif // BLOCK_H_

