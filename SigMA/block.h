/*VI. Conditions of Use
 *
 * Xiaoyu Chen, Amol Prakash, and Martin Tompa give permission for Dent Earl and his institution to use the StatSigmMA - w software developed at the University of Washington for research purposes, on the following conditions:
 *
 * 1. Dent Earl and his institutional colleagues may modify the StatSigmMA - w software and distribute the resulting modified software for research purposes, provided (a) that it is distributed  together with these Conditions of Use, (b) that Martin Tompa receive a copy of the finalized modified software, and (c) that Xiaoyu Chen, Amol Prakash, and Martin Tompa are credited with the authorship of the software.
 *  2. The StatSigmMA - w software will be used by you and/or your institution solely for noncommercial purposes, except with express permission from the authors.
 * 3. Any risk associated with using the StatSigmMA - w software at your institution is with you and your institution.
 * 4. StatSigmMA - w will be cited in any publication(s) reporting on data obtained from it as:
 *    Amol Prakash and Martin Tompa, "Measuring the Accuracy of Genome - Size Multiple Alignments". Genome Biology, vol. 8, issue 6, June 2007, R124.
 *    Xiaoyu Chen and Martin Tompa, "Comparative assessment of methods for aligning multiple genome sequences". Nature Biotechnology, vol. 28, no. 6, June 2010, 567 - 572.
 */
#ifndef BLOCK_H_
#define BLOCK_H_

#include "alpha.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string.h>
#include <stdlib.h>
#include "common.h"

using namespace std;
using namespace __gnu_cxx;

int name2species(char* name) { 
    for (int i = 0; i < g_NUM_SPECIES; i++)
        if (strncmp(name, g_SPECIES_NAMES[i], strlen(g_SPECIES_NAMES[i])) == 0)
            return i;
    return -1;
}
class single_block {
 public:
    int** block;
    int length;
    int start;
    int align_start;
    int human_ref;
    list<double> block_score;
    single_block() {
        length = 0;
        start = 0;
        align_start = 0;
        human_ref = -1;
        block = new int*[g_options.maxBlockSize];
        for (int i = 0; i < g_options.maxBlockSize; i++) {
            block[i] = new int[g_NUM_SPECIES];
            for (int j = 0; j < g_NUM_SPECIES; j++)
                block[i][j] = c_RANDOM;
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
    int total_length;
    blocks_type() {total_length = 0;}
    bool read_from_stream(char* filename);
    bool read_single_from_maf(ifstream& ifs);
};
bool blocks_type::read_single_from_maf(ifstream& ifs) {
    single_block* seq1;
    single_block* seq2;
    if (chain.size() == 0) {
        seq1 = new single_block;
        chain.push_back(seq1);
        seq2 = new single_block;
        chain.push_back(seq2);
    } else {
        seq1 = chain[0];
        seq2 = chain[1];
        for (int i = 0; i < g_options.maxBlockSize; i++)
            for (int j = 0; j < g_NUM_SPECIES; j++)
                seq1->block[i][j] = seq2->block[i][j];
        seq1->align_start += seq1->length;
        seq1->length = seq2->length;
        seq1->block_score = seq2->block_score;
        seq1->start = seq2->start;
        seq1->human_ref = seq2->human_ref;
    }
    total_length = g_CHR_LEN; // Every time make block_length = chr_len (used to compute statistics)
    char *buffer = new char[g_options.maxBlockSize + 1];
    buffer[0] = '\0';
    bool isSuccess = 0;
    int pos = -1;
    do {
        for (int i = 0; i < g_options.maxBlockSize; i++)
            for (int j = 0; j < g_NUM_SPECIES; j++)
                seq2->block[i][j] = c_RANDOM;
        seq2->human_ref = -1;
        seq2->length = 0;
        seq2->start = 0;
        while (ifs && (buffer[0] != 'a')) {
            buffer[0] = '\0';
            ifs.getline(buffer, g_options.maxBlockSize);
        }
        char* str1 = strtok(buffer, "=");
        if (str1) {
            str1 = strtok(NULL, "=");
            seq2->block_score.clear();
            seq2->block_score.push_back(atof(str1));
        }
        buffer[0] = '\0';
        ifs.getline(buffer, g_options.maxBlockSize);
        while ((ifs) && (buffer[0] != 'a')) {
            if (buffer[0] == 's') {
                int sp = name2species(buffer + 2);
                isSuccess = 1;
                if (0 <= sp) {
                    char* str = strtok(buffer," "); // s
                    str = strtok(NULL," "); // name
                    str = strtok(NULL," "); // start
                    int coord = atoi(str);
                    if (strncmp(buffer + 2, g_options.refSpecies, strlen(g_options.refSpecies)) == 0) {
                        seq2->start = coord;
                        seq2->human_ref = sp;
                    }
                    str = strtok(NULL," "); // length
                    str = strtok(NULL," "); // strand
                    str = strtok(NULL," "); // source length
                    str = strtok(NULL," "); // alignment field
                    int len = strlen(str);
                    for (int m = 0; m < len; m++)
                        seq2->block[m][sp] = alpha2int(toupper(str[m]));
                    seq2->length = len;
                }
            }
            buffer[0] = '\0';
            // record the position of pointer
            pos = ifs.tellg();
            ifs.getline(buffer, g_options.maxBlockSize);
        }
        int chars = 0;
        for (int i = 0; i < seq1->length; i++)
            if ((seq1->human_ref != -1) && (seq1->block[i][seq1->human_ref] != c_GAP))
                chars++;
        debug("# Reading sc=%d\t"
              "old_non_gap_char=%d\t"
              "new_seq_start=%d\t"
              "new_seq_length=%d\t"
              "old_human=%d\t"
              "new_human=%d\n", *(seq2->block_score.begin()), chars, seq2->start, 
              seq2->length, seq1->human_ref, seq2->human_ref);
        if (((chars == 0) || (chars + seq1->start == seq2->start)) 
            && (seq1->length + seq2->length < g_options.maxBlockSize)
            && (0 < seq2->start)) {
            for (int i = 0; i < seq2->length; i++)
                for (int j = 0; j < g_NUM_SPECIES; j++)
                    seq1->block[i + seq1->length][j] = seq2->block[i][j];
            seq1->length += seq2->length;
            seq1->block_score.push_back(*(seq2->block_score.begin()));
            seq1->human_ref = seq2->human_ref;
            if (chars == 0)
                seq1->start = seq2->start;
            seq2->length = 0;
            debug("# Joined to the previous segment with total length = %d\n", seq1->length);
        }
    } while (ifs && (seq2->length == 0));
    if (ifs && (pos != -1)) {
        // reset the pointer back a line (starting with 'a')
        // get ready for reading the next SigMA block
        ifs.seekg(pos);
    }
    return isSuccess;
}

bool blocks_type::read_from_stream(char* filename) {
    ifstream ifs(filename);
    char *buffer = new char[g_options.maxBlockSize + 1];
    bool isSuccess = 0;
    if (strstr(filename,".maf")) {
        debug("# Reading .maf file %s\n", filename);
        while (ifs) {
            buffer[0]='\0';
            while (ifs && (ifs.getline(buffer, g_options.maxBlockSize))
                   &&(buffer[0] != 's'))
                continue;
            single_block* seq = new single_block;
            while ((ifs) && (buffer[0] == 's')) {
                int sp = name2species(buffer + 2);
                isSuccess = 1;
                if (0 <= sp) {
                    char* str = strtok(buffer," ");
                    str = strtok(NULL," ");
                    str = strtok(NULL," ");
                    // int start = atoi(str);
                    str = strtok(NULL," ");
                    str = strtok(NULL," ");
                    str = strtok(NULL," ");
                    str = strtok(NULL," ");
                    int len = strlen(str);
                    for (int m = 0;m<len;m++)
                        seq->block[m][sp]=alpha2int(toupper(str[m]));
                    seq->length = len;
                }
                buffer[0]='\0';
                ifs.getline(buffer, g_options.maxBlockSize);
            }
            if (0 < seq->length) {
                chain.push_back(seq);
                total_length += seq->length;
            }
        }
    
        return isSuccess;
    } else if (strstr(filename,".fa")) {
        debug("# Reading .fasta file %s\n", filename);
        int sp = -1;
        single_block* seq = new single_block;
        while (ifs) {
            ifs.getline(buffer, g_options.maxBlockSize);
            if (buffer[0] == '>') {
                char* a = strtok(buffer + 1, " ");
                while(!a)
                    a = strtok(NULL, " ");
                sp = name2species(a);
                seq->length = 0;
            }
            else if (sp>=0) {
                isSuccess = 1;
                for (unsigned i = 0; i < strlen(buffer); i++)
                    seq->block[seq->length + i][sp] = alpha2int(toupper(buffer[i]));
                seq->length+=strlen(buffer);
            }
            buffer[0] = '\0';
        }
        if (seq->length > 0) {
            chain.push_back(seq);
            total_length += seq->length;
        }
        return isSuccess;
    } else if (strstr(filename, ".msf")) {
        debug("# Reading .msf file %s\n", filename);
        int cur_length = 0;
        int firstspecies=-1;
        single_block* seq = new single_block;
        while (ifs) {
            int sp=-1;
            ifs.getline(buffer, g_options.maxBlockSize);
            char* name = strtok(buffer, " ");
            if (name)
                sp = name2species(name);
            if (sp >= 0) {
                isSuccess = 1;
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
            total_length += seq->length;
        }
        return isSuccess;
    } else {
        cerr << "Alignment file format not identified" << endl;
        exit(EXIT_FAILURE);
    }
}


#endif // BLOCK_H_

