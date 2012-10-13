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
#ifndef PVALUE_H_
#define PVALUE_H_
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "tree.h"
#include "block.h"
#include "blast.h"
//#include <pair.h>

using namespace std;
using namespace __gnu_cxx;

struct final_seg_type {
    double score; // segment score
    int start;
    int end;
    int align_start;
    int align_end;
    double pv;
    double pv2;
    int block_num; // Block number from reading maf
    int nearest_left; // position of left nearest significant segment
    int nearest_right; // position of right nearest significant segment
    int block_score; //  list < double> block_score;
    int left_end;
    int right_end;
    final_seg_type() {pv = -1;}
};
int cmpPairInts(const void *a, const void *b) {
    double aa = ((pair <int, int>*) a)->second;
    double bb = ((pair <int, int>*) b)->second;
    if (aa > bb) return -1;
    if (aa < bb) return 1;
    else return 0;
}
bool compute_branch_parameters(Tree *node, Tree *parent, const blocks_type& blocks,
                               double &K, double &lambda, double &H) {
    debug("Computing branch parameters for node start number: %d end number: %d\n",
          node->leaf_startnum, node->leaf_endnum);
    int max_score = 10 * kScoreBins;
    int min_score = -10 * kScoreBins;
    double *score_dist = new double[max_score - min_score + 1];
    if (score_dist == NULL) {
        cerr << "Unable to create score_dist when requesting %d doubles." << max_score - min_score + 1
             << endl;
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < (max_score - min_score + 1); ++i)
        score_dist[i] = 0;
    double num_tuples = pow(kAlphabetSize, g_NUM_SPECIES);
    int total_computations = 0;
    debug("Comparison of number of tuples, all=%d, block_len=%d\n",
          num_tuples, blocks.chain.at(0)->length);
    if (num_tuples < blocks.chain.at(0)->length * blocks.chain.at(0)->length) {
        total_computations = 10000; // default value
        debug("Computing all scores %d\n", num_tuples);
        for (int num = 0; num < num_tuples; ++num) {
            for (int k = 0; k < g_NUM_SPECIES; ++k)
                list_of_leaves[k]->nucleotide = ((static_cast<int>(num / pow(kAlphabetSize, g_NUM_SPECIES - k - 1))) 
                                                 % (kAlphabetSize) + 1);
            root->compute_prob();
            double ortho_pv = root->return_prob();
            double non_ortho_pv;
            if (parent->max_prob() == kGap)
                non_ortho_pv = ortho_pv;
            else {
                node->non_ortho = true;
                root->compute_prob();
                non_ortho_pv = root->return_prob();
                node->non_ortho = false;
            }
            int score = 0;
            if ((ortho_pv == 0) || (non_ortho_pv == 0))
                score = 0;
            else
                score = static_cast<int>(round(log(ortho_pv / non_ortho_pv) * kScoreBins));
            if (score < min_score)
                score = min_score;
            if (score > max_score)
                score = max_score;
            score_dist[score - min_score] += non_ortho_pv;
        }
    } else {
        /*
          for (int i = 0; i < max_score - min_score + 1; ++i)
          score_dist[i] += double(PSEUDOCOUNT_NORMALISATION)/(blocks.chain[0]->length*(max_score - min_score));
        */
        debug("Computing sequence specific scores max(N^2, %d) where "
              "N=%d (again based on random sampling)\n", g_options.totalNumberTuples,
              blocks.chain[0]->length);
        int bl1 = 0;
        int bl2 = 0;
        int* good_tuple1 = new int[g_options.maxBlockSize];
        int* good_tuple2 = new int[g_options.maxBlockSize];
        if (good_tuple1 == NULL) {
            cerr << "Unable to allocate memory to create good_tuple1." << endl;
            exit(EXIT_FAILURE);
        }
        if (good_tuple2 == NULL) {
            cerr << "Unable to allocate memory create good_tuple2." << endl;
            exit(EXIT_FAILURE);
        }
        int good1 = 0, good2 = 0;
        for (int i = 0; i < g_options.maxBlockSize; ++i)
            good_tuple1[i] = good_tuple2[i] = -1;
        for (int l1 = 0; l1 < blocks.chain[bl1]->length; ++l1) {
            for (int i = 0; i < node->leaf_startnum; ++i) {
                if (blocks.chain[bl1]->block[l1][i] != kRandom) {
                    good_tuple1[good1] = l1;
                }
            }
            for (int i = node->leaf_endnum + 1; i < g_NUM_SPECIES; ++i) {
                if (blocks.chain[bl1]->block[l1][i] != kRandom) {
                    good_tuple1[good1] = l1;
                }
            }
            if (good_tuple1[good1] == l1) {
                ++good1;
            }
        }
        for (int l2 = 0; l2 < blocks.chain[bl2]->length; ++l2) {
            for (int i = node->leaf_startnum; i <= node->leaf_endnum; ++i) {
                if (blocks.chain[bl2]->block[l2][i] != kRandom) {
                    good_tuple2[good2] = l2;
                }
            }
            if (good_tuple2[good2] == l2) {
                ++good2;
            }
        }
        debug("Good1 %d\n", good1);
        debug("Good2 %d\n", good2);
        if ((good1 > 0) && (good2 > 0)) {
            for (int tp = 0; tp < g_options.totalNumberTuples; ++tp) {
                // random choose a tuple from good_tuple1 and good_tuple2
                int goodID1 = static_cast<int>(rand() / (((double)RAND_MAX + 1) / good1));
                int goodID2 = static_cast<int>(rand() / (((double)RAND_MAX + 1) / good2));
                int l1 = good_tuple1[goodID1];
                int l2 = good_tuple2[goodID2];
                if ((total_computations + 1) % 10000 == 0)
                    debug("Tuples read = %d\n",total_computations + 1);
                ++total_computations;
                for (int i = 0; i < node->leaf_startnum; ++i)
                    list_of_leaves[i]->nucleotide = blocks.chain[bl1]->block[l1][i];
                for (int i = node->leaf_startnum; i<=node->leaf_endnum; ++i)
                    list_of_leaves[i]->nucleotide = blocks.chain[bl2]->block[l2][i];
                for (int i = node->leaf_endnum + 1; i < g_NUM_SPECIES; ++i)
                    list_of_leaves[i]->nucleotide = blocks.chain[bl1]->block[l1][i];
                root->compute_prob();
                double ortho_pv = root->return_prob();
                double non_ortho_pv;
                if (parent->max_prob() == kGap)
                    non_ortho_pv = ortho_pv;
                else {
                    node->non_ortho = true;
                    root->compute_prob();
                    non_ortho_pv = root->return_prob();
                    node->non_ortho = false;
                }
                int score = 0;
                if ((ortho_pv == 0) || (non_ortho_pv == 0))
                    score = 0;
                else
                    score = static_cast<int>(log(ortho_pv / non_ortho_pv) * kScoreBins);
                // cout << "\t" << score << "  " << ortho_pv << " " << non_ortho_pv << endl;
                if (score < min_score)
                    score = min_score;
                if (score > max_score)
                    score = max_score;
                ++score_dist[score - min_score];
            }
        } else
            debug("No valid tuples: Good1=%d Good2=%d\n", good1, good2);
        delete [] good_tuple1;
        delete [] good_tuple2;
        debug("Finished score sampling\n");
        debug("%s", timeStamp());
    }
    int new_max = 0;
    int new_min = 0;
    for (int i = 0; i < max_score - min_score + 1; ++i)
        if (score_dist[i] != 0) {
            if (i + min_score > new_max)
                new_max = i + min_score;
            if (i + min_score < new_min)
                new_min = i + min_score;
        }
    for (int i = 0; i < new_max - new_min + 1; ++i)
        score_dist[i] = score_dist[i - min_score + new_min];
    debug("Finished score distribution\n");
    debug("%s", timeStamp());
    if (total_computations < 5) {
        K = 1;
        lambda = 0.0;
        H = 1;
        debug("Too few good tuples\n");
        return false;
    } else if (! karlin(new_min, new_max, score_dist, &lambda, &K, &H)) {
        cout << "Karlin Altschul estimation Failed \n";
        debug("Finished Karlin\n");
        debug("%s", timeStamp());
        return false;
    }
    //  lambda *= kScoreBins;          b'coz all scores are normalized
    debug("Parameters K=%d lambda=%f H=%f\n", K, lambda, H);
    debug("%s", timeStamp());
    delete [] score_dist;
    return true;
}
void scan_alignFile(char *mafFile) {
    /* reads through the mafFile once and sets global variables
     */
    g_CHR_START = -1;
    g_CHR_LEN = 0;
    g_ALIGN_LEN = 0;
    int tempSize = 10 * g_options.maxBlockSize;
    int last_pos = 0;
    int cur_len = 0;
    char *buffer = new char[tempSize];
    if (buffer == NULL) {
        cerr << "Error, unable to create new buffer." << endl;
        exit(EXIT_FAILURE);
    }
    bool firstBlock = true;
    int linesRead = 1;
    ifstream ifs(mafFile);
    if (!ifs.is_open()) {
        cerr << "Error, unable to open maf file " << mafFile << endl;
        exit(EXIT_FAILURE);
    }
    debug("Opening maf %s\n", mafFile);
    while (ifs) {
        buffer[0] = '\0';
        ifs.getline(buffer, tempSize);
        if (ifs.fail() && !ifs.eof()) {
            cerr << "Error, failbit has been set while scanning " << mafFile << endl;
            exit(EXIT_FAILURE);
        }
        ++linesRead;
        if ((strncmp(buffer, "a score=", 8) == 0) || (!ifs)) {
            // make sure g_options.maxBlockSize >= actual alignment block size
            if (g_options.maxBlockSize < cur_len) {
                g_options.maxBlockSize = cur_len;
            }
            g_ALIGN_LEN += cur_len;
            cur_len = 0;
            // sequence length = last position - first position + 1
            if (last_pos > 0) {
                g_CHR_LEN = last_pos - g_CHR_START + 1;
            }
        } else {
            char* a;
            if (buffer[0] == 's') {
                // Checking for reference
                bool is_reference = false;
                if (strncmp(buffer + 2, g_options.refSpecies, strlen(g_options.refSpecies)) == 0) {
                    is_reference = true;
                }
                if (!is_reference) {
                    continue;
                }
                a = strtok(buffer," "); // s
                a = strtok(NULL," "); // name
                a = strtok(NULL," "); // start
                if (firstBlock) {
                    g_CHR_START = atoi(a);
                    firstBlock = false;
                }
                last_pos = atoi(a);
                a = strtok(NULL," "); // length
                // get the last genomic position
                last_pos += atoi(a) - 1;
                a = strtok(NULL," "); // strand
                a = strtok(NULL," "); // source length
                a = strtok(NULL," "); // alignment field
                cur_len = strlen(a);
            }
        }
    }
    debug("Lines read = %d\n", linesRead);
    ifs.close();
    delete [] buffer;
}
void printBlock(blocks_type *block) {
    // function for debugging.
    printf("##############################\n");
    for (int i = 0; i < g_options.maxBlockSize; ++i) {
        for (int j = 0; j < g_NUM_SPECIES; ++j) {
            if (block->chain[0]->block[i][j] != 0) {
                printf("%d%s", block->chain[0]->block[i][j], (j == g_NUM_SPECIES - 1) ? "\n" : ", ");
            } else {
                printf(" %s", (j == g_NUM_SPECIES - 1) ? "\n" : ", ");
            }
        }
    }
    printf("##############################\n");
    for (int i = 0; i < g_options.maxBlockSize; ++i) {
        for (int j = 0; j < g_NUM_SPECIES; ++j) {
            if (block->chain[1]->block[i][j] != 0) {
                printf("%d%s", block->chain[1]->block[i][j], (j == g_NUM_SPECIES - 1) ? "\n" : ", ");
            } else {
                printf(" %s", (j == g_NUM_SPECIES - 1) ? "\n" : ", "); 
            }
        }
    }
    printf("##############################\n");
}
double compute_pvalue_branch(Tree *node, Tree *parent, blocks_type &blocks, char *mafFile) {
    /* computes the p-value for one branch in the tree
     */
    int max_seg = 6000000;
    final_seg_type *all_seg = new final_seg_type[max_seg + 1];
    scan_alignFile(mafFile);
    debug("Reference chromosome start position = %d\n", g_CHR_START);
    debug("Reference sequence(non-gap) length = %d\n", g_CHR_LEN);
    debug("Alignment length = %d\n", g_ALIGN_LEN);
    debug("Max block size = %d\n", g_options.maxBlockSize);
    double K = 0.29;
    double lambda = 0.0011;
    double H = 0.68;
    // compute_branch_parameters(node, parent, blocks, K,lambda, H);
    double *K_arr = new double[g_options.totalIterateParam];
    double *lambda_arr = new double[g_options.totalIterateParam];
    double *H_arr = new double[g_options.totalIterateParam];
    blocks_type *blocks_param = new blocks_type;
    ifstream ifs1(mafFile);
    // int runs = 0;
    float rate = g_ALIGN_LEN / static_cast<float>(g_options.totalIterateParam);
    debug("Align length = %d\n", g_ALIGN_LEN);
    int sim_count = 0;
    int block_count = 0;
    int lineNumber = -1;
    // printf("going to try to read\n");
    while ((blocks_param->read_single_from_maf(ifs1)) && (sim_count < g_options.totalIterateParam)) {
        lineNumber = blocks_param->lineNumber;
        if (!blocks_param->containsReference) {
            // reference not found in this block, skip
            ++g_NUM_BLOCKS_SKIPPED;
            continue;
        }
        int currBlockSize = blocks_param->chain[0]->length;
        int currRepNum = static_cast<int>(currBlockSize / rate);
        // if failures, need try more times
        int maxTry = currRepNum * 2;
        int blockSimCount = 0;
        for (int r = 0; r < maxTry; ++r) {	
            if (sim_count >= g_options.totalIterateParam) {
                continue;
            }
            if (compute_branch_parameters(node, parent, *blocks_param, K, lambda, H)) {
                K_arr[sim_count] = K;
                lambda_arr[sim_count] = lambda;
                H_arr[sim_count] = H;
                ++sim_count;
                ++blockSimCount;
                // get enough estimation from current block 
                if (blockSimCount == currRepNum) {
                    break;
                }
            } else {
                debug("That block did not work, try again\n");
            }
        }
        ++block_count;
    }
    // printf("sim_count = %d\n", sim_count);
    delete blocks_param;
    ifs1.close();
    cout << "## StatSigMAw " << kVersion << " " << timeStamp();
    cout << "# rseed: " << g_options.rseed << endl;
    cout << "# number of blocks skipped due to missing reference= " << g_NUM_BLOCKS_SKIPPED << endl;
    cout << "# number of sigMA blocks= " << block_count << endl;
    for (int i = 0; i < sim_count; ++i) {
        K_arr[0] += K_arr[i];
        lambda_arr[0] += lambda_arr[i];
        H_arr[0] += H_arr[i];
    }
    if (sim_count == 0) {
        cerr << "Warning, sim_count == 0, no reference found in maf block?" << endl;
        cerr << "Last read line number was " << lineNumber - 1 << endl;
        // exit(EXIT_FAILURE);
    }
    K = K_arr[0] / sim_count;
    lambda = lambda_arr[0] / sim_count;
    H = H_arr[0] / sim_count;
    delete [] K_arr;
    delete [] lambda_arr;
    delete [] H_arr;
    // do not change the format of this next out line, it is exactly checked for in combine.cpp
    cout << "# number of sets of parameter computed= " << sim_count << endl;
    int cur_seg = 0;
    ifstream ifs(mafFile); // reset mafFile
    pair <int, int> *limit = new pair <int, int>[g_options.maxSegments];
    pair <double, double> *score_left_right = new pair <double, double>[g_options.maxSegments];
    int count = 0;
    // printf("reading differently\n");
    while (blocks.read_single_from_maf(ifs) && (cur_seg < max_seg)) {
        // if (!blocks.containsReference) {
            // printf("blocks.containsReference: %s\n", blocks.containsReference ? "true" : "false");
            // continue;
        // }
        int all_good_segments = 0;
        int block_good_segments = 0;
        int bl = 0;
        double cur_score = 0;
        for (int l = 0; l < blocks.chain[bl]->length; ++l) {
            for (int i = 0; i < g_NUM_SPECIES; ++i)
                list_of_leaves[i]->nucleotide = blocks.chain[bl]->block[l][i];
            root->compute_prob();
            double ortho_pv = root->return_prob();
            double non_ortho_pv;
            int score;
            if (parent->max_prob() == kGap) {
                non_ortho_pv = ortho_pv;
                score = 0;
            } else {
                node->non_ortho = true;
                root->compute_prob();
                non_ortho_pv = root->return_prob();
                node->non_ortho = false;
                int max_score = 10 * kScoreBins;
                int min_score= -10 * kScoreBins;
                if ((ortho_pv == 0) || (non_ortho_pv == 0))
                    score = 0;
                else
                    score = static_cast<int>(round(log(ortho_pv / non_ortho_pv) * kScoreBins));
                if (score < min_score)
                    score = min_score;
                if (score > max_score)
                    score = max_score;
                // Adjusting for cutting branch for which there is no data, and approximate
                // errors in ortho_pv/non_ortho_pv which may not be exactly zero
                if (root->all_random(node) || node->all_random(NULL))
                    //	if ((score <= 0) && (score >= -50))
                    score = 10000 * min_score;
            }
            cur_score += score;
            if (score > 0) {
                int max_j = -1;
                for (int i = (block_good_segments - 1); i >= 0; --i)
                    if ((score_left_right[all_good_segments + i].first < cur_score - score) && (i > max_j))
                        max_j = i;
                if ((max_j == -1) || (score_left_right[all_good_segments + max_j].second >= cur_score)) {
                    score_left_right[all_good_segments + block_good_segments].first = cur_score - score;
                    score_left_right[all_good_segments + block_good_segments].second = cur_score;
                    limit[all_good_segments + block_good_segments].first = l;
                    limit[all_good_segments + block_good_segments].second = l;
                    if (block_good_segments + all_good_segments < g_options.maxSegments - 1)
                        ++block_good_segments;
                } else {
                    score_left_right[all_good_segments + max_j].second = cur_score;
                    limit[all_good_segments + max_j].second = l;
                    block_good_segments = max_j + 1;
                    int new_max_j = -1;
                    bool done = false;
                    do {
                        for (int i = 0; i < block_good_segments - 1; ++i)
                            if ((score_left_right[all_good_segments + i].first
                                 < score_left_right[all_good_segments + block_good_segments - 1].first)
                                && (i > new_max_j))
                                new_max_j = i;
                        if ((new_max_j == -1)
                            || (score_left_right[all_good_segments + new_max_j].second
                                >= score_left_right[all_good_segments + block_good_segments - 1].second))
                            done = true;
                        else {
                            score_left_right[all_good_segments + new_max_j].second = score_left_right[all_good_segments + block_good_segments - 1].second;
                            limit[all_good_segments + new_max_j].second = limit[all_good_segments + block_good_segments - 1].second;
                            block_good_segments = new_max_j + 1;
                        }
                    } while (!done);
	
                }
            }
        }
        //    all_good_segments += block_good_segments;
        //    block_good_segments = 0;
        for (int i = 0; i < block_good_segments; ++i)
            if (cur_seg < max_seg) {
                all_seg[cur_seg].score = score_left_right[i].second - score_left_right[i].first;
                all_seg[cur_seg].block_num = count;
                all_seg[cur_seg].block_score = blocks.chain[bl]->block_score.size();
                int chars = 0;
                for (int j = 0; j < limit[i].first; ++j) {
                    int p = blocks.chain[bl]->block[j][blocks.chain[bl]->reference];
                    if ((p >= 1) && (p <= 4))
                        ++chars;
                }
                all_seg[cur_seg].start = blocks.chain[bl]->start + chars;
                all_seg[cur_seg].align_start = blocks.chain[bl]->align_start + limit[i].first;
                for (int j = limit[i].first; j < limit[i].second; ++j) {
                    int p = blocks.chain[bl]->block[j][blocks.chain[bl]->reference];
                    if ((p >= 1) && (p <= 4))
                        ++chars;
                }
                all_seg[cur_seg].end = blocks.chain[bl]->start + chars;
                all_seg[cur_seg].align_end = blocks.chain[bl]->align_start + limit[i].second;
                ++cur_seg;
            }
        count += blocks.chain[bl]->block_score.size();
    }
    ifs.close();
    delete [] limit;
    delete [] score_left_right;
    debug("num_segments = %d\n", cur_seg);
    for (int i = 0; i < cur_seg; ++i)
        debug("Segment %d\t%f\t%d\t%d\t%d\t%d\n", i, all_seg[i].score, all_seg[i].start,
              all_seg[i].end, all_seg[i].align_start, all_seg[i].align_end);
    pair <int, int>* best_scores = new pair <int, int>[cur_seg];
    for (int i = 0; i < cur_seg; ++i) {
        best_scores[i].first = i;
        best_scores[i].second = static_cast<int>(all_seg[i].score);
    }
    qsort(best_scores, cur_seg, sizeof(best_scores[0]), cmpPairInts);
    double best_pvalue = 1.0;
    int cur_best_index = 0;
    if (cur_best_index < cur_seg) {
        do {
            int index = best_scores[cur_best_index].first;
            int totalLength = blocks.totalLength;
            double mean_length = log(K * totalLength * totalLength) / H;
            double normalized_score = (lambda * all_seg[index].score
                                       - log(K * (totalLength - mean_length)
                                             * (totalLength - mean_length)));
            if (normalized_score < kMinTotalScore)
                break;
            if (normalized_score > kMaxTotalScore) normalized_score = kMaxTotalScore;
            best_pvalue = g_MULTI_PV[1][static_cast<int>(normalized_score - kMinTotalScore)];
            all_seg[index].pv = best_pvalue;
            all_seg[index].nearest_left = index;
            all_seg[index].nearest_right = -1;
            all_seg[index].left_end = 0;
            all_seg[index].right_end = cur_seg;
            all_seg[index].pv2 = all_seg[index].pv;
            debug("pv %d\t%f\t%d\t%d\t%f\n", index, all_seg[index].score, all_seg[index].start,
                  all_seg[index].end, all_seg[index].pv);
            ++cur_best_index;
        } while ((best_pvalue < g_PV_THRESH_DVDCQ) && (cur_best_index < cur_seg));
    }
    // End coordinates
    if (cur_best_index < cur_seg) {
        do {
            int index = best_scores[cur_best_index].first;
            int find_left = index;
            while ((find_left >= 0) && (all_seg[find_left].pv == -1))
                --find_left;
            int find_right = index;
            while ((find_right < cur_seg) && (all_seg[find_right].pv == -1))
                ++find_right;
            int totalLength = blocks.totalLength;
            if (find_right < cur_seg)
                totalLength = all_seg[find_right].align_start;
            if (find_left >= 0)
                totalLength -= all_seg[find_left].align_end;
            double mean_length = log(K * totalLength * totalLength) / H;
            double normalized_score = (lambda * all_seg[index].score
                                       - log(K * (totalLength - mean_length)
                                             * (totalLength - mean_length)));
            if (normalized_score < kMinTotalScore)
                break;
            if (normalized_score > kMaxTotalScore) normalized_score = kMaxTotalScore;
            best_pvalue = g_MULTI_PV[1][static_cast<int>(normalized_score - kMinTotalScore)];
            if ((0 <= find_left) && (best_pvalue < all_seg[find_left].pv))
                best_pvalue = all_seg[find_left].pv;
            if ((find_right < cur_seg) && (best_pvalue < all_seg[find_right].pv))
                best_pvalue = all_seg[find_right].pv;
            all_seg[index].pv = best_pvalue;
            all_seg[index].nearest_left = index;
            all_seg[index].nearest_right = -1;
            all_seg[index].left_end = 0;
            all_seg[index].right_end = cur_seg;
            all_seg[index].pv2 = all_seg[index].pv;
            debug("pv %d\t%f\t%d\t%d\t%f\n", index, all_seg[index].score, all_seg[index].start,
                  all_seg[index].end, all_seg[index].pv);
            ++cur_best_index;
        } while ((best_pvalue < g_PV_THRESH_DVDCQ) && (cur_best_index < cur_seg));
        debug("out while due to best_pvalue = %f\n", best_pvalue);
    }
    pair <int, int>* intermediate_scores = new pair <int, int>[cur_seg];
    if (cur_best_index < cur_seg) {
        do {
            debug("cur_seg: %d cur_best_index: %d\n", cur_seg, cur_best_index);
            debug("%s", timeStamp());
            int index = best_scores[cur_best_index].first;
            int find_left = index;
            while ((find_left >= 0) && (all_seg[find_left].pv == -1))
                --find_left;
            debug("index %d find left: %d\n", index, find_left);
            int find_right = index;
            while ((find_right < cur_seg) && (all_seg[find_right].pv == -1))
                ++find_right;
            debug(" find right: %d\n", find_right);
            // find_left is the leftmost that is not included
            // find_right is the rightmost that is not included
            bool take_left = true;
            bool take_right = true;
            if ((0 <= find_left) && (0.01 < all_seg[find_left].pv) &&
                (all_seg[find_left].right_end <= index)) {
                find_left = all_seg[find_left].left_end;
                take_left = false;
            }
            if ((find_right < cur_seg) && (0.01 < all_seg[find_right].pv) &&
                (all_seg[find_right].left_end >= index)) {
                find_right = all_seg[find_right].right_end;
                take_right = false;
            }
            debug("left: %d %d, right: %d %d\n", find_left, take_left, find_right, take_right);
            double max_pvalue = 0;
            // subsampling kMaxNumberContexts of all possible contexts
            int leftNum = index - find_left;
            int rightNum = find_right - index;
            if (leftNum * rightNum > kMaxNumberContexts) {
                int contextCount = 0;
                while (contextCount < kMaxNumberContexts) {
                    // random number [0, index - find_left - 1]
                    // i1: [find_left + 1, index]
                    int i1 = find_left + 1 + static_cast<int>(rand() / (((double)RAND_MAX) / (leftNum - 1)));
                    // random number [0, find_right - index - 1]
                    // i2: [index, find_right - 1]
                    int i2 = index + static_cast<int>(rand() / (((double)RAND_MAX) / (rightNum - 1)));
                    int i_real = 0;
                    for (int i = 0; i < i2 - i1 + 1; ++i)
                        if (i1 + i != index)
                            intermediate_scores[i_real++].second = static_cast<int>(all_seg[i1 + i].score);
                    qsort(intermediate_scores, i2 - i1, sizeof(intermediate_scores[0]), cmpPairInts);
                    for (int i = i2 - i1 - 1; i >= 0; --i)
                        intermediate_scores[i + 1].second = intermediate_scores[i].second;
                    intermediate_scores[0].second = static_cast<int>(all_seg[index].score);
                    double total_score = 0;
                    double fact = 1;
                    double best_pvalue1 = 1;
                    int totalLength = blocks.totalLength;
                    if (i2 < cur_seg - 1)
                        totalLength = all_seg[i2 + 1].align_start;
                    if (i1 > 0)
                        totalLength -= all_seg[i1 - 1].align_end;
                    double mean_length = log(K*totalLength*totalLength)/H;
                    for (int i = 0; i < i2 - i1 + 1; ++i) {
                        total_score += (lambda * intermediate_scores[i].second
                                        - log(K * (totalLength - mean_length)
                                              * (totalLength - mean_length)));
                        fact *= (i + 1);
                        double total_normalized_score = total_score + log(fact);
                        if (total_normalized_score > kMaxTotalScore)
                            total_normalized_score = kMaxTotalScore;
                        // If less than min_score or large number of segments
                        if ((total_normalized_score < kMinTotalScore) ||
                            (g_MULTI_PV[i + 1][static_cast<int>(total_normalized_score - kMinTotalScore)] < 0))
                            break;
                        if (g_MULTI_PV[i + 1][static_cast<int>(total_normalized_score - kMinTotalScore)] *
                            pow(2, i + 1) < best_pvalue1)
                            best_pvalue1 = (g_MULTI_PV[i + 1][static_cast<int>(total_normalized_score - kMinTotalScore)]
                                            * pow(2, i + 1));
                    }
                    if (best_pvalue1 > max_pvalue) {
                        max_pvalue = best_pvalue1;
                        all_seg[index].left_end = i1 - 1;
                        all_seg[index].right_end = i2 + 1;
                    }
                    ++contextCount;
                } // end while loop for kMaxNumberContexts
            } else {
                // if true, enumerate all context
                // for possible context containing the current segment
                for (int i1 = find_left + 1; i1 <= index; ++i1)   // start from the left of i1
                    for (int i2 = index; i2 < find_right; ++i2) { // to the right of i2
                        int i_real = 0;
                        for (int i = 0; i < i2 - i1 + 1; ++i) {
                            if (i1 + i != index) {
                                intermediate_scores[i_real++].second = static_cast<int>(all_seg[i1 + i].score);
                            }
                        }
                        qsort(intermediate_scores, i2 - i1, sizeof(intermediate_scores[0]), cmpPairInts);
                        for (int i = i2 - i1 - 1; i >= 0; --i) {
                            intermediate_scores[i + 1].second = intermediate_scores[i].second;
                        }
                        intermediate_scores[0].second = static_cast<int>(all_seg[index].score);
                        double total_score = 0;
                        double fact = 1;
                        double best_pvalue1 = 1;
                        int totalLength = blocks.totalLength;
                        if (i2 < cur_seg - 1) {
                            totalLength = all_seg[i2 + 1].align_start;
                        }
                        if (i1 > 0) {
                            totalLength -= all_seg[i1 - 1].align_end;
                        }
                        double mean_length = log(K * totalLength * totalLength) / H;
                        for (int i = 0; i < i2 - i1 + 1; ++i) {
                            total_score += (lambda * intermediate_scores[i].second
                                            - log(K * (totalLength - mean_length)
                                                  * (totalLength - mean_length)));
                            fact *= (i + 1);
                            double total_normalized_score = total_score + log(fact);
                            if (total_normalized_score > kMaxTotalScore) {
                                total_normalized_score = kMaxTotalScore;
                            }
                            // If less than min_score or large number of segments
                            if ((total_normalized_score < kMinTotalScore) ||
                                (g_MULTI_PV[i + 1][static_cast<int>(total_normalized_score - kMinTotalScore)] < 0)) {
                                break;
                            }
                            if ((g_MULTI_PV[i + 1][static_cast<int>(total_normalized_score - kMinTotalScore)]
                                 * pow(2, i + 1)) < best_pvalue1) {
                                best_pvalue1 = (g_MULTI_PV[i + 1][static_cast<int>(total_normalized_score - kMinTotalScore)]
                                                * pow(2, i + 1));
                            }
                        }
                        if (best_pvalue1 > max_pvalue) {
                            max_pvalue = best_pvalue1;
                            all_seg[index].left_end = i1 - 1;
                            all_seg[index].right_end = i2 + 1;
                        }
                    } // end of loop for possible context (i1 & i2)
            }
            all_seg[index].nearest_left = index;
            all_seg[index].nearest_right = -1;
            all_seg[index].pv2 = max_pvalue;
            if ((find_left >= 0) && (max_pvalue < all_seg[find_left].pv) && take_left) {
                max_pvalue = all_seg[find_left].pv;
                all_seg[index].nearest_left = find_left;
                all_seg[index].nearest_right = 0;
                all_seg[index].left_end = all_seg[find_left].left_end;
                all_seg[index].right_end = all_seg[find_left].right_end;
            }
            if ((find_right < cur_seg) && (max_pvalue < all_seg[find_right].pv) && take_right) {
                max_pvalue = all_seg[find_right].pv;
                all_seg[index].nearest_left = find_right;
                all_seg[index].nearest_right = 1;
                all_seg[index].left_end = all_seg[find_right].left_end;
                all_seg[index].right_end = all_seg[find_right].right_end;
            }
            all_seg[index].pv = max_pvalue;
            debug("pv %d\t%f\t%d\t%d\t%f\t%d\t%d\n", index, all_seg[index].score, all_seg[index].start,
                  all_seg[index].end, all_seg[index].pv, find_left, find_right);
            ++cur_best_index;
        } while (cur_best_index < cur_seg);
    }
    cout << "# num_segments= " << cur_seg << endl;
    cout << "# XXXXX num\tblock_num\tblock_score\tscore start\tend\tpvalue 1\t"
         << "nearest_left\tnearest_right\tleft_end\tright_end\tpvalue 2\talign_start\talign_end" << endl;
    int grandtotal = 0;
    int prev_good = 0;
    int minus_one = -1;
    for (int i = 0; i < cur_seg; ++i) {
        // Before first
        if ((i == 0) && (all_seg[i].align_start != 0)) {
            cout << "# Seg_fis " << i << "\t" << all_seg[i].block_num << "\t"
                 << all_seg[i].block_score << "\t" << all_seg[i].score << "\t" << 0 << "\t"
                 << all_seg[i].start - 1 << "\t" << "1" << "\t" << minus_one << "\t";
            if (all_seg[i].nearest_right == -1)
                cout << all_seg[i].start;
            else
                cout << all_seg[i].nearest_right;
            cout << "\t" << "-\t" << "-\t" << "-\t" << "0 \t" << all_seg[i].align_start - 1 << endl;
        }
        // Segment
        cout << "# Segment " << i << "\t" << all_seg[i].block_num << "\t" << all_seg[i].block_score
             << "\t" << all_seg[i].score << "\t" << all_seg[i].start << "\t" << all_seg[i].end
             << "\t" << all_seg[i].pv << "\t" << all_seg[i].nearest_left << "\t"
             << all_seg[i].nearest_right << "\t" << all_seg[i].left_end << "\t"
             << all_seg[i].right_end << "\t" << all_seg[i].pv2 << "\t" << all_seg[i].align_start
             << "\t" << all_seg[i].align_end << endl;
        if (all_seg[i].pv < 0.1) {
            if (all_seg[i].start - prev_good>=50)
                grandtotal += all_seg[i].start - prev_good;
            prev_good = all_seg[i].end;
        }
        // Check next
        if ((i < cur_seg - 1) && (all_seg[i].align_end + 1 < all_seg[i + 1].align_start)) {
            cout << "# Seg_nex " << i << "\t" << all_seg[i].block_num << "\t" << all_seg[i].block_score
                 << "\t" << all_seg[i].score << "\t";
            if (all_seg[i].end + 1 < all_seg[i + 1].start) {
                cout << all_seg[i].end + 1 << "\t" << all_seg[i + 1].start - 1 << "\t" << "1" << "\t";
            // When it's gaps in reference
            } else {
                cout << all_seg[i + 1].start << "\t" << all_seg[i + 1].start << "\t" << "1" << "\t";
            }
            if (all_seg[i].nearest_left == -1) {
                cout << all_seg[i].end << "\t";
            } else {
                cout << all_seg[i].nearest_left << "\t";
            }
            if (all_seg[i + 1].nearest_right == -1) {
                cout << all_seg[i + 1].start;
            } else {
                cout << all_seg[i + 1].nearest_right;
            }
            cout << "\t" << "-\t" << "-\t" << "-\t";
            cout << all_seg[i].align_end + 1 << "\t" << all_seg[i + 1].align_start - 1 << endl;
        }
        // Check end
        if ((i == cur_seg - 1) && ((all_seg[i].align_end) < (g_ALIGN_LEN - 1))) {
            cout << "# Seg_end " << i << "\t" << all_seg[i].block_num << "\t"
                 << all_seg[i].block_score << "\t" << all_seg[i].score << "\t"
                 << all_seg[i].end + 1 << "\t" << g_CHR_START + g_CHR_LEN - 1 << "\t" << "1" << "\t";
            if (all_seg[i].nearest_left == -1) {
                cout << all_seg[i].end - g_CHR_LEN << "\t" << minus_one;
            } else { 
                cout << all_seg[i].nearest_left << "\t" << minus_one;
            }
            cout << "\t" << "-\t" << "-\t" << "-\t";
            cout << all_seg[i].align_end + 1 << "\t" << g_ALIGN_LEN - 1 << endl;
        }
    }
    cout << "# grandtotal = " << grandtotal << endl;
    delete[] all_seg;
    return 0;
}
void compute_pvalue(Tree* node, blocks_type& blocks, double* branch_pvalue,
                    int& branch_num, Tree* branch_ptr[], int branchIndexToProcess, char* mafFile) {
    /* walks the tree recursively, calling compute_pvalue_branch() for each branch in the tree
     */
    if (node->is_leaf_num < 0) {
        Tree* ll = node->left_subtree;
        Tree* rr = node->right_subtree;
        if (node == root) {
            branch_ptr[branch_num] = ll;
            if (branch_num == branchIndexToProcess) {
                debug("compute_pvalue_branch() for branch index %d\n", branch_num);
                branch_pvalue[branch_num++] = compute_pvalue_branch(ll, node, blocks, mafFile);
            } else {
                branch_pvalue[branch_num++] = 0;
            }
        } else {
            branch_ptr[branch_num] = ll;
            if (branch_num == branchIndexToProcess) {
                debug("compute_pvalue_branch() for branch index %d\n", branch_num);
                branch_pvalue[branch_num++] = compute_pvalue_branch(ll, node, blocks, mafFile);
            } else {
                branch_pvalue[branch_num++] = 0;
            }
            branch_ptr[branch_num] = rr;
            if (branch_num == branchIndexToProcess) {
                debug("compute_pvalue_branch() for branch index %d\n", branch_num);
                branch_pvalue[branch_num++] = compute_pvalue_branch(rr, node, blocks, mafFile);
            } else {
                branch_pvalue[branch_num++] = 0;
            }
        }
        if (ll->is_leaf_num < 0)
            compute_pvalue(ll, blocks, branch_pvalue, branch_num,
                           branch_ptr, branchIndexToProcess, mafFile);
        if (rr->is_leaf_num < 0)
            compute_pvalue(rr, blocks, branch_pvalue, branch_num,
                           branch_ptr, branchIndexToProcess, mafFile);
    }
}
#endif // PVALUE_H_
