/*VI. Conditions of Use
 *
 * Xiaoyu Chen, Amol Prakash, and Martin Tompa give permission for Dent Earl and his institution to use the StatSigmMA - w software developed at the University of Washington for research purposes, on the following conditions:
 *
 * 1. Dent Earl and his institutional colleagues may modify the StatSigmMA - w software and distribute the resulting modified software for research purposes, provided (a) that it is distributed  together with these Conditions of Use, (b) that Martin Tompa receive a copy of the finalized modified software, and (c) that Xiaoyu Chen, Amol Prakash, and Martin Tompa are credited with the authorship of the software.
 *  2. The StatSigmMA - w software will be used by you and / or your institution solely for noncommercial purposes, except with express permission from the authors.
 * 3. Any risk associated with using the StatSigmMA - w software at your institution is with you and your institution.
 * 4. StatSigmMA - w will be cited in any publication(s) reporting on data obtained from it as:
 *    Amol Prakash and Martin Tompa, "Measuring the Accuracy of Genome - Size Multiple Alignments". Genome Biology, vol. 8, issue 6, June 2007, R124.
 *    Xiaoyu Chen and Martin Tompa, "Comparative assessment of methods for aligning multiple genome sequences". Nature Biotechnology, vol. 28, no. 6, June 2010, 567 - 572.
 */
#ifndef TREE_H_
#define TREE_H_

#include <set>
#include "alpha.h"
#include <fstream>
#include <iostream>
#include <string.h>

using namespace std;
using namespace __gnu_cxx;

class Tree {
 public:
    double left_branchlength;
    double right_branchlength;
    Tree* left_subtree;
    Tree* right_subtree;
    Tree* parent;
    double prob[c_ALPHABET_SIZE + 1];
    int is_leaf_num;
    int leaf_startnum;
    int leaf_endnum;
    int nucleotide;
    bool non_ortho;
    char tree_name[1000];
    Tree(int sp_num, char* sp_name);
    Tree(Tree* l, Tree* r, double l_br, double r_br);
    ~Tree();
    void compute_prob();
    double return_prob();
    int max_prob();
    bool all_equal_prob();
    bool all_random(Tree* t);
};
bool Tree::all_random(Tree* t) {
    if (this == t) return true;
    else if (is_leaf_num >= 0)
        return (nucleotide == c_RANDOM);
    return (left_subtree->all_random(t) && right_subtree->all_random(t));
}
bool Tree::all_equal_prob() {
    int a1 = int(prob[1] * 1e5);
    int a2 = int(prob[2] * 1e5);
    int a3 = int(prob[3] * 1e5);
    int a4 = int(prob[4] * 1e5);
    int a5 = int(prob[5] * 1e5);
    if ((a1 == a2) && (a1 == a3) && (a1 == a4) && (a1!=0) && ((a1 == a5) || ((a5 == 1e5) && (a1 > 9900))))
        return true;
    else
        return false;
}
int Tree::max_prob() {
    int max_val = -1;
    double max_prob = 0;
    for (int i = 1; i < c_ALPHABET_SIZE; i++)
        if (prob[i] >= max_prob) {
            max_prob = prob[i];
            max_val = i;
        }
    // This subtree is a parent
    if (all_equal_prob() || ((is_leaf_num < 0) && 
                             (left_subtree->all_equal_prob() || right_subtree->all_equal_prob())))
        {}
    else if (max_prob * 0.9 < prob[c_GAP])
        max_val = c_GAP;
    return max_val;
}
Tree::Tree(int a, char* sp_name) {
    left_subtree = NULL;
    right_subtree = NULL;
    is_leaf_num = a;
    strcpy(tree_name, sp_name);
    leaf_startnum = leaf_endnum = a;
    non_ortho = false;
    parent = NULL;
}
Tree::Tree(Tree* l, Tree* r, double l_br, double r_br) {
    left_subtree = l;
    l->parent = this;
    right_subtree = r;
    r->parent = this;
    left_branchlength = l_br;
    right_branchlength = r_br;
    is_leaf_num = -1;
    leaf_startnum = l->leaf_startnum;
    leaf_endnum = r->leaf_endnum;
    non_ortho = false;
    sprintf(tree_name, "(%s %s)", l->tree_name, r->tree_name);
    parent = NULL;
}
Tree::~Tree() {
    if (left_subtree) delete left_subtree;
    if (right_subtree) delete right_subtree;
}
void Tree::compute_prob() {
    if (is_leaf_num >= 0) {
        for (int i = 1; i <= c_ALPHABET_SIZE; i++)
            prob[i] = 0;
        if (nucleotide == c_RANDOM)
            for (int i = 1; i <= c_ALPHABET_SIZE; i++)
                prob[i] = 1;
        else
            prob[nucleotide] = 1;
    } else {
        left_subtree->compute_prob();
        right_subtree->compute_prob();    
        for (int i = 1; i <= c_ALPHABET_SIZE; i++) {
            double left_prod = 0;
            double right_prod = 0;
            if (!left_subtree->non_ortho)
                for (int j = 1; j <= c_ALPHABET_SIZE; j++)
                    left_prod += left_subtree->prob[j] * evolution_prob(i, j, left_branchlength);
            else
                for (int j = 1; j <= c_ALPHABET_SIZE; j++)
                    left_prod += left_subtree->prob[j] * g_back[j];
            if (!right_subtree->non_ortho)
                for (int j = 1; j <= c_ALPHABET_SIZE; j++)
                    right_prod += right_subtree->prob[j] * evolution_prob(i, j, right_branchlength);
            else
                for (int j = 1; j <= c_ALPHABET_SIZE; j++)
                    right_prod += right_subtree->prob[j] * g_back[j];
            prob[i] = left_prod * right_prod;
            // Make sure this is right
            //      if ((left_subtree->nucleotides.find(NON_ORTHO)!=left_subtree->nucleotides.end()) && (right_subtree->nucleotides.find(NON_ORTHO)!=right_subtree->nucleotides.end()))
            //	nucleotides.insert(NON_ORTHO);
        }
    }
    //  cout << is_leaf_num << "\t" << prob[1] << "\t" << prob[2] << "\t" << prob[3] << "\t" << prob[4] << "\t" << prob[5] << endl;
}

double Tree::return_prob() {
    //  cout << tree_name << endl;
    //  cout << "returning " << is_leaf_num << "\t" << prob[1] << "\t" << prob[2] << "\t" << prob[3] << "\t" << prob[4] << "\t" << prob[5] << endl;
    if (is_leaf_num<0) {
        if (left_subtree->all_equal_prob() || 
            right_subtree->all_equal_prob() || 
            left_subtree->non_ortho || 
            right_subtree->non_ortho)
            return (right_subtree->return_prob())*(left_subtree->return_prob());
    }
    double ff = 0;
    for (int i = 1; i <= c_ALPHABET_SIZE; i++)
        ff += prob[i] * g_back[i];
    //  cout << "returned " << is_leaf_num << "\t" << prob[1] << "\t" << prob[2] << "\t" << prob[3] << "\t" << prob[4] << "\t" << prob[5] << "\t" << ff << endl;
    return ff;
}

#endif // TREE_H_
