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

#include <iostream>
#include "tree.h"
#include "alpha.h"
#include "global.h"
#include "block.h"
#include <math.h>


Tree* root;
Tree* list_of_leaves[MAX_NUM_SPECIES];

#include "pvalue.h"

// Helper Functions
Tree* read_tree(char newick_format[],int& sp_num,int& pos);
void generate_and_read_newick_tree();
void read_multi_pv(char* filename);
void init_background();



// Main program

int main(int ARGC,char* ARGV[]) {

    BRANCH_MULTIPLIER = atof(ARGV[3]);
    char* pos = ARGV[0];
    while (strstr(pos,"SigMA")!=NULL)
        pos = strstr(pos,"SigMA")+1;
    pos--;
    char binary_path[1000];
    strncpy(binary_path,ARGV[0],pos-ARGV[0]);
    binary_path[int(pos-ARGV[0])]='\0';
    // Generates and reads a newick tree using the FASTA file in arg1
    generate_and_read_newick_tree();

#ifdef DEBUG
    cout << "# Generated Tree " << root->tree_name << " " << "#species = " << species_num << endl;
#endif

    // Read the pvalue distribution of the various scores and number of segments
    char filename[1000];
    sprintf(filename,"%smulti_segment_pvalue.txt",binary_path);
    read_multi_pv(filename);
#ifdef DEBUG
    cout << "# Read multipv " << filename << endl;
#endif

#ifdef PROTEIN
    sprintf(filename,"%s/pmb.dat",binary_path);
    read_aa_matrix(filename);
#ifdef DEBUG
    cout << "# Read AA matrix " << filename << endl;
#endif
#endif

    init_background();  

    blocks_type blocks;
    /*
      blocks.read_from_stream(ARGV[2]);
      #ifdef DEBUG
      cout << "# Read #blocks=" << blocks.chain.size() << " ; last blocklength= " << blocks.chain[blocks.chain.size()-1]->length << endl;
      #endif
    */

    double pvalue[2*MAX_NUM_SPECIES];
    Tree* branch_ptr[2*MAX_NUM_SPECIES];
    int branch=0;
    compute_pvalue(root,blocks,pvalue,branch,branch_ptr,atoi(ARGV[2]),ARGV[1]);

    double max_pvalue=-1;
    Tree* weakest_branch=NULL;
    for (int i=0;i<branch;i++) {
        if (pvalue[i]>max_pvalue) {
            max_pvalue = pvalue[i];
            weakest_branch = branch_ptr[i];
        }
    }

    cout << "p-value= " << max_pvalue << "\tbranch= " << weakest_branch->tree_name << "\ttree= " << root->tree_name << endl;

    return 0;
}






/*--------------------- Helper Functions ------------------ */


Tree* read_tree(char newick_format[],int& sp_num,int& pos) {
    while (isblank(newick_format[pos]))
        pos++;
    if (newick_format[pos]=='(') {
        pos++;
        Tree* left=read_tree(newick_format,sp_num,pos);
        while (newick_format[pos]!=':')
            pos++;
        pos++;
        while (isblank(newick_format[pos]))
            pos++;
        int start_pos=pos;
        while (newick_format[pos]!=',')
            pos++;
        char branch_length[100];
        strncpy(branch_length,newick_format+start_pos,pos-start_pos);
        branch_length[pos-start_pos]='\0';
        float left_branch=atof(branch_length);
        if (left_branch<0.001) left_branch=0.001;
        left_branch *= BRANCH_MULTIPLIER;
        if (left_branch>MAX_BRANCH_LENGTH) left_branch=MAX_BRANCH_LENGTH;
    
        pos++;
        Tree* right=read_tree(newick_format,sp_num,pos);
        while (newick_format[pos]!=':')
            pos++;
        pos++;
        while (isblank(newick_format[pos]))
            pos++;
        start_pos=pos;
        while (newick_format[pos]!=')')
            pos++;
        strncpy(branch_length,newick_format+start_pos,pos-start_pos);
        branch_length[pos-start_pos]='\0';
        float right_branch=atof(branch_length);
        if (right_branch<0.001) right_branch=0.001;
        right_branch *= BRANCH_MULTIPLIER;
        if (right_branch>MAX_BRANCH_LENGTH) right_branch=MAX_BRANCH_LENGTH; 
    
        pos++;
        return (new Tree(left,right,left_branch,right_branch));
    }
    else {
        int start_pos=pos;
        while (newick_format[pos]!=':')
            pos++;
        strncpy(species_name[sp_num],newick_format+start_pos,pos-start_pos);
        species_name[sp_num][pos-start_pos]='\0';    
        Tree* node=new Tree(sp_num,species_name[sp_num]);
        list_of_leaves[sp_num]=node;
        sp_num++;
        return node;
    }
}


void generate_and_read_newick_tree() {

    int pos=0;
    species_num=0;
    root=read_tree(PHYLOGENY,species_num,pos);
}


void read_multi_pv(char* filename) {
    multi_pv = new double*[MAX_SEGMENTS+1];
    for (int i=0;i<MAX_SEGMENTS+1;i++)
        multi_pv[i] = new double[int(MAX_TOTAL_SCORE-MIN_TOTAL_SCORE+1)];

    for (int i=0;i<MAX_SEGMENTS+1;i++)
        for (int j=0;j<MAX_TOTAL_SCORE-MIN_TOTAL_SCORE+1;j++)
            multi_pv[i][j]=-1;
  
    for (int j=MIN_TOTAL_SCORE;j<=MAX_TOTAL_SCORE;j++)
        multi_pv[1][j-MIN_TOTAL_SCORE]=1-exp(-1*exp(-1*j));
    ifstream ifs(filename);
    int a,b;
    double c;
    if (!ifs) {
        cout << "multi_segment_pvalue.txt not found ... quitting\n";
        exit(0);
    }
    while (ifs) {
        ifs >> a >> b >> c;
        if ((a<=MAX_SEGMENTS)&&(b<=MAX_TOTAL_SCORE)&&(b>=MIN_TOTAL_SCORE)) {
            //      if (b<0.5) c=1; 
            multi_pv[a][b-MIN_TOTAL_SCORE]=c;
        }
    }
    ifs.close();
}

void init_background() {
    back[0]=0;
#ifdef PROTEIN
    for (int i=1;i<ALPHABET_SIZE;i++)
        back[i]=aa_back[i-1]*(1-GAP_BACKGRD);
#else
    for (int i=1;i<ALPHABET_SIZE;i++)
        back[i]=(1-GAP_BACKGRD-back[0])/float(ALPHABET_SIZE-1);
#endif
    back[GAP]=GAP_BACKGRD;
}
