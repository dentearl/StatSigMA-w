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
#ifndef PVALUE_H
#define PVALUE_H

#include "tree.h"
#include "block.h"
#include "blast.h"
//#include <pair.h>
#include <stdlib.h>
#include <math.h>

using namespace std;
using namespace __gnu_cxx;

int compare(const void* a, const void* b) {
    double aa = ((pair<int,int>*) a)->second;
    double bb = ((pair<int,int>*) b)->second;
    if (aa > bb) return -1;
    if (aa < bb) return 1;
    else return 0;
}

bool compute_branch_parameters(Tree* node, Tree* parent, const blocks_type& blocks, 
                               double& K, double& lambda, double& H) {
  
#ifdef DEBUG
    cout << "# Computing branch parameters for node number " << node->leaf_startnum << " " << node->leaf_endnum << endl;
#endif

    srand(time(NULL));

    int max_score = 10 * SCORE_BINS;
    int min_score = -10 * SCORE_BINS;
    double *score_dist = (double*) malloc(max_score - min_score + 1);
    for (int i = 0; i < max_score - min_score + 1; i++)
        score_dist[i] = 0;
  
    double num_tuples = pow(ALPHABET_SIZE, species_num);
    int total_computations = 0;
#ifdef DEBUG
    cout << "# Comparison of number of tuples all=" << num_tuples << " " << "block_len= " << blocks.chain[0]->length << endl;
#endif
    if (num_tuples < blocks.chain[0]->length * blocks.chain[0]->length) {
        total_computations=10000; // default value
#ifdef DEBUG
        cout << "# Computing all scores " << num_tuples << "\n";
#endif
    
        for (int num=0;num<num_tuples;num++) {
            for (int k=0;k<species_num;k++)
                list_of_leaves[k]->nucleotide = (int(num/pow(ALPHABET_SIZE,species_num-k-1)))%(ALPHABET_SIZE) +1;
      
            root->compute_prob();
            double ortho_pv=root->return_prob();
            double non_ortho_pv;
            if (parent->max_prob()==GAP)
                non_ortho_pv = ortho_pv;
            else {
                node->non_ortho = true;
                root->compute_prob();
                non_ortho_pv=root->return_prob();
                node->non_ortho = false;
            }
            int score=((ortho_pv==0)||(non_ortho_pv==0)) ? 0 : int(round(log(ortho_pv/non_ortho_pv)*SCORE_BINS));
      
            if (score<min_score) score=min_score;
            if (score>max_score) score=max_score;
            score_dist[score-min_score] += non_ortho_pv;
        }
    }
  
    else {
        /*
          for (int i=0;i<max_score-min_score+1;i++)
          score_dist[i] += double(PSEUDOCOUNT_NORMALISATION)/(blocks.chain[0]->length*(max_score-min_score));
        */
#ifdef DEBUG
        cout << "# Computing sequence specific scores max(N^2," << TOTAL_NUM_TUPLES << ") where N=" << blocks.chain[0]->length << " (again based on random sampling)\n";
#endif
    
        int bl1=0;
        int bl2=0;
    
        int* good_tuple1 = (int *)malloc(sizeof(int) * MAX_BLOCK_SIZE);
        int* good_tuple2 = (int *)malloc(sizeof(int) * MAX_BLOCK_SIZE);
        int good1 = 0, good2 = 0;
        for (int i=0;i<MAX_BLOCK_SIZE;i++)
            good_tuple1[i]=good_tuple2[i]=-1;

        for (int l1=0;l1<blocks.chain[bl1]->length;l1++) {
            for (int i=0;i<node->leaf_startnum;i++)
                if (blocks.chain[bl1]->block[l1][i]!=RANDOM)
                    good_tuple1[good1] = l1;
            for (int i=node->leaf_endnum+1;i<species_num;i++)
                if (blocks.chain[bl1]->block[l1][i]!=RANDOM)
                    good_tuple1[good1] = l1;
            if (good_tuple1[good1] == l1) good1++;
      
        }
        for (int l2=0;l2<blocks.chain[bl2]->length;l2++) {
            for (int i=node->leaf_startnum;i<=node->leaf_endnum;i++)
                if (blocks.chain[bl2]->block[l2][i]!=RANDOM)
                    good_tuple2[good2] = l2;
            if (good_tuple2[good2] == l2) good2++;
        }

#ifdef DEBUG
        cout<<"Good1 "<<good1<<endl;
        cout<<"Good2 "<<good2<<endl;
#endif

        if ((good1 > 0) && (good2 > 0)) {      
            for (int tp=0; tp<TOTAL_NUM_TUPLES; tp++) {
                // random choose a tuple from good_tuple1 and good_tuple2
                int goodID1 = (int)(rand() / (((double)RAND_MAX+1) / good1));
                int goodID2 = (int)(rand() / (((double)RAND_MAX+1) / good2));
                int l1 = good_tuple1[goodID1];
                int l2 = good_tuple2[goodID2];

#ifdef DEBUG
                //if (tp < 100)
                //  cout << l1 << "\t" << l2 << endl;
                if (total_computations%10000==0)
                    cout << "# Tuples read = " << total_computations+1 << endl;
#endif
                total_computations++;
      
                for (int i=0;i<node->leaf_startnum;i++)
                    list_of_leaves[i]->nucleotide=blocks.chain[bl1]->block[l1][i];
                for (int i=node->leaf_startnum;i<=node->leaf_endnum;i++)
                    list_of_leaves[i]->nucleotide=blocks.chain[bl2]->block[l2][i];
                for (int i=node->leaf_endnum+1;i<species_num;i++)
                    list_of_leaves[i]->nucleotide=blocks.chain[bl1]->block[l1][i];
      
                root->compute_prob();
                double ortho_pv=root->return_prob();
      
                double non_ortho_pv;
                if (parent->max_prob()==GAP)
                    non_ortho_pv = ortho_pv;
                else {
                    node->non_ortho = true;
                    root->compute_prob();
                    non_ortho_pv=root->return_prob();
                    node->non_ortho = false;
                }
                int score=((ortho_pv==0)||(non_ortho_pv==0)) ? 0 : int(log(ortho_pv/non_ortho_pv)*SCORE_BINS);
                // cout << "\t" << score << "  " << ortho_pv << " " << non_ortho_pv << endl;
                if (score<min_score) score=min_score;
                if (score>max_score) score=max_score;
                score_dist[score-min_score]++;      
            }
        }
        else {
#ifdef DEBUG
            cout<<"No valid tuples: Good1="<<good1<<" Good2="<<good2<<endl;
#endif
        }
    
        free(good_tuple1);
        free(good_tuple2);

#ifdef DEBUG
        cout<<"Score sampling finished"<<endl;
        cout<<timeStamp()<<endl;
#endif

        /*
          bool good_tuple1[MAX_BLOCK_SIZE];
          bool good_tuple2[MAX_BLOCK_SIZE];
          for (int i=0;i<MAX_BLOCK_SIZE;i++)
          good_tuple1[i]=good_tuple2[i]=false;
          int good1=0;
          int good2=0;
          for (int l1=0;l1<blocks.chain[bl1]->length;l1++) {
          for (int i=0;i<node->leaf_startnum;i++)
          if (blocks.chain[bl1]->block[l1][i]!=RANDOM)
          good_tuple1[l1]=true;
          for (int i=node->leaf_endnum+1;i<species_num;i++)
          if (blocks.chain[bl1]->block[l1][i]!=RANDOM)
          good_tuple1[l1]=true;
          if (good_tuple1[l1]) good1++;
          }
          for (int l2=0;l2<blocks.chain[bl2]->length;l2++) {
          for (int i=node->leaf_startnum;i<=node->leaf_endnum;i++)
          if (blocks.chain[bl2]->block[l2][i]!=RANDOM)
          good_tuple2[l2]=true;
          if (good_tuple2[l2]) good2++;
          }

          double pval=0;
          if ((good1>0) &&(good2>0))
          pval=float(TOTAL_NUM_TUPLES)/(float(good1) * good2);

          for (int l1=0;l1<blocks.chain[bl1]->length;l1++) {
          if (good_tuple1[l1])
          for (int l2=0;l2<blocks.chain[bl2]->length;l2++)
          if (good_tuple2[l2])
          if (pval >= float(rand())/RAND_MAX) {
          #ifdef DEBUG
	      if (total_computations%10000==0)
          cout << "# Tuples read = " << total_computations+1 << endl;
          #endif
	      total_computations++;
	      
	      for (int i=0;i<node->leaf_startnum;i++)
          list_of_leaves[i]->nucleotide=blocks.chain[bl1]->block[l1][i];
	      for (int i=node->leaf_startnum;i<=node->leaf_endnum;i++)
          list_of_leaves[i]->nucleotide=blocks.chain[bl2]->block[l2][i];
	      for (int i=node->leaf_endnum+1;i<species_num;i++)
          list_of_leaves[i]->nucleotide=blocks.chain[bl1]->block[l1][i];
	      
	      root->compute_prob();
	      double ortho_pv=root->return_prob();
	      
	      double non_ortho_pv;
	      if (parent->max_prob()==GAP)
          non_ortho_pv = ortho_pv;
	      else {
          node->non_ortho = true;
          root->compute_prob();
          non_ortho_pv=root->return_prob();
          node->non_ortho = false;
	      }
	      int score=((ortho_pv==0)||(non_ortho_pv==0)) ? 0 : int(log(ortho_pv/non_ortho_pv)*SCORE_BINS);
	      //	  cout << "\t" << score << "  " << ortho_pv << " " << non_ortho_pv << endl;
	      if (score<min_score) score=min_score;
	      if (score>max_score) score=max_score;
	      score_dist[score-min_score]++;
	      
	      if (total_computations>TOTAL_NUM_TUPLES) {
          l2 = blocks.chain[bl2]->length;
          l1 = blocks.chain[bl1]->length;
	      }
          }
          }
        */
    }

    /*
      #ifdef DEBUG
      for (int i=0;i<max_score-min_score+1;i++)
      if (score_dist[i]!=0)
      cout << "# Score_dist " << i << " " << score_dist[i] << endl;
      #endif
    */
    int new_max=0;
    int new_min=0;
    for (int i=0;i<max_score-min_score+1;i++)
        if (score_dist[i]!=0) {
            if (i+min_score>new_max)
                new_max=i+min_score;
            if (i+min_score<new_min)
                new_min=i+min_score;
        }
    for (int i=0;i<new_max-new_min+1;i++)
        score_dist[i]=score_dist[i-min_score+new_min];
  

#ifdef DEBUG
    cout<<"Finish score distribution"<<endl;
    cout<<timeStamp()<<endl;
#endif

    if (total_computations<5) {
        K=1;
        lambda=0;
        H=1;
#ifdef DEBUG
        cout << "# Too few good tuples\n";
#endif
        return false;
    }
    //  else if (! karlin(min_score,max_score,score_dist,&lambda,&K,&H))
    else if (! karlin(new_min,new_max,score_dist,&lambda,&K,&H)) {
        cout << "Karlin Altschul estimation Failed\n";

#ifdef DEBUG
        cout<<"Finished karlin"<<endl;
        cout<<timeStamp()<<endl;
#endif

        return false;
        exit(0);
    }
    //  lambda *= SCORE_BINS;          b'coz all scores are normalized

#ifdef DEBUG
    cout << "# Parameters K=" << K << "  lambda=" << lambda << " H=" << H << endl;
    cout<<timeStamp()<<endl;
#endif
    return true;
}



struct final_seg_type {
    double score;       // segment score
    int start;
    int end;
    int align_start;
    int align_end;
    double pv;
    double pv2;
    int block_num;       // Block number from reading maf
    int nearest_left;    // position of left nearest significant segment
    int nearest_right;   // position of right nearest significant segment
    int block_score;     //  list<double> block_score;
    int left_end;
    int right_end;
    final_seg_type() {pv=-1;}
};



void scan_alignFile(char* mafFile) {

    CHR_START = -1;
    CHR_LEN = 0;
    ALIGN_LEN = 0;
    int tempSize = 10 * MAX_BLOCK_SIZE;

    int last_pos = 0;
    int cur_len = 0;
    char *buffer = (char *) malloc(tempSize);
    bool firstBlock = true;

    ifstream ifs(mafFile);
    while (ifs) {
        buffer[0]='\0';
        ifs.getline(buffer,tempSize);
        if ((strncmp(buffer,"a score=",8)==0) || (!ifs)) {
            // make sure MAX_BLOCK_SIZE >= actual alignment block size
            if (cur_len > MAX_BLOCK_SIZE)
                MAX_BLOCK_SIZE = cur_len;

            ALIGN_LEN += cur_len;
            cur_len = 0;

            // sequence length = last position - first position + 1
            if (last_pos > 0)
                CHR_LEN = last_pos - CHR_START + 1;
        }

        else {
            char* a;
            if (buffer[0]=='s') {
                // Checking for human                                                   
                bool is_human=false;
                if (strncmp(buffer+2,REF_SPECIES,strlen(REF_SPECIES))==0)
                    is_human=true;

                if (is_human) {
                    a = strtok(buffer," ");
                    a = strtok(NULL," ");
                    a = strtok(NULL," ");
                    if (firstBlock) {
                        CHR_START = atoi(a);
                        firstBlock = false;
                    }
                    last_pos = atoi(a);

                    a = strtok(NULL," ");
                    // get the last genomic position
                    last_pos += atoi(a)-1;

                    a = strtok(NULL," ");
                    a = strtok(NULL," ");
                    a = strtok(NULL," ");
                    cur_len = strlen(a);
	  
                }
            }
        }
    }
    ifs.close();
}



double compute_pvalue_branch(Tree* node,Tree* parent,blocks_type& blocks,char* filename) {

    int max_seg = 6000000;
    final_seg_type* all_seg = new final_seg_type[max_seg+1];

    scan_alignFile(filename);
#ifdef DEBUG
    cout<<"Chrmosome start position = "<<CHR_START<<endl;
    cout<<"Sequence(non-gap) length = "<<CHR_LEN<<endl;
    cout<<"Alignment length = "<<ALIGN_LEN<<endl;
    cout<<"Max block size = "<<MAX_BLOCK_SIZE<<endl;
#endif


    double K=0.29;
    double lambda=0.0011;
    double H=0.68;

    //   compute_branch_parameters(node,parent,blocks,K,lambda,H);
    srand(time(NULL));
    double K_arr[TOTAL_ITERATE_PARAM];
    double lambda_arr[TOTAL_ITERATE_PARAM];
    double H_arr[TOTAL_ITERATE_PARAM];
    blocks_type* blocks_param = new blocks_type;
    ifstream ifs1(filename);
    // int runs = 0;
  
    float rate = ALIGN_LEN / (float)TOTAL_ITERATE_PARAM;

#ifdef DEBUG
    cout<<"align length "<<ALIGN_LEN<<endl;
#endif

    int sim_count = 0;
    int block_count = 0;
    while ((blocks_param->read_single_from_maf(ifs1))&&(sim_count<TOTAL_ITERATE_PARAM)) {
        int currBlockSize = blocks_param->chain[0]->length;
        int currRepNum = (int)(currBlockSize / rate);      
        // if failures, need try more times
        int maxTry = currRepNum * 2;
        int blockSimCount = 0;
        for (int r=0; r<maxTry; r++) {	
            if (compute_branch_parameters(node,parent,*blocks_param,K,lambda,H)) {
                K_arr[sim_count] = K;
                lambda_arr[sim_count] = lambda;
                H_arr[sim_count] = H;
                sim_count++;
                blockSimCount++;
                // get enough estimation from current block
                if (blockSimCount == currRepNum) break;
            }
#ifdef DEBUG
            else
                cout << "# That block didn't work, try again\n";
#endif
        }
      
        block_count++;
    }
    delete blocks_param;
    ifs1.close();
    cout<<"# number of sigMA blocks = "<<block_count<<endl;
  
    for (int i=0;i<sim_count;i++) {
        K_arr[0] += K_arr[i];
        lambda_arr[0] += lambda_arr[i];
        H_arr[0] += H_arr[i];
    }
    K = K_arr[0]/sim_count;
    lambda = lambda_arr[0]/sim_count;
    H = H_arr[0]/sim_count;
  
    cout << "# number of sets of parameter computed= " << sim_count << endl;


    int cur_seg=0;

    ifstream ifs(filename);
    pair<int,int> limit[MAX_SEGMENTS];
    pair<double,double> score_left_right[MAX_SEGMENTS];
    int count=0;
  
    while (blocks.read_single_from_maf(ifs)&&(cur_seg<max_seg)) {
    
        int all_good_segments=0;
        int block_good_segments=0;
        int bl=0;
    
        //    cout << "New block " << blocks.chain[bl]->start << " " << blocks.chain[bl]->length << endl;

        //    if (blocks.chain[bl]->start > 64800000) {

        double cur_score=0;
        for (int l=0;l<blocks.chain[bl]->length;l++) {
            for (int i=0;i<species_num;i++)
                list_of_leaves[i]->nucleotide = blocks.chain[bl]->block[l][i];

            root->compute_prob();
            double ortho_pv=root->return_prob();
            double non_ortho_pv;
            int score;
            if (parent->max_prob()==GAP) {
                non_ortho_pv = ortho_pv;
                score=0;
            }
            else {
                node->non_ortho = true;
                root->compute_prob();
                non_ortho_pv=root->return_prob();
                node->non_ortho = false;
                int max_score=10*SCORE_BINS;
                int min_score=-10*SCORE_BINS;
                score=((ortho_pv==0)||(non_ortho_pv==0)) ? 0 : int(round(log(ortho_pv/non_ortho_pv)*SCORE_BINS));
                if (score<min_score) score=min_score;
                if (score>max_score) score=max_score;
                // Adjusting for cutting branch for which there is no data, and approximate errors in ortho_pv/non_ortho_pv which may not be exactly zero
                //	cout << " I am here2\n";

                if (root->all_random(node) || node->all_random(NULL))
                    //	if ((score<=0)&&(score>=-50))
                    score=10000*min_score;
            }

            /*
              if (blocks.chain[bl]->start < 2000000) {
              int chars=0;
              for (int jj=0;jj<l;jj++) {
              int p = blocks.chain[bl]->block[jj][blocks.chain[bl]->human_ref];
              if ((p>=1)&&(p<=4))
              chars++;
              }
              cout << "# " << blocks.chain[bl]->start+chars << " " << score << endl;
              }
              cout << "pos="<<blocks.chain[bl]->start+chars<<"\t";
              for (int i=0;i<species_num;i++)
              cout << list_of_leaves[i]->nucleotide;
              cout << "\tscore " << score << "\t" << ortho_pv << " " << non_ortho_pv << endl;
            */

            cur_score += score;
            if (score>0) {
                int max_j=-1;
                for (int i=block_good_segments-1;i>=0;i--)
                    if ((score_left_right[all_good_segments+i].first < cur_score-score)&&(i>max_j))
                        max_j=i;
                if ((max_j==-1)||(score_left_right[all_good_segments+max_j].second >= cur_score)) {
                    score_left_right[all_good_segments+block_good_segments].first = cur_score-score;
                    score_left_right[all_good_segments+block_good_segments].second = cur_score;
                    limit[all_good_segments+block_good_segments].first = l;
                    limit[all_good_segments+block_good_segments].second = l;
                    if (block_good_segments+all_good_segments<MAX_SEGMENTS-1) {
                        block_good_segments++;
                        /*
                          #ifdef DEBUG

                          cout << "# New segment " << bl << " " << l << " " << score << endl;
                          #endif
                        */
                    }
                }
                else {
                    score_left_right[all_good_segments+max_j].second = cur_score;
                    limit[all_good_segments+max_j].second = l;
                    block_good_segments=max_j+1;
	  
                    int new_max_j=-1;
                    bool done=false;
                    do {
                        for (int i=0;i<block_good_segments-1;i++)
                            if ((score_left_right[all_good_segments+i].first<score_left_right[all_good_segments+block_good_segments-1].first)&&(i>new_max_j))
                                new_max_j=i;
                        if ((new_max_j==-1)||(score_left_right[all_good_segments+new_max_j].second >= score_left_right[all_good_segments+block_good_segments-1].second))
                            done=true;
                        else {
                            score_left_right[all_good_segments+new_max_j].second = score_left_right[all_good_segments+block_good_segments-1].second;
                            limit[all_good_segments+new_max_j].second = limit[all_good_segments+block_good_segments-1].second;
                            block_good_segments=new_max_j+1;
                        }
                    } while (!done);
	  
                }
            }
        }
    
        //    all_good_segments += block_good_segments;
        //    block_good_segments=0;
    
        for (int i=0;i<block_good_segments;i++)
            if (cur_seg<max_seg) {
                all_seg[cur_seg].score = score_left_right[i].second - score_left_right[i].first;
                all_seg[cur_seg].block_num = count;
                all_seg[cur_seg].block_score = blocks.chain[bl]->block_score.size();
                int chars=0;
                for (int j=0;j<limit[i].first;j++) {
                    int p = blocks.chain[bl]->block[j][blocks.chain[bl]->human_ref];
                    if ((p>=1)&&(p<=4))
                        chars++;
                }
                all_seg[cur_seg].start = blocks.chain[bl]->start + chars;
                all_seg[cur_seg].align_start = blocks.chain[bl]->align_start + limit[i].first;
                for (int j=limit[i].first;j<limit[i].second;j++) {
                    int p = blocks.chain[bl]->block[j][blocks.chain[bl]->human_ref];
                    if ((p>=1)&&(p<=4))
                        chars++;
                }
                all_seg[cur_seg].end = blocks.chain[bl]->start + chars;
                all_seg[cur_seg].align_end = blocks.chain[bl]->align_start + limit[i].second;

                //	cout << "segment " << all_seg[cur_seg].score << "\t" << all_seg[cur_seg].start << "\t" << all_seg[cur_seg].end << endl;
                cur_seg++;
                //	cout << blocks.chain[bl]->start << "\t" << chars << "\t" << limit[i].first << "\t" << limit[i].second << endl;
            }
        count+= blocks.chain[bl]->block_score.size();
        //    }
    }
    ifs.close();

    //  cout <<"# Me here2\n";

#ifdef DEBUG
    cout << "# num_segments=" << cur_seg << endl;
    for (int i=0;i<cur_seg;i++)
        cout << "# Segment " << i << "\t" << all_seg[i].score << "\t" << all_seg[i].start << "\t" << all_seg[i].end << "\t" << all_seg[i].align_start << "\t" << all_seg[i].align_end << endl;
#endif


    pair<int,int>* best_scores = new pair<int,int>[cur_seg];
    for (int i=0;i<cur_seg;i++) {
        best_scores[i].first = i;
        best_scores[i].second = (int)all_seg[i].score;
    }

    qsort(best_scores,cur_seg,sizeof(pair<int,int>),compare);

  
    double best_pvalue;
    int cur_best_index=0;
    if (cur_best_index<cur_seg) {
        do {
            int index=best_scores[cur_best_index].first;
      
            int total_length = blocks.total_length;
            double mean_length = log(K*total_length*total_length)/H;
            double normalized_score = lambda*all_seg[index].score-log(K*(total_length-mean_length)*(total_length-mean_length));
      
            if (normalized_score<MIN_TOTAL_SCORE)
                break;
            if (normalized_score > MAX_TOTAL_SCORE) normalized_score=MAX_TOTAL_SCORE;
            best_pvalue = multi_pv[1][int(normalized_score-MIN_TOTAL_SCORE)];
            all_seg[index].pv = best_pvalue;
      
            all_seg[index].nearest_left = index;
            all_seg[index].nearest_right = -1;
      
            all_seg[index].left_end = 0;
            all_seg[index].right_end = cur_seg;
            all_seg[index].pv2 = all_seg[index].pv;

#ifdef DEBUG      
            cout << "pv " << index << "\t" << all_seg[index].score << "\t" << all_seg[index].start << "\t" << all_seg[index].end << "\t" << all_seg[index].pv << endl;
#endif

            cur_best_index++;
        } while ((best_pvalue<PV_THRESH_DVDCQ)&&(cur_best_index<cur_seg));
    }
    // End coordinates

#ifdef DEBUG
    cout <<"# Me here3\n";
#endif

    if (cur_best_index<cur_seg) {
        do {
            int index=best_scores[cur_best_index].first;
      
            int find_left=index;
            while ((find_left>=0)&&(all_seg[find_left].pv==-1))
                find_left--;

            int find_right=index;
            while ((find_right<cur_seg)&&(all_seg[find_right].pv==-1))
                find_right++;

            int total_length=blocks.total_length;
            if (find_right<cur_seg)
                total_length = all_seg[find_right].align_start;
            if (find_left>=0)
                total_length -= all_seg[find_left].align_end;
      
            double mean_length = log(K*total_length*total_length)/H;
            double normalized_score = lambda*all_seg[index].score-log(K*(total_length-mean_length)*(total_length-mean_length));
            if (normalized_score<MIN_TOTAL_SCORE)
                break;
            if (normalized_score > MAX_TOTAL_SCORE) normalized_score=MAX_TOTAL_SCORE;
            best_pvalue = multi_pv[1][int(normalized_score-MIN_TOTAL_SCORE)];
      
            if ((find_left>=0)&&(best_pvalue < all_seg[find_left].pv))
                best_pvalue = all_seg[find_left].pv;
            if ((find_right<cur_seg)&&(best_pvalue < all_seg[find_right].pv))
                best_pvalue = all_seg[find_right].pv;
            all_seg[index].pv = best_pvalue;
      
            all_seg[index].nearest_left = index;
            all_seg[index].nearest_right = -1;
      
            all_seg[index].left_end = 0;
            all_seg[index].right_end = cur_seg;
            all_seg[index].pv2 = all_seg[index].pv;

#ifdef DEBUG
            cout << "pv " << index << "\t" << all_seg[index].score << "\t" << all_seg[index].start << "\t" << all_seg[index].end << "\t" << all_seg[index].pv << endl;
#endif

            cur_best_index++;      
        } while ((best_pvalue<PV_THRESH_DVDCQ)&&(cur_best_index<cur_seg));

#ifdef DEBUG
        cout<<"out while due to best_pvalue = "<<best_pvalue<<endl;
#endif
    }

#ifdef DEBUG
    cout <<"# Me here4"<<endl;  
#endif

    pair<int,int>* intermediate_scores = new pair<int,int>[cur_seg];

    if (cur_best_index < cur_seg) {
        do {
#ifdef DEBUG
            cout <<"cur_seg:"<<cur_seg<<" cur_best_index:"<<cur_best_index<<endl;
            cout<<timeStamp()<<endl;
#endif

            int index=best_scores[cur_best_index].first;
      
            int find_left=index;
            while ((find_left>=0)&&(all_seg[find_left].pv==-1))
                find_left--;

#ifdef DEBUG
            cout << "index " << index << endl;
            cout <<"find left: "<<find_left<<endl;
#endif
      
            int find_right=index;
            while ((find_right<cur_seg)&&(all_seg[find_right].pv==-1))
                find_right++;

#ifdef DEBUG
            cout <<"find right: "<<find_right<<endl;
#endif
            
            // find_left is the leftmost that is not included
            // find_right is the rightmost that is not included
      
            bool take_left=true;
            bool take_right=true;
            if ((find_left>=0)&&(all_seg[find_left].pv > 0.01)&&(all_seg[find_left].right_end <= index)) {
                find_left=all_seg[find_left].left_end;
                take_left=false;
            }
      
            if ((find_right<cur_seg)&&(all_seg[find_right].pv > 0.01)&&(all_seg[find_right].left_end >= index)) {
                find_right=all_seg[find_right].right_end;
                take_right=false;
            }
      
#ifdef DEBUG
            cout<<"left:"<<find_left<<" "<<take_left<<", right:"<<find_right<<" "<<take_right<<endl;
#endif

            double max_pvalue=0;
            // subsampling MAX_NUM_CONTEXT of all possible contexts
            int leftNum = index - find_left;
            int rightNum = find_right - index;
            if (leftNum * rightNum > MAX_NUM_CONTEXT) {
                int contextCount = 0;

                while (contextCount < MAX_NUM_CONTEXT) { 
                    // random number [0, index-find_left-1]
                    // i1: [find_left+1, index]
                    int i1 = find_left + 1 + (int)(rand() / (((double)RAND_MAX) / (leftNum-1)));
                    // random number [0, find_right-index-1]
                    // i2: [index, find_right-1]
                    int i2 = index + (int)(rand() / (((double)RAND_MAX) / (rightNum-1)));
	  
                    int i_real=0;
                    for (int i=0;i<i2-i1+1;i++)
                        if (i1+i!=index)
                            intermediate_scores[i_real++].second = (int)all_seg[i1+i].score;
	  
                    qsort(intermediate_scores,i2-i1,sizeof(pair<int,int>),compare);
	  
                    for (int i=i2-i1-1;i>=0;i--)
                        intermediate_scores[i+1].second = intermediate_scores[i].second;
                    intermediate_scores[0].second=(int)all_seg[index].score;
	  
                    double total_score=0;
                    double fact=1;
                    double best_pvalue1=1;
                    int total_length=blocks.total_length;
                    if (i2<cur_seg-1)
                        total_length = all_seg[i2+1].align_start;
                    if (i1>0)
                        total_length -= all_seg[i1-1].align_end;
                    double mean_length = log(K*total_length*total_length)/H;
	  
                    for (int i=0;i<i2-i1+1;i++) {
                        total_score += lambda*intermediate_scores[i].second-log(K*(total_length-mean_length)*(total_length-mean_length));
                        fact *= (i+1);
                        double total_normalized_score = total_score + log(fact);
                        if (total_normalized_score > MAX_TOTAL_SCORE) total_normalized_score=MAX_TOTAL_SCORE;
	    
                        // If less than min_score or large number of segments
                        if ((total_normalized_score<MIN_TOTAL_SCORE)||(multi_pv[i+1][int(total_normalized_score-MIN_TOTAL_SCORE)]<0))
                            break;
                        if (multi_pv[i+1][int(total_normalized_score-MIN_TOTAL_SCORE)]*pow(2,i+1) < best_pvalue1)
                            best_pvalue1 = multi_pv[i+1][int(total_normalized_score-MIN_TOTAL_SCORE)]*pow(2,i+1);
                    }

                    if (best_pvalue1>max_pvalue) {
                        max_pvalue=best_pvalue1;
                        all_seg[index].left_end = i1-1;
                        all_seg[index].right_end = i2+1;
                    }

                    contextCount++;
                } // end while loop for MAX_NUM_CONTEXT
            }

            // if true, enumerate all context
            else {
                // for possible context containing the current segment
                for (int i1=find_left+1;i1<=index;i1++)   // start from the left of i1
                    for (int i2=index; i2<find_right;i2++) { // to the right of i2
	  
                        int i_real=0;
                        for (int i=0;i<i2-i1+1;i++)
                            if (i1+i!=index)
                                intermediate_scores[i_real++].second = (int)all_seg[i1+i].score;
	  
                        qsort(intermediate_scores,i2-i1,sizeof(pair<int,int>),compare);
	  
                        for (int i=i2-i1-1;i>=0;i--)
                            intermediate_scores[i+1].second = intermediate_scores[i].second;
                        intermediate_scores[0].second=(int)all_seg[index].score;
	  
                        double total_score=0;
                        double fact=1;
                        double best_pvalue1=1;
                        int total_length=blocks.total_length;
                        if (i2<cur_seg-1)
                            total_length = all_seg[i2+1].align_start;
                        if (i1>0)
                            total_length -= all_seg[i1-1].align_end;
                        double mean_length = log(K*total_length*total_length)/H;
	  
                        for (int i=0;i<i2-i1+1;i++) {
                            total_score += lambda*intermediate_scores[i].second-log(K*(total_length-mean_length)*(total_length-mean_length));
                            fact *= (i+1);
                            double total_normalized_score = total_score + log(fact);
                            if (total_normalized_score > MAX_TOTAL_SCORE) total_normalized_score=MAX_TOTAL_SCORE;
	    
                            // If less than min_score or large number of segments
                            if ((total_normalized_score<MIN_TOTAL_SCORE)||(multi_pv[i+1][int(total_normalized_score-MIN_TOTAL_SCORE)]<0))
                                break;
                            if (multi_pv[i+1][int(total_normalized_score-MIN_TOTAL_SCORE)]*pow(2,i+1) < best_pvalue1)
                                best_pvalue1 = multi_pv[i+1][int(total_normalized_score-MIN_TOTAL_SCORE)]*pow(2,i+1);
                            /*
                              #ifdef DEBUG
                              cout << "scores " << i << "\t" << intermediate_scores[i].second << "\t" << total_normalized_score << "\t" << multi_pv[i+1][int(total_normalized_score-MIN_TOTAL_SCORE)] << endl;
                              #endif*/
                        }
                        /*
                          #ifdef DEBUG       
                          cout << index << "\t" << find_left << "\t" << find_right << "\t" << i1 << "\t" << i2 << "\t" << total_length << "\t" << mean_length << "\t" << best_pvalue1 << endl;
                          #endif*/

                        if (best_pvalue1>max_pvalue) {
                            max_pvalue=best_pvalue1;
                            all_seg[index].left_end = i1-1;
                            all_seg[index].right_end = i2+1;
                            /*
                              #ifdef DEBUG
                              cout << "Improving max_pvalue " << max_pvalue << "\t" << index << "\t" << i1 << "\t" << i2 << endl;
                              #endif
                            */
                        }
                    } // end of loop for possible context (i1 & i2)
            }

            all_seg[index].nearest_left = index;
            all_seg[index].nearest_right = -1;
            all_seg[index].pv2 = max_pvalue;
      
            if ((find_left>=0)&&(max_pvalue < all_seg[find_left].pv)&&(take_left)) {
                max_pvalue = all_seg[find_left].pv;
                all_seg[index].nearest_left = find_left;
                all_seg[index].nearest_right = 0;
                all_seg[index].left_end = all_seg[find_left].left_end;
                all_seg[index].right_end = all_seg[find_left].right_end;
            }
      
            if ((find_right<cur_seg)&&(max_pvalue < all_seg[find_right].pv)&&(take_right)) {
                max_pvalue = all_seg[find_right].pv;
                all_seg[index].nearest_left = find_right;
                all_seg[index].nearest_right = 1;
                all_seg[index].left_end = all_seg[find_right].left_end;
                all_seg[index].right_end = all_seg[find_right].right_end;
            }
      
            all_seg[index].pv = max_pvalue;

#ifdef DEBUG      
            cout << "pv " << index << "\t" << all_seg[index].score << "\t" << all_seg[index].start << "\t" << all_seg[index].end << "\t" << all_seg[index].pv << "\t" << find_left  << "\t" << find_right << endl;
#endif

            cur_best_index++;
        } while (cur_best_index<cur_seg);
    }

#ifdef DEBUG
    cout <<"# Me here5\n";
#endif

    /*
      for (int i=0;i<cur_seg;i++)
      if (all_seg[i].pv>0.01) {
      int j=i;
      while ((j>=0)&&(all_seg[j].pv>0.01))
      j--;

      if (j>=0)
      all_seg[i].nearest_left=all_seg[j].end;
      else
      all_seg[i].nearest_left=0;

      j=i;
      while ((j<cur_seg)&&(all_seg[j].pv>0.01))
      j++;
      if (j<cur_seg)
      all_seg[i].nearest_right=all_seg[j].start;
      else
      all_seg[i].nearest_right=cur_seg;
      }
      else
      all_seg[i].nearest_left=all_seg[i].nearest_right=-1;   // -1 means position is undefined
    */


    cout << "# XXXXX num block_num block_score score start end pvalue "
         << "nearest_left nearest_right left_end right_end pv2 align_start align_end\n";

    cout << "# num_segments=" << cur_seg << endl;

    int grandtotal=0;
    int prev_good=0;
    for (int i=0;i<cur_seg;i++) {

        // Before first
        if ((i==0)&&(all_seg[i].align_start!=0)) {
            cout << "# Seg_fis " << i << "\t" << all_seg[i].block_num << "\t" 
                 << all_seg[i].block_score << "\t" << all_seg[i].score << "\t" << 0 << "\t" 
                 << all_seg[i].start-1 << "\t" << "1" << "\t" << "-1" << "\t";
    
            if (all_seg[i].nearest_right ==-1)
                cout << all_seg[i].start;
            else
                cout << all_seg[i].nearest_right;
            cout << "\t" << "-\t" << "-\t" << "-\t" << "0\t" << all_seg[i].align_start-1 << endl;
        }
    
        // Segment
        cout << "# Segment " << i << "\t" << all_seg[i].block_num << "\t" << all_seg[i].block_score 
             << "\t" << all_seg[i].score << "\t" << all_seg[i].start << "\t" << all_seg[i].end 
             << "\t" << all_seg[i].pv << "\t" << all_seg[i].nearest_left << "\t" 
             << all_seg[i].nearest_right << "\t" << all_seg[i].left_end << "\t" 
             << all_seg[i].right_end << "\t" << all_seg[i].pv2 << "\t" << all_seg[i].align_start 
             << "\t" << all_seg[i].align_end << endl;

        if (all_seg[i].pv<0.1) {
            if (all_seg[i].start-prev_good>=50)
                grandtotal+=all_seg[i].start-prev_good;
            prev_good=all_seg[i].end;
        }


        // Check next
        if ((i<cur_seg-1)&&(all_seg[i].align_end+1<all_seg[i+1].align_start)) {
            cout << "# Seg_nex " << i << "\t" << all_seg[i].block_num << "\t" << all_seg[i].block_score 
                 << "\t" << all_seg[i].score << "\t";

            if (all_seg[i].end+1<all_seg[i+1].start)
                cout<< all_seg[i].end+1 << "\t" << all_seg[i+1].start-1 << "\t" << "1" << "\t";
            // When it's gaps in reference 
            else 
                cout<< all_seg[i+1].start << "\t" << all_seg[i+1].start << "\t" << "1" << "\t";

            if (all_seg[i].nearest_left ==-1)
                cout << all_seg[i].end << "\t";
            else
                cout << all_seg[i].nearest_left << "\t";
      
            if (all_seg[i+1].nearest_right ==-1)
                cout << all_seg[i+1].start;
            else
                cout << all_seg[i+1].nearest_right;
            cout << "\t" << "-\t" << "-\t" << "-\t";
            cout << all_seg[i].align_end+1 << "\t" << all_seg[i+1].align_start-1 << endl;
        }


        // Check end
        if ((i==cur_seg-1)&&((all_seg[i].align_end) < (ALIGN_LEN-1))) {
            cout << "# Seg_end " << i << "\t" << all_seg[i].block_num << "\t" 
                 << all_seg[i].block_score << "\t" << all_seg[i].score << "\t" 
                 << all_seg[i].end+1 << "\t" << CHR_START+CHR_LEN-1 << "\t" << "1" << "\t";
            if (all_seg[i].nearest_left ==-1)
                cout << all_seg[i].end-CHR_LEN << "\t-1";
            else
                cout << all_seg[i].nearest_left << "\t-1";
            cout << "\t" << "-\t" << "-\t" << "-\t";
            cout << all_seg[i].align_end+1 << "\t" << ALIGN_LEN-1 <<endl;
        }
    }
    cout << "#grandtotal= " << grandtotal << endl;


    /*
      double sc_dist[1000][100];
      double len_dist[1000][100];
      double sc_dist_n[1000][100];
      double len_dist_n[1000][100];
      for (int i=0;i<1000;i++)
      for (int j=0;j<100;j++) {
      sc_dist[i][j]=0;
      len_dist[i][j]=0;
      sc_dist_n[i][j]=0;
      len_dist_n[i][j]=0;
      }
      for (int i=0;i<cur_seg;i++) {
      if (all_seg[i].pv>0.1) {
      int dist=all_seg[i].nearest;
      if (dist>999) dist=999;
      int score=int(all_seg[i].phast*100);
      if (score>99) score=99;
      for (int k1=0;k1<=dist;k1++)
      for (int k2=0;k2<=score;k2++) {
      sc_dist[k1][k2]++;
	  sc_dist_n[k1][k2]+= all_seg[i].end - all_seg[i].start;
      }
      }
      int dist=all_seg[i].end - all_seg[i].start;
      if (dist>999) dist=999;
      int score;
      if (all_seg[i].pv==0) score=99;
      else score = int(-10*log10(all_seg[i].pv));
      if (score>99) score=99;
      for (int k1=0;k1<=dist;k1++)
      for (int k2=0;k2<=score;k2++) {
      len_dist[k1][k2]++;
      len_dist_n[k1][k2]+= all_seg[i].end - all_seg[i].start;
      }
      }

      for (int i=999;i>=0;i--)
      for (int j=99;j>=0;j--) {
      sc_dist[i][j] /= sc_dist[0][0];
      sc_dist_n[i][j] /= sc_dist_n[0][0];
      len_dist[i][j] /= len_dist[0][0];
      len_dist_n[i][j] /= len_dist_n[0][0];
      }

      for (int i=0;i<1000;i++) {
      for (int j=0;j<100;j++)
      cout << i << "\t" << j/10.0 << "\t" << len_dist[i][j] << "\t" << sc_dist[i][j] << "\t" << len_dist_n[i][j] << "\t" << sc_dist_n[i][j] << endl;
      cout << endl;
      }

    */  
   
    delete[] all_seg;
  
    return 0;
} 



void compute_pvalue(Tree* node, blocks_type& blocks, double* branch_pvalue, 
                    int& branch_num, Tree* branch_ptr[], int doit, char* filename) {
    if (node->is_leaf_num < 0) {
        Tree* ll = node->left_subtree;
        Tree* rr = node->right_subtree;
        if (node == root) {
            branch_ptr[branch_num] = ll;
            if (branch_num == doit)
                branch_pvalue[branch_num++]  = compute_pvalue_branch(ll, node, blocks, filename);
            else
                branch_pvalue[branch_num++]  = 0;
        }else{
            branch_ptr[branch_num] = ll;
            if (branch_num == doit)
                branch_pvalue[branch_num++]  = compute_pvalue_branch(ll, node, blocks, filename);
            else
                branch_pvalue[branch_num++]  = 0;

            branch_ptr[branch_num] = rr;
            if (branch_num == doit)
                branch_pvalue[branch_num++]  = compute_pvalue_branch(rr, node, blocks, filename);
            else
                branch_pvalue[branch_num++]  = 0;
        }
        if (ll->is_leaf_num < 0) compute_pvalue(ll, blocks, branch_pvalue, branch_num, 
                                                branch_ptr, doit, filename);
        if (rr->is_leaf_num < 0) compute_pvalue(rr, blocks, branch_pvalue, branch_num, 
                                                branch_ptr, doit, filename);
    }
}

#endif
