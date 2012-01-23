#include "global_combine.h"
#include "tree_combine.h"


float get_PV(int i) {
  return cur_pv[i];
}


void set_PV(int i, float pv) { 
  cur_pv[i]=pv;
}


bool get_Present(int i, int sp) { 
  return (present[sp][i]); 
}

 
bool is_NonePresent(int i) { 
  bool res = true;
  for (int s=0; s<branch_num; s++)
    if (present[s][i]) {
      res = false;
      break;
    }
  return res; 
}


void set_Present(int i, int sp, bool cur_present) { 
  present[sp][i] = cur_present;
}


bool get_Gap(int i, int sp) { 
  return (gaps[sp][i]); 
}


void set_Gap(int i, int sp, bool cur_gaps) {
  gaps[sp][i] = cur_gaps;
}


// Initialize worst_set for all species at position i.
void init_WorstSet(int i) {
  for (int s=0; s<(branch_num+1); s++)
    worst_set[s][i] = false;
}


void add_WorstBranch(int i,int sp) {
  worst_set[sp][i] = true;
  if (sp<branch_num)
    for (int s=0; s<branch_num; s++)
      // check if sp is on the path from s to human
      if ((get_Present(i,s)) && (is_OnPathToHuman(s,sp)))
	worst_set[s][i] = true;
}


bool isin_WorstSet(int i,int sp) {
  return (worst_set[sp][i]);
}


void read_present() {
  // Initialize
  for (int i=0; i<ALIGN_LEN; i++)
    for (int s=0; s<branch_num; s++)
      set_Present(i, s, false);
  int cur_pos = 0;
  int cur_len = 0;
  bool cur_present[branch_num];
  for (int s=0; s<branch_num; s++)
    cur_present[s] = false; 

  ifstream ifs(mafFile);

  char buffer[MAX_BLOCK_SIZE];
  bool firstBlock = true;
  int start_pos = 0;
  while (ifs) {
    buffer[0]='\0';
    ifs.getline(buffer,MAX_BLOCK_SIZE);
    if ((strncmp(buffer,"a score=",8)==0) || (!ifs)) {
      insert_inner_branches(cur_present);
      for (int i=0; i<cur_len; i++)
	for (int s=0; s<branch_num; s++)
	  set_Present(cur_pos+i, s, cur_present[s]);

      // Reset
      cur_pos += cur_len;
      cur_len=0;
      for (int s=0; s<branch_num; s++)
	cur_present[s] = false;
    }
    else {
      char* a;
      if (buffer[0]=='s') {
	// Checking for human
	bool is_human=false;
	if (strncmp(buffer+2,REF_SPECIES,strlen(REF_SPECIES))==0)
	  is_human=true;

	a = strtok(buffer," ");
	a = strtok(NULL," ");
	int species_line=name2branch(a);

	a = strtok(NULL," ");
	a = strtok(NULL," ");
	a = strtok(NULL," ");
	a = strtok(NULL," ");
	a = strtok(NULL," ");
	if (is_human)
	  cur_len = strlen(a);
	if (a[0]!='N')
	  cur_present[species_line] = true;
      }
    }
  }
  ifs.close();
}


void read_gaps() {
  for (int i=0; i<ALIGN_LEN; i++)
    for (int s=0; s<branch_num; s++)
      set_Gap(i, s, false);
  
  int cur_pos=0;
  int cur_len=0;
  char buffer[MAX_BLOCK_SIZE];
  bool **cur_gaps;
  cur_gaps = (bool **)malloc(sizeof(bool *) * branch_num);
  for (int s=0; s<branch_num; s++)
    cur_gaps[s] = (bool *)malloc(sizeof(bool) * MAX_BLOCK_SIZE);
  for (int s=0; s<branch_num; s++)
    for (int i=0; i<MAX_BLOCK_SIZE; i++)
      cur_gaps[s][i]=true;

  ifstream ifs(mafFile);
  while (ifs) {
    buffer[0]='\0';
    ifs.getline(buffer,MAX_BLOCK_SIZE);
    if ((strncmp(buffer,"a score=",8)==0) || (!ifs)) {      
      bool tmp_gaps[branch_num];
      for (int i=0; i<cur_len; i++) {
	for (int s=0; s<branch_num; s++)
	  tmp_gaps[s] = cur_gaps[s][i];

	check_gapped_branches(tmp_gaps);
	for (int s=0; s<branch_num; s++)
	  set_Gap(cur_pos+i, s, tmp_gaps[s]);
      }

      // Reset
      cur_pos += cur_len;
      for (int s=0; s<branch_num; s++)
	for (int i=0; i<MAX_BLOCK_SIZE; i++)
	  cur_gaps[s][i]=true;
      cur_len=0;
    }
    else {
      char* a;
      if (buffer[0]=='s') {
	// Checking for reference species
	bool is_human=false;
	if (strncmp(buffer+2,REF_SPECIES,strlen(REF_SPECIES))==0)
	  is_human=true;

	a = strtok(buffer," ");
	a = strtok(NULL," ");
	int species_line = name2branch(a);

	a = strtok(NULL," ");
	a = strtok(NULL," ");
	a = strtok(NULL," ");
	a = strtok(NULL," ");
	a = strtok(NULL," ");
	if (is_human)
	  cur_len = strlen(a);
	// a includes gaps in the reference species here.
	for (int i=0; i<strlen(a); i++)
	  if ((a[i]!='-') && (a[i]!='N'))
	    cur_gaps[species_line][i] = false;
      }
    }
  }

  ifs.close();
  for (int s=0; s<branch_num; s++)
    free(cur_gaps[s]);
  free(cur_gaps);
}


bool read_segment(ifstream &ifs,int& start,int& end,int& left_nearest,int& right_nearest,double& pv) {
  char text[10000];
  int segment,block,block_size;
  double score;
  if (!ifs) return false;

  char errorMsg[50] = "# number of sets of parameter computed= 0";
  text[0] = text[2] = '\0';
  ifs.getline(text,10000);
  while ((ifs) && (text[2]!='S')) {
    // if no parameter estimated, abort this file.
    if (strcmp(text, errorMsg) == 0) {
      return false;
    }
    text[0]=text[2]='\0';
    ifs.getline(text,10000);
  }

  if (text[2]=='S') {
    char* a = strtok(text+2," \t");
    a = strtok(NULL," \t");
    segment = atoi(a);
    a = strtok(NULL," \t");
    block = atoi(a);
    a = strtok(NULL," \t");
    block_size = atoi(a);
    a = strtok(NULL," \t");
    score = atof(a);
    a = strtok(NULL," \t");
    a = strtok(NULL," \t");
    a = strtok(NULL," \t");
    pv = atof(a);
    a = strtok(NULL," \t");
    left_nearest = atoi(a);
    a = strtok(NULL," \t");
    right_nearest = atoi(a);
    a = strtok(NULL," \t");
    a = strtok(NULL," \t");
    a = strtok(NULL," \t");
    a = strtok(NULL," \t");
    start = atoi(a);
    a = strtok(NULL," \t");
    end = atoi(a);
    if (end>ALIGN_LEN) end=ALIGN_LEN;

    return true;
  }
  return false;
}


// Identifying good regions  (removed gap_check, because now gaps are counted after outputting good regions)
void output_good_regions() {
  int start = -1;
  for (int i=0; i<ALIGN_LEN; i++) {
    if (get_PV(i) < PVTHRESH_GOOD) {
      if (start == -1)
	start = i;
    }
    else if (start != -1) {
      int currLen = i - start;

      if (currLen >= MIN_SEG_SIZE) {
	// the first and the last human coordinates in current region                
	int c = start;
	while ((hgCoor[c] == -1) && (c < ALIGN_LEN - 1))
	  c++;
	int firstHgCoor = hgCoor[c];
	
	c = start + currLen - 1;
	while ((hgCoor[c] == -1) && (c > 0))
	  c--;
	int lastHgCoor = hgCoor[c];
	
	cout << "my_good_region " << " " << start << " " << currLen <<  " " <<hgChr<< " " << firstHgCoor+1 << " " << lastHgCoor+1 << endl;
      }

      start = -1;
    }
  }
}


void read_output_pv(const char* branch_multiplier) {
  for (int sp=0;sp<branch_num;sp++) {
    int start,end,left_nearest,right_nearest;
    double pv;
    char filename[1000];
    sprintf(filename,"%s/%s.%d_%s.out",
	    SigMAwOutDir, SigMAwOutPrefix, sp, branch_multiplier);
    //cout<<filename<<endl;

    ifstream ifs(filename);
    while (read_segment(ifs,start,end,left_nearest,right_nearest,pv))
      for (int i=start; i<=end; i++)
	// get the max PValue for all branches
	if ((get_Present(i,sp)) && (get_PV(i) < pv))
	  set_PV(i,pv);
    ifs.close();
  }
}


void read_output_worstbranch(const char* branch_multiplier) {
  // get the tree ID for current tree
  char treeID;
  if (strcmp(branch_multiplier, "100") == 0) treeID = '2';
  if (strcmp(branch_multiplier, "1") == 0) treeID = '1';
  if (strcmp(branch_multiplier, "0.01") == 0) treeID = '0';

  for (int sp=0;sp<branch_num;sp++) {
    int start,end,left_nearest,right_nearest;
    double pv;
    char filename[1000];

    sprintf(filename,"%s/%s.%d_%s.out",
	    SigMAwOutDir, SigMAwOutPrefix, sp, branch_multiplier);
    ifstream ifs(filename);
    while (read_segment(ifs,start,end,left_nearest,right_nearest,pv))
      for (int i=start; i<=end; i++) {
	// when the branch is present at the position, 
	// the pvalue of the branch equals to min(max pvalue), 
	// and current tree is the reponsible tree for such pvalue.
	if (get_Present(i,sp) && (get_PV(i) == pv) && (tree_pv[i] == treeID))
	  add_WorstBranch(i,sp);
      }
    ifs.close();
  }
}


void read_pv_dist(char* branch_multiplier,int sp) {

  double pv_dist[10];
  for (int i=0;i<10;i++)
    pv_dist[i]=0;
  int total=0;

  int start,end,left_nearest,right_nearest;
  double pv;
  char filename[1000];
  sprintf(filename,"%s/%s.%d_%s.out",
	  SigMAwOutDir, SigMAwOutPrefix, sp, branch_multiplier);

  ifstream ifs(filename);
  while (read_segment(ifs,start,end,left_nearest,right_nearest,pv))
    if (end-start >= 0) {
      int total_start_end=0;
      for (int j=start;j<end;j++)
	if (get_Present(j,sp))
	  total_start_end++;
      if (pv<1e-9) pv=1e-9;
      pv = -1*log10(pv);
      pv_dist[int(pv)] += total_start_end;
      total+=total_start_end;
    }
  ifs.close();
  for (int i=0;i<10;i++)
    cout << "my_pv " << sp << " " << i << " " << pv_dist[i]/float(total) << " " << pv_dist[i] << endl; 
}


void scan_alignFile() {

  ALIGN_LEN = 0;
  int cur_len = 0;
  char buffer[MAX_BLOCK_SIZE];
  bool firstBlock = true;
  int currStart = 0;
  int currNtNum = 0;

  ifstream ifs(mafFile);
  while (ifs) {
    buffer[0]='\0';
    ifs.getline(buffer,MAX_BLOCK_SIZE);
    if ((strncmp(buffer,"a score=",8)==0) || (!ifs)) {
      ALIGN_LEN += cur_len;
      cur_len=0;
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
	  currStart = atoi(a);
	  if (firstBlock) {
	    CHR_START = atoi(a);
	    firstBlock = false;
	  }
	  a = strtok(NULL," ");
	  currNtNum = atoi(a);
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



void map_humanCoordinates() {
    int cumStart = 0;
    int cur_len = 0;
    char buffer[MAX_BLOCK_SIZE];
    bool firstBlock = true;
    char *orient;
    int currHgPos = 0;
    int chrLen = 0;

    ifstream ifs(mafFile);
    while (ifs) {
        buffer[0]='\0';
        ifs.getline(buffer,MAX_BLOCK_SIZE);
        if ((strncmp(buffer,"a score=",8)==0) || (!ifs)) {
            cumStart += cur_len;
            cur_len=0;
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
                    // get the chr info
                    if (firstBlock) {
                        char* temp = strstr(a, "chr");
                        strcpy(hgChr, temp);
                        firstBlock = false;
                    }

                    a = strtok(NULL," ");
                    currHgPos = atoi(a);

                    a = strtok(NULL," ");
                    orient = strtok(NULL," ");
	  
                    a = strtok(NULL," ");
                    chrLen = atoi(a);

                    a = strtok(NULL," ");
                    cur_len = strlen(a);

                    if (strcmp(orient, "-") == 0) {
                        cout<<"Warning: there is negative strand in human!"<<endl;
                        currHgPos = chrLen - currHgPos - 1;
                    }

                    for (int i=0; i<cur_len; i++) {
                        // 'N' is actually genomic positions
                        if (a[i]=='-')
                            hgCoor[cumStart+i] = -1;
                        else {
                            hgCoor[cumStart+i] = currHgPos;
                            currHgPos++;
                        }
                    }
                }
            }
        }
    }
    ifs.close();
}



int main(int argc, char* argv[]) {

  strcpy(mafFile, argv[1]);
  strcpy(SigMAwOutDir, argv[2]);
  strcpy(SigMAwOutPrefix, argv[3]);
  
  scan_alignFile();
  cout<<"alignmnt length: "<<ALIGN_LEN<<endl;

  cur_pv = (float *)malloc(sizeof(float) * ALIGN_LEN);
  temp_pv = (float *)malloc(sizeof(float) * ALIGN_LEN);
  // record the tree ID that has min(max PValue)
  tree_pv = (char *)malloc(sizeof(char) * ALIGN_LEN);

  read_newick_tree();
  assert(branch_num * ALIGN_LEN < 1e9);
  cout<<"read tree"<<endl;  

  // memory allocate
  present = (bool **)malloc(sizeof(bool *) * branch_num);
  for (int s=0; s<branch_num; s++) 
    present[s] = (bool *)malloc(sizeof(bool) * ALIGN_LEN);

  // read in present[] for each branch-position,
  // indicating if the branch is "present" at the  position
  read_present();
  cout<<"read present"<<endl;

  // read in curr_pv[], the min(max pvalue over branches) over trees.
  // Firstly, identify the max pvalue among different branches. 
  // Then, identify the min pvalue among the three trees.
  // Using temp_pv as the structure to store temporary pvalue (as the 3 trees are read)
  
  // initialize
  for (int i=0;i<ALIGN_LEN;i++) {
    set_PV(i,0);
    temp_pv[i] = 1;
    tree_pv[i] = 'n';
  }

  // for tree "100", ID is 2.
  read_output_pv("100");
  for (int i=0; i<ALIGN_LEN; i++) {
    if ((temp_pv[i] > cur_pv[i])&&(!is_NonePresent(i))) {
      temp_pv[i] = cur_pv[i];
      tree_pv[i] = '2';
    }
    set_PV(i,0);
  }

  // for tree "0.01", ID is 0
  read_output_pv("0.01");
  for (int i=0; i<ALIGN_LEN; i++) {
    if ((temp_pv[i] > cur_pv[i])&&(!is_NonePresent(i))) {
      temp_pv[i] = cur_pv[i];
      tree_pv[i] = '0';
    }
    set_PV(i,0);
  }

  // for tree "1", ID is 1.
  read_output_pv("1");
  for (int i=0; i<ALIGN_LEN; i++) {
    if (temp_pv[i] < cur_pv[i])
      cur_pv[i]=temp_pv[i];
    else
      tree_pv[i] = '1';
  }
  free(temp_pv);
  cout<<"identified pvalues"<<endl;

  // Identify the worst set that are branches having pvalue == curr_pv at each position.
  // memory allocate
  worst_set = (bool **)malloc(sizeof(bool *) * (branch_num+1));
  for (int s=0; s<(branch_num+1); s++)
    worst_set[s] = (bool *)malloc(sizeof(bool) * ALIGN_LEN);

  for (int i=0; i<ALIGN_LEN; i++)
    init_WorstSet(i);
  read_output_worstbranch("100");
  read_output_worstbranch("0.01");
  read_output_worstbranch("1");
  cout<<"identified worst_set"<<endl;

  // Bonferonni Correction for  the curr_pv (due to three trees different)
  for (int i=0; i<ALIGN_LEN; i++) {
    float f = get_PV(i);
    if (f>0.33)
      set_PV(i,1);
    else
      set_PV(i,f*3);
  }
  

  int count[branch_num+1][3];
  for (int sp=0;sp<branch_num+1;sp++) {
    count[sp][0]=count[sp][1]=count[sp][2]=0;
    for (int i=0; i<ALIGN_LEN; i++)
      if (sp<branch_num) {
	if (get_Present(i,sp))
	  count[sp][0]++;
      }
      else {
	if (!is_NonePresent(i)) 
	  count[sp][0]++;
      }
  }
    

  // memory allocate
  for (int s=0; s<branch_num; s++)
    free(present[s]);
  free(present);
  gaps = (bool **)malloc(sizeof(bool *) * branch_num);
  for (int s=0; s<branch_num; s++)
    gaps[s] = (bool *)malloc(sizeof(bool) * ALIGN_LEN);

  // Reading gaps
  read_gaps();
  cout<<"read gaps"<<endl;


  // map alignment coordinates to human coordinates
  hgCoor = (int *)malloc(sizeof(int) * ALIGN_LEN);
  map_humanCoordinates();
  cout<<"got human coordinates map"<<endl;

  // Identigying good regions
  output_good_regions();

  // Identifying bad regions
  int start;
  for (int sp=0; sp<branch_num; sp++) {
    start=-1;

    for (int i=0; i<ALIGN_LEN; i++) {
      // when current position is a bad one
      if ((get_PV(i) >= PVTHRESH_BAD) && (isin_WorstSet(i,sp))) {
	// starting position
	if (start == -1)
	  start = i;
      }
      // when current position is a good one, output the region before it if the region is long enough.
      else if (start != -1) {
	int currLen = i-start;

	if (currLen >= MIN_SEG_SIZE) {
	  // count the number of gaps in the region.
	  int gap_count=0;
	  for (int j=start; j<i; j++)
	    if (get_Gap(j, sp))
	      gap_count++;

	  // the first and the last human coordinates in current region
	  int c = start;
	  while ((hgCoor[c] == -1) && (c < ALIGN_LEN - 1))
	    c++;
	  int firstHgCoor = hgCoor[c];
	  
	  c = start + currLen - 1;
	  while ((hgCoor[c] == -1) && (c > 0))
	    c--;
	  int lastHgCoor = hgCoor[c];

	  // output
	  if (gap_count < 0.5*currLen) {
	    // ucsc browser coordinate start with 1, hgCoor+1
	    cout << "my_bad_region " << start << " " << sp << " " << currLen << " " <<hgChr<< " " << firstHgCoor+1 << " " << lastHgCoor+1 << endl; 
	    count[sp][1] += currLen;
	    count[sp][2]++;
	    for (int my_pos=start; my_pos<i; my_pos++) {
	      // using 'branch_num' in the 'worst_set' to indicate this position is a bad one.
	      add_WorstBranch(my_pos,branch_num);
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
  start=-1;
  for (int i=0; i<ALIGN_LEN; i++) {
    if (isin_WorstSet(i,branch_num)) {
      if (start==-1)
	start=i;
    }
    else if (start!=-1) {
      count[branch_num][1] += i-start;
      count[branch_num][2]++;
      start=-1;
    }
  }
  
  // Output summary
  for (int sp=0; sp<branch_num+1; sp++)
    printf("my_count %d %d %d(%4.2f%%) %d\n",sp,count[sp][0],count[sp][1],count[sp][1]*100.0/float(count[sp][0]),count[sp][2]);
  cout<<"printed count"<<endl;  
}
