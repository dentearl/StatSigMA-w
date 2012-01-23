#include "global_combine.h"
#include "tree_combine.h"


int main(int argc, char* argv[]) {

  read_newick_tree();
  cerr<<"read tree"<<endl;  

  for (int i=0; i<53; i++) {
    
    int br = name2branch(species_name[i]);

    if (br > -1)
      cout<<br<<" "<<species_name[i]<<endl;
  }
}
