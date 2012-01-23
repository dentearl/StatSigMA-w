#ifndef ALPHA_H
#define ALPHA_H

#include "global.h"
#include <math.h>


#ifdef PROTEIN
  const int ALPHABET_SIZE=21;     // A=1  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V=20
  const int GAP=21;
#include "protmat.h"

#else
  const int ALPHABET_SIZE=5;     // A=1, C=2, G=3, T=4, gap=5
  const int GAP=5;
#endif


const int RANDOM=0;
const int NON_ORTHO=100;
float back[ALPHABET_SIZE+1];



float evolution_prob(int node1, int node2, float branchlength) {
  /*
  if (node1==GAP)
    return (node1==node2) ? 1 : 0;
  else if (node1==node2)
    return (exp(-SUBSTITUTION_PER_SITE*branchlength)+(1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[node2])+
 (exp(-SUBSTITUTION_PER_SITE*branchlength)-exp(-0.01*SUBSTITUTION_PER_SITE*branchlength))*GAP_BACKGRD;
  else if (node2==GAP)
    return (1-exp(-0.01*SUBSTITUTION_PER_SITE*branchlength))*back[node2];
  else
    return (1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[node2];
  */

#ifdef PROTEIN

  if (node1==GAP)
    return (node1==node2) ? 1 : 0;
  if (node2==GAP)
    return (1-exp(-0.01*branchlength))*back[GAP];
  float m=0;
  float t=1;
  for (int i=0;i<MAX_POWERS;i++) {
    m+=aa_rate_matrix_pow[i][node1-1][node2-1]*t;
    t*=branchlength;
  }
  return (1-((1-exp(-0.01*branchlength))*back[GAP]))*m;


#else

  // Amol's old code
  if (node1==GAP)
    return (node1==node2) ? 1 : 0;
  else if (node1==node2)
    return (exp(-SUBSTITUTION_PER_SITE*branchlength)+(1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[node2]);
  else
    return (1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[node2];
  
  /*
  // try to allow insertion, not successful!
  // Lead to problems in KA parameter estimation
  if (node1==GAP)
    // gap -> gap
    if (node1==node2) 
      return 1 - (4*(1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[GAP]);
    // gap -> A/C/G/T
    else
      return (1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[GAP];
  // A/C/G/T -> itself
  else if (node1==node2)
    return (exp(-SUBSTITUTION_PER_SITE*branchlength)+(1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[node2]);
  // A/C/G/T -> other nucleotide or gap
  else
    return (1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[node2];
  */

#endif
}


#ifdef PROTEIN
int alpha2int(char a) {
  switch (a) {
  case 'A' : return 1; break;
  case 'R' : return 2; break;
  case 'N' : return 3; break;
  case 'D' : return 4; break;
  case 'C' : return 5; break;
  case 'Q' : return 6; break;
  case 'E' : return 7; break;
  case 'G' : return 8; break;
  case 'H' : return 9; break;
  case 'I' : return 10; break;
  case 'L' : return 11; break;
  case 'K' : return 12; break;
  case 'M' : return 13; break;
  case 'F' : return 14; break;
  case 'P' : return 15; break;
  case 'S' : return 16; break;
  case 'T' : return 17; break;
  case 'W' : return 18; break;
  case 'Y' : return 19; break;
  case 'V' : return 20; break;
  case '-' : return GAP; break;
  case '.' : return GAP; break;
  default: return RANDOM; break;
  }
}
#else
int alpha2int(char a) {
  switch (a) {
  case 'A' : return 1; break;
  case 'C' : return 2; break;
  case 'G' : return 3; break;
  case 'T' : return 4; break;
  case '-' : return GAP; break;  
  case '.' : return GAP; break;  
  default: return RANDOM; break;
  }
}    
#endif







#endif
