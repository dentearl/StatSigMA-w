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
#ifndef ALPHA_H_
#define ALPHA_H_

#include "global.h"
#include <math.h>

#ifdef PROTEIN
  const int kAlphabetSize = 21;     // A=1  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V=20
  const int kGap = 21;
#include "protmat.h"
#else
  const int kAlphabetSize = 5;     // A=1, C=2, G=3, T=4, gap=5
  const int kGap = 5;
#endif // PROTEIN

const int kRandom = 0;
const int kNonOrtho = 100;
float g_back[kAlphabetSize + 1];
float evolution_prob(int node1, int node2, float branchlength) {
  /*
  if (node1 == GAP)
    return (node1 == node2) ? 1 : 0;
  else if (node1 == node2)
    return (exp(-SUBSTITUTION_PER_SITE*branchlength)+(1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[node2])+
 (exp(-SUBSTITUTION_PER_SITE*branchlength)-exp(-0.01*SUBSTITUTION_PER_SITE*branchlength))*GAP_BACKGRD;
  else if (node2 == GAP)
    return (1-exp(-0.01*SUBSTITUTION_PER_SITE*branchlength))*back[node2];
  else
    return (1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[node2];
  */
#ifdef PROTEIN
  if (node1 == kGap)
    return (node1 == node2) ? 1 : 0;
  if (node2 == kGap)
    return (1 - exp(-0.01 * branchlength)) * g_back[kGap];
  float m = 0.0;
  float t = 1.0;
  for (int i = 0; i < MAX_POWERS; i++) {
    m += aa_rate_matrix_pow[i][node1 - 1][node2 - 1] * t;
    t *= branchlength;
  }
  return (1 - ((1 - exp(-0.01 * branchlength)) * g_back[kGap])) * m;
#else
  // Amol's old code
  if (node1 == kGap)
    return (node1 == node2) ? 1 : 0;
  else if (node1 == node2)
    return (exp( -kSubstitutionPerSite * branchlength) + 
            (1 - exp(-kSubstitutionPerSite * branchlength)) * 
            g_back[node2]);
  else
    return (1 - exp(-kSubstitutionPerSite * branchlength)) * g_back[node2];
  /*
  // try to allow insertion, not successful!
  // Lead to problems in KA parameter estimation
  if (node1 == GAP)
    // gap -> gap
    if (node1 == node2) 
      return 1 - (4*(1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[GAP]);
    // gap -> A/C/G/T
    else
      return (1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[GAP];
  // A/C/G/T -> itself
  else if (node1 == node2)
    return (exp(-SUBSTITUTION_PER_SITE*branchlength)+(1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[node2]);
  // A/C/G/T -> other nucleotide or gap
  else
    return (1-exp(-SUBSTITUTION_PER_SITE*branchlength))*back[node2];
  */
#endif // PROTEIN
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
  case '-' : return kGap; break;
  case '.' : return kGap; break;
  default: return kRandom; break;
  }
}
#else
int alpha2int(char a) {
  switch (a) {
  case 'A' : return 1; break;
  case 'C' : return 2; break;
  case 'G' : return 3; break;
  case 'T' : return 4; break;
  case '-' : return kGap; break;  
  case '.' : return kGap; break;  
  default: return kRandom; break;
  }
}    
#endif // PROTEIN
#endif // ALPHA_H_
