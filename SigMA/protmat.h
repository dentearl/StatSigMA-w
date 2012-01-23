#ifndef PROTMAT_H
#define PROTMAT_H

#include <fstream.h>
#include <iostream.h>

const int MAX_POWERS=10;

double** aa_rate_matrix;
double* aa_back;
double** aa_rate_matrix_pow[MAX_POWERS];


void matrix_multiply(double** target, double** src1, double** src2) {
  for (int i=0;i<20;i++)
    for (int j=0;j<20;j++) {
      target[i][j]=0;
      for (int k=0;k<20;k++)
	target[i][j] += src1[i][k]*src2[k][j];
    }
}


void compute_powers() {
  for (int i=0;i<MAX_POWERS;i++) {
    if (i==0) {
      for (int k1=0;k1<20;k1++)
	for (int k2=0;k2<20;k2++)
	  if (k1==k2)
	    aa_rate_matrix_pow[i][k1][k2]=1;
	  else
	    aa_rate_matrix_pow[i][k1][k2]=0;
    }
    else {
      matrix_multiply(aa_rate_matrix_pow[i],aa_rate_matrix_pow[i-1],aa_rate_matrix);
      for (int k1=0;k1<20;k1++)
	for (int k2=0;k2<20;k2++)
	  aa_rate_matrix_pow[i][k1][k2]/=(i+1);      
    }
  }
}


void read_aa_matrix(char* filename) {

  aa_back = new double[20];
  aa_rate_matrix = new double*[20];
  for (int i=0;i<20;i++)
    aa_rate_matrix[i] = new double[20];
  for (int i=0;i<MAX_POWERS;i++) {
    aa_rate_matrix_pow[i] = new double*[20];
    for (int j=0;j<20;j++)
      aa_rate_matrix_pow[i][j] = new double[20];
  }

  ifstream ifs(filename);
  int i=1;
  int j=0;
  while ((ifs)&&(i<20)) {
    ifs >> aa_rate_matrix[i][j];
    aa_rate_matrix[j][i] = aa_rate_matrix[i][j];
    j++;
    if (i==j) {
      i++;
      j=0;
    }
  }

  i=0;
  while (ifs) {
    ifs >> aa_back[i];
    i++;
  }

  for (i=0;i<20;i++) {
    double total=0;
    for (j=0;j<20;j++)
      if (i!=j) total+=aa_rate_matrix[i][j];
    aa_rate_matrix[i][i]=-1*total;
  }

  for (i=0;i<20;i++)
    for (j=0;j<20;j++)
      aa_rate_matrix[i][j]*=aa_back[j];

  double total=0;
  for (i=0;i<20;i++)
    total += aa_rate_matrix[i][i]*aa_back[i];

  for (i=0;i<20;i++)
    for (j=0;j<20;j++)
      aa_rate_matrix[i][j] /= (-1*total);

  compute_powers();
}




void transition_probability(double dist) {
  for (int k1=0;k1<20;k1++) {
    for (int k2=0;k2<20;k2++) {
      float m=0;
      float t=1;
      for (int i=0;i<MAX_POWERS;i++) {
	m+=aa_rate_matrix_pow[i][k1][k2]*t;
	t*=dist;
      }
      cout << m << " ";
    }
    cout << endl;
  }

  cout << "_________________" << endl;
  for (int k1=0;k1<20;k1++) {
    double tot=0;
    for (int k2=0;k2<20;k2++) {
      float m=0;
      float t=1;
      for (int i=0;i<MAX_POWERS;i++) {
	m+=aa_rate_matrix_pow[i][k1][k2]*t;
	t*=dist;
      }
      tot += m;
    }
    cout << tot << endl;
  }

  double tot=0;
  for (int k1=0;k1<20;k1++)
    tot += aa_back[k1];
  cout << "FReq " << tot << endl;
}

#endif
