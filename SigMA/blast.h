#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAXIT 500       /* Maximum number of iterations used in calculating K */


char* timeStamp() {
  time_t rawTime;
  struct tm *timeStr;

  time ( &rawTime );
  timeStr = localtime ( &rawTime );

  return asctime(timeStr);
}



int gcd(int a,int b)
{
        int c;

        if (b<0) b= -b;
        if (b>a) { c=a; a=b; b=c; }
        for (;b;b=c) { c=a%b; a=b; }
        return a;
}


int karlin(int low,int high,double* pr,double* lambda,double* K,double* H)
{
        int     i,j,range,lo,hi,first,last;
        double  up,neww,sum,Sum,av,beta,ftemp;
        double  *p,*P,*ptrP,*ptr1,*ptr2; //,exp();
	//        char    *calloc();

        /* Check that scores and their associated probabilities are valid     */

#ifdef DEBUG 
        cout <<"Start karlin: low=" <<low<<" high="<<high<<endl; 
        cout<<timeStamp()<<endl;
#endif

        if (low>=0) {
#ifdef DEBUG 
	  cout << "# Lowest score must be negative.\n";
#endif
                return 0;
        }

        for (i=range=high-low;i> -low && !pr[i];--i);

        if (i<= -low) {
#ifdef DEBUG 
	  cout << "# A positive score must be possible.\n";
#endif
                return 0;
        }

        for (sum=i=0;i<=range;sum+=pr[i++]) if (pr[i]<0) {
#ifdef DEBUG 
	  cout << "# Negative probabilities not allowed.\n";
#endif
                return 0;
        }

	/*
#ifdef DEBUG
        if (sum<0.99995 || sum>1.00005) 
	  cout << "Probabilities sum to " << sum << " (ideally should be 1) Normalizing.\n";
#endif
	*/
        p= (double *) calloc(range+1,sizeof(double));
        for (Sum=low,i=0;i<=range;++i) Sum+=i*(p[i]=pr[i]/sum);


	//////////////////////////////////////////////////////////////
        double* p_sparse= (double *) calloc(range+1,sizeof(double));
        int* sc_sparse= (int *) calloc(range+1,sizeof(int));
	int num_sparse=0;
        for (i=0;i<=range;i++)
	  if (p[i]>0) {
	    sc_sparse[num_sparse]=i;
	    p_sparse[num_sparse++]=p[i];
	  }

        if (Sum>=0.1) {
#ifdef DEBUG 
	  cout << "Invalid (non-negative) expected score: " << Sum << endl;
#endif
                return 0;
        }

        /* Calculate the parameter lambda */

        up=0.125*(1.0/SCORE_BINS);
        do {
                up*=2;
                ptr1=p;
		//		beta=exp(up);
		//		ftemp=exp(up*(low-1));
                for (sum=i=0;i<=range;++i) sum+= *ptr1++ * exp(up*(low+i));
        }
        while (sum<1.0);

#ifdef DEBUG
        cout << "Finish loop1" << endl;
        cout << timeStamp() << endl;
#endif

        for (*lambda=j=0;j<100;++j) {
                neww=(*lambda+up)/2.0;
		//		beta=exp(neww);
		//		ftemp=exp(neww*(low-1));

                ptr1=p;
                for (sum=i=0;i<=range;++i) 
		  sum+= *ptr1++ * exp(neww*(low+i));
                if (sum<1.0) *lambda=neww;
                else up=neww;
        }

#ifdef DEBUG
        cout << "Finish loop2" << endl;
        cout << timeStamp() << endl;
#endif
        /* Calculate the pamameter K */

	//        ftemp=exp(*lambda*(low-1));
        ptr1=p;
        for (av=0,i=low;i<=high;++i) av+= *ptr1++ *i*exp((*lambda)*i);
         *H= *lambda*av/log(2.0);


        Sum=lo=hi=0;
        P= (double *) calloc(MAXIT*range+1,sizeof(double));
	//////////////////////////////////////////////////////////////////
	for (i=0;i<MAXIT*range+1;i++)
	  P[i]=0;

        for (*P=sum=j=1;j<=MAXIT && sum>0.0001;Sum+=sum/=j++) {
	  /*
                first=last=range;
                for (ptrP=P+(hi+=high)-(lo+=low);ptrP>=P;*ptrP-- =sum) {
                        ptr1=ptrP-first;
			ptr2=p+first;
			for (sum=0,i=first;i<=last;++i) sum+= *ptr1-- * *ptr2++;
                        if (first) --first;
                        if (ptrP-P<=range) --last;
                }
	  */
	  //////////////////////////////////////////////////////////////////

		ptrP=P;
		for (int jj=hi-low;jj>=0;jj--) {
		  double mm=ptrP[jj];
		  if (mm>0) {
		    ptrP[jj]=0;
		    for (int kk=0;kk<num_sparse;kk++)
		      ptrP[jj+sc_sparse[kk]] += mm*p_sparse[kk];
		  }
                }
		hi+=high;
		lo+=low;
		/*
		for (int jj=0;jj<=hi-lo;jj++)
		  if (P[jj]>0)
		    cout << "# lambda_ptrP " << jj << "\t" << P[jj] << endl;
		*/

		//                ftemp=exp(*lambda*(lo-1));
                for (sum=0,i=lo;i;++i) sum+= *++ptrP * exp((*lambda)*i);
                for (;i<=hi;++i) sum+= *++ptrP;
        }

#ifdef DEBUG
        cout << "Finish loop3, # of iterations = " <<j << endl;
        cout << timeStamp() << endl;
#endif

        for (i=low;!p[i-low];++i);
        for (j= -i;i<high && j>1;) if (p[++i-low]) j=gcd(j,i);
        *K = (j*exp(-2*Sum))/(av*(1.0-exp(- *lambda*j)));

        free(p);
        free(P);
        free(p_sparse);
        free(sc_sparse);
        return 1;               /* Parameters calculated successfully */
}

