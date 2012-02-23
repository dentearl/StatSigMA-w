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
#include <stdlib.h>
#include <math.h>
#include <time.h>

const int kMaxIterations = 500;       /* Maximum number of iterations used in calculating K */

char* timeStamp() {
    time_t rawTime;
    struct tm *timeStr;
    time (&rawTime);
    timeStr = localtime(&rawTime);
    return asctime(timeStr);
}
int gcd(int a, int b) {
    int c;
    if (b < 0) 
        b = -b;
    if (a < b) { 
        c = a;
        a = b;
        b = c;
    }
    while(1) {
        c = a % b; 
        if (c == 0)
            return b;
        a = b;
        b = c;
    }
}

int karlin(int low, int high, double* pr, double* lambda, double* K, double* H) {
    int i, j, range, lo, hi;
    double up, neww, sum, Sum, av;
    double *p, *P, *ptrP, *ptr1;
    /* Check that scores and their associated probabilities are valid     */
    debug("Start Karlin: low=%d high=%d\n", low, high);
    debug("%s", timeStamp());
    if (0 <= low) {
        debug("Lowest score must be negative.\n");
        return 0;
    }
    for (i = range = high - low; i > -low && !pr[i]; --i)
        continue;
    if (i <= -low) {
        debug("A positive score must be possible.\n");
        return 0;
    }
    for (sum = i = 0; i <= range ; sum += pr[i++]) 
        if (pr[i] < 0) {
            debug("Negative probabilities not allowed.\n");
            return 0;
        }
    p = new double[range + 1];
    for (Sum = low , i = 0; i <= range; ++i) 
        Sum += i * (p[i] = pr[i] / sum);
	//////////////////////////////////////////////////////////////
    double *p_sparse = new double[range + 1];
    int *sc_sparse = new int[range + 1];
    for (i = 0; i < range + 1; ++i) {
        p_sparse[i] = 0.0;
        sc_sparse[i] = 0;
    }
	int num_sparse = 0;
    for (i = 0; i <= range; i++)
        if (p[i] > 0) {
            sc_sparse[num_sparse] = i;
            p_sparse[num_sparse++] = p[i];
        }
    if (Sum >= 0.1) {
        debug("Invalid (non-negative) expected score: %f\n", Sum);
        return 0;
    }
    /* Calculate the parameter lambda */
    up = 0.125 * (1.0 / kScoreBins);
    do {
        up *= 2;
        ptr1 = p;
        for (sum = i = 0; i <= range; ++i) 
            sum += *ptr1++ * exp(up * (low + i));
    } while (sum < 1.0);

    debug("Finished loop1\n");
    debug("%s", timeStamp());
    for (*lambda = j = 0; j < 100; ++j) {
        neww = (*lambda + up) / 2.0;
        ptr1 = p;
        for (sum = i = 0; i <= range; ++i) 
            sum += *ptr1++ * exp(neww * (low + i));
        if (sum < 1.0) 
            *lambda = neww;
        else 
            up = neww;
    }
    debug("Finished loop2\n");
    debug("%s", timeStamp());
    /* Calculate the pamameter K */
    ptr1 = p;
    for (av = 0, i = low; i <= high; ++i) 
        av += *ptr1++ *i * exp((*lambda) * i);
    *H = *lambda * av / log(2.0);
    Sum = lo = hi = 0;
    P = new double[kMaxIterations * range + 1];
	//////////////////////////////////////////////////////////////////
	for (i = 0; i < kMaxIterations * range + 1; i++)
        P[i] = 0;
    for (*P = sum = j = 1; j <= kMaxIterations && sum > 0.0001; Sum += sum /= j++) {
		ptrP = P;
		for (int jj = hi - low; jj >= 0;jj--) {
            double mm = ptrP[jj];
            if (mm > 0) {
                ptrP[jj] = 0;
                for (int kk = 0; kk < num_sparse; kk++)
                    ptrP[jj + sc_sparse[kk]] += mm * p_sparse[kk];
            }
        }
		hi += high;
		lo += low;
        for (sum = 0, i = lo; i; ++i) 
            sum += *++ptrP * exp((*lambda) * i);
        for ( ; i <= hi; ++i) sum += *++ptrP;
    }
    debug("Finished loop3, number of iterations=%d\n", j);
    debug("%s", timeStamp());
    for (i = low; !p[i - low]; ++i)
        continue;
    for (j = -i; i < high && j > 1; ) 
        if (p[++i - low]) 
            j = gcd(j, i);
    *K = (j * exp(-2 * Sum)) / (av * (1.0 - exp(- *lambda * j)));
    delete [] p;
    delete [] P;
    delete [] p_sparse;
    delete [] sc_sparse;
    return 1;               /* Parameters calculated successfully */
}

