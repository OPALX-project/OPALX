/* Sort.C
   Sorting routines

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    ------------
   10/03/2006    Copied from "numerical recipies"                Rene Bakker

   Last Revision:
   $Id: sort.C 85 2007-05-03 15:55:00Z bakker $
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Algorithms/bet/BetError.h"
#include "Algorithms/bet/math/sort.h"


/* ivector
   allocate an array of integers [0..n-1] and check memory,
   reallocate old vector if applicable */
static int *ivector(int n, int *old = NULL) {
    int *b;

    b = (int *) malloc(sizeof(int) * n);
    if(!b) {
        //writeBetError(errModeAll,errFatal,"Insufficient memory malloc %d bytes (sort.C)",sizeof(int)*n);
        writeBetError("Insufficient memory malloc in sort.C");
    }

    if(old) {
        int
        i,
        iMax = sizeof(old) / sizeof(int);

        for(i = 0; (i < iMax) && (i < n); i++) b[i] = old[i];
        free(old);
    }
    return b;
} /* ivector */


#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp
#define ISWAP(ia,ib) itemp=(ia);(ia)=(ib);(ib)=itemp
#define PSWAP(ap,bp) tempp=(ap);(ap)=(bp);(bp)=tempp
#define M 7

/* sort()
   Sorts an array arr[0..n-1] into ascending numerical order using the Quicksort
   algorithm. n is input; arr is replaced on output by its sorted rearrangement.
*/
void sortN(double *x, double **y, int n) {
    int
    nstack = 50;

    int
    i, ir = n - 1,
       j, k, l = 0;
    int
    *istack, jstack = 0;
    double
    a, temp,
     *bp = NULL, *tempp;

    istack = ivector(nstack);

    for(;;) {
        if(ir - l < M) {
            for(j = l + 1; j <= ir; j++) {
                a = x[j];
                if(y) bp = y[j];
                for(i = j - 1; i >= 0; i--) {
                    if(x[i] <= a) break;
                    x[i+1] = x[i];
                    if(y) y[i+1] = y[i];
                }
                x[i+1] = a;
                if(y) y[i+1] = bp;
            }
            if(!jstack) {
                free(istack);
                return;
            }
            ir = istack[jstack];
            l = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l + ir) >> 1;
            SWAP(x[k], x[l+1]);
            if(y) {PSWAP(y[k], y[l+1]);}
            if(x[l+1] > x[ir]) {
                SWAP(x[l+1], x[ir]);
                if(y) {PSWAP(y[l+1], y[ir]);}
            }
            if(x[l] > x[ir]) {
                SWAP(x[l], x[ir]);
                if(y) {PSWAP(y[l], y[ir]);}
            }
            if(x[l+1] > x[l]) {
                SWAP(x[l+1], x[l]);
                if(y) {PSWAP(y[l+1], y[l]);}
            }
            i = l + 1;
            j = ir;
            a = x[l];
            if(y) bp = y[l];
            for(;;) {
                do i++;
                while(x[i] < a);
                do j--;
                while(x[j] > a);
                if(j < i) break;
                SWAP(x[i], x[j]);
                if(y) {PSWAP(y[i], y[j]);}
            }
            x[l] = x[j];
            x[j] = a;
            if(y) {
                y[l] = y[j];
                y[j] = bp;
            }
            jstack += 2;
            if(jstack > nstack - 1) { // make stack bigger
                nstack += 50;
                istack = ivector(nstack, istack);
            }
            if(ir - i + 1 >= j - l) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j - 1;
            } else {
                istack[jstack] = j - 1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
} /* sortN() */

void sort2(double *x, double *y, int n) {
    int
    nstack = 50;
    int
    i, ir = n - 1,
       j, k, l = 0;
    int
    *istack, jstack = 0;
    double
    a, b =  0.0, temp;

    istack = ivector(nstack);

    for(;;) {
        if(ir - l < M) {
            for(j = l + 1; j <= ir; j++) {
                a = x[j];
                b = y[j];
                for(i = j - 1; i >= 0; i--) {
                    if(x[i] <= a) break;
                    x[i+1] = x[i];
                    y[i+1] = y[i];
                }
                x[i+1] = a;
                y[i+1] = b;
            }
            if(!jstack) {
                free(istack);
                return;
            }
            ir = istack[jstack];
            l = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l + ir) >> 1;
            SWAP(x[k], x[l+1]);
            SWAP(y[k], y[l+1]);
            if(x[l+1] > x[ir]) {
                SWAP(x[l+1], x[ir]);
                SWAP(y[l+1], y[ir]);
            }
            if(x[l] > x[ir]) {
                SWAP(x[l], x[ir]);
                SWAP(y[l], y[ir]);
            }
            if(x[l+1] > x[l]) {
                SWAP(x[l+1], x[l]);
                SWAP(y[l+1], y[l]);
            }
            i = l + 1;
            j = ir;
            a = x[l];
            if(y) b = y[l];
            for(;;) {
                do i++;
                while(x[i] < a);
                do j--;
                while(x[j] > a);
                if(j < i) break;
                SWAP(x[i], x[j]);
                SWAP(y[i], y[j]);
            }
            x[l] = x[j];
            x[j] = a;
            y[l] = y[j];
            y[j] = b;
            jstack += 2;
            if(jstack > nstack - 1) { // make stack bigger
                nstack += 50;
                istack = ivector(nstack, istack);
            }
            if(ir - i + 1 >= j - l) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j - 1;
            } else {
                istack[jstack] = j - 1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
} /* sort2() */

void isort2(double *x, int *y, int n) {
    int
    nstack = 50;

    int
    i, ir = n - 1,
       j, k, l = 0;
    int
    ib = -1, itemp,
      *istack, jstack = 0;
    double
    a, temp;

    istack = ivector(nstack);

    for(;;) {
        if(ir - l < M) {
            for(j = l + 1; j <= ir; j++) {
                a = x[j];
                if(y) ib = y[j];
                for(i = j - 1; i >= 0; i--) {
                    if(x[i] <= a) break;
                    x[i+1] = x[i];
                    if(y) y[i+1] = y[i];
                }
                x[i+1] = a;
                if(y) y[i+1] = ib;
            }
            if(!jstack) {
                free(istack);
                return;
            }
            ir = istack[jstack];
            l = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l + ir) >> 1;
            SWAP(x[k], x[l+1]);
            if(y) {ISWAP(y[k], y[l+1]);}
            if(x[l+1] > x[ir]) {
                SWAP(x[l+1], x[ir]);
                if(y) {ISWAP(y[l+1], y[ir]);}
            }
            if(x[l] > x[ir]) {
                SWAP(x[l], x[ir]);
                if(y) {ISWAP(y[l], y[ir]);}
            }
            if(x[l+1] > x[l]) {
                SWAP(x[l+1], x[l]);
                if(y) {ISWAP(y[l+1], y[l]);}
            }
            i = l + 1;
            j = ir;
            a = x[l];
            if(y) ib = y[l];
            for(;;) {
                do i++;
                while(x[i] < a);
                do j--;
                while(x[j] > a);
                if(j < i) break;
                SWAP(x[i], x[j]);
                if(y) {ISWAP(y[i], y[j]);}
            }
            x[l] = x[j];
            x[j] = a;
            if(y) {
                y[l] = y[j];
                y[j] = ib;
            }
            jstack += 2;
            if(jstack > nstack - 1) { // make stack bigger
                nstack += 50;
                istack = ivector(nstack, istack);
            }
            if(ir - i + 1 >= j - l) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j - 1;
            } else {
                istack[jstack] = j - 1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
} /* isort2() */

void sort1(double *x, int n) {
    int nstack = 50;

    int
    i, ir = n - 1,
       j, k, l = 0;
    int jstack = 0, *istack;
    double a, temp;

    istack = ivector(nstack);
    for(;;) {
        if(ir - l < M) {
            for(j = l + 1; j <= ir; j++) {
                a = x[j];
                for(i = j - 1; i >= 0; i--) {
                    if(x[i] <= a) break;
                    x[i+1] = x[i];
                }
                x[i+1] = a;
            }
            if(jstack == 0) break;
            ir = istack[jstack--];
            l = istack[jstack--];
        } else {
            k = (l + ir) >> 1;
            SWAP(x[k], x[l+1]);
            if(x[l+1] > x[ir]) {
                SWAP(x[l+1], x[ir]);
            }
            if(x[l] > x[ir]) {
                SWAP(x[l], x[ir]);
            }
            if(x[l+1] > x[l]) {
                SWAP(x[l+1], x[l]);
            }
            i = l + 1;
            j = ir;
            a = x[l];
            for(;;) {
                do i++;
                while(x[i] < a);
                do j--;
                while(x[j] > a);
                if(j < i) break;
                SWAP(x[i], x[j]);
            }
            x[l] = x[j];
            x[j] = a;
            jstack += 2;
            if(jstack > nstack - 1) { // make stack bigger
                nstack += 50;
                istack = ivector(nstack, istack);
            }
            if(ir - i + 1 >= j - l) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j - 1;
            } else {
                istack[jstack] = j - 1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
    free(istack);
}

#undef M
#undef SWAP
#undef ISWAP
#undef PSWAP

