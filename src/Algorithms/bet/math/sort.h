/* sort.h
   sorting routines

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   09-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: sort.h 85 2007-05-03 15:55:00Z bakker $
*/


#ifndef _SORT_DEF
#define _SORT_DEF

/* sort()
   Sorts an array arr[0..n-1] into ascending numerical order using the Quicksort 
   algorithm. n is input; arr is replaced on output by its sorted rearrangement.
*/
void sort1(
  double *,            // array to sort
  int);      // number of elements in array

void sort2(
  double *,            // array to sort
  double *,            // array to sort along
  int);      // number of elements in array

void isort2(
  double *,            // array to sort
  int *,               // index array to sort along
  int);      // number of elements in array

void sortN(
  double *,            // array to sort
  double **,           // matrix to sort with x
  int);      // number of elements in array

#endif
