/* Driver for routine sort */

#include <stdio.h>
#include <stdlib.h>

#include "error.h"
#include "sort.h"

#define MAXSTR 80
#define NP 100
#define NQ   4

int main(void)
{
  char txt[MAXSTR];
  int i,j,*id;
  double *a,*b,**c;
  FILE *fp;
  
  initErrorMsg();

  a  = (double *) malloc(sizeof(double)*NP);
  b  = (double *) malloc(sizeof(double)*NP);
  c  = (double **) malloc(sizeof(double*)*NP);
  id = (int *) malloc(sizeof(int)*NP);

  for (i=0; i<NP; i++) 
    c[i] = (double *) malloc(sizeof(double)*NQ);

  if ((fp = fopen("tarray.dat","r")) == NULL)
    writeError(errModeAll,errFatal,"Data file tarray.dat not found\n");
  fgets(txt,MAXSTR,fp);
  for (i=0;i<NP;i++) fscanf(fp,"%lf",&a[i]);
  fclose(fp);
  printf("\noriginal array:\n");
  for (i=0;i<=9;i++) {
    for (j=0;j<=9;j++) printf("%7.2lf",a[10*i+j]);
    printf("\n");
  }
  sort1(a,NP);
  printf("\nsorted array:\n");
  for (i=0;i<=9;i++) {
    for (j=0;j<=9;j++) printf("%7.2lf",a[10*i+j]);
    printf("\n");
  }

  if ((fp = fopen("tarray.dat","r")) == NULL)
    writeError(errModeAll,errFatal,"Data file tarray.dat not found\n");
  fgets(txt,MAXSTR,fp);
  for (i=0;i<NP/2;i++) fscanf(fp,"%lf %lf",&a[i],&b[i]);
  fclose(fp);
  printf("\noriginal array:\n");
  for (i=0;i<=9;i++) {
    for (j=0;j<=3;j++) printf("[%7.2lf,%7.2lf]",a[5*i+j],b[5*i+j]);
    printf("\n");
  }
  sort2(a,b,NP/2);
  printf("\nsorted array:\n");
  for (i=0;i<=9;i++) {
    for (j=0;j<=3;j++) printf("[%7.2lf,%7.2lf]",a[5*i+j],b[5*i+j]);
    printf("\n");
  }


  if ((fp = fopen("tarray.dat","r")) == NULL)
    writeError(errModeAll,errFatal,"Data file tarray.dat not found\n");
  fgets(txt,MAXSTR,fp);
  for (i=0;i<NP/2;i++) {
    fscanf(fp,"%lf %lf",&a[i],&b[i]);
    b[i]  = (double) i;
    id[i] = i;
  }
  fclose(fp);
  printf("\noriginal array:\n");
  for (i=0;i<=9;i++) {
    for (j=0;j<=3;j++) printf("[%7.2lf,%3d]",a[5*i+j],id[5*i+j]);
    printf("\n");
  }
  isort2(a,id,NP/2);
  printf("\nsorted array 1:\n");
  for (i=0;i<=9;i++) {
    for (j=0;j<=3;j++) printf("[%7.2lf,%7.2lf]",a[5*i+j],b[5*i+j]);
    printf("\n");
  }
  printf("\nsorted array 2:\n");
  for (i=0;i<=9;i++) {
    for (j=0;j<=3;j++) printf("[%7.2lf,%3d]",a[5*i+j],id[5*i+j]);
    printf("\n");
  }


  if ((fp = fopen("tarray.dat","r")) == NULL)
    writeError(errModeAll,errFatal,"Data file tarray.dat not found\n");
  fgets(txt,MAXSTR,fp);
  for (i=0;i<NP/5;i++) {
    fscanf(fp,"%lf",&a[i]);
    for (j=0; j<4; j++) fscanf(fp,"%lf",&c[i][j]);
  }
  fclose(fp);
  printf("\noriginal array:\n");
  for (i=0;i<=19;i++) {
    printf("[%7.2lf",a[i]);
    for (j=0;j<4;j++) printf(",%7.2lf",c[i][j]);
    printf("]\n");
  }
  sortN(a,c,NP/5);
  printf("\nsorted array:\n");
  for (i=0;i<=19;i++) {
    printf("[%7.2lf",a[i]);
    for (j=0;j<4;j++) printf(",%7.2lf",c[i][j]);
    printf("]\n");
  }


  free(a);
  free(b);
  for (i=0; i<NP; i++) free(c[i]);
  free(c);
  return 0;
}
