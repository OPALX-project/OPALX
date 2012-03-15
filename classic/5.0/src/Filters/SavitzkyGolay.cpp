#include "Filters/SavitzkyGolay.h"
#include "fftw3.h"
#include <fstream>

SavitzkyGolayFilter::SavitzkyGolayFilter(int np, int nl, int nr, int m):
  NumberPoints(np),
  NumberPointsLeft(nl),
  NumberPointsRight(nr),
  PolynomialOrder(m),
  Coefs(np,0.0),
  CoefsDeriv(np,0.0)
{
  savgol(Coefs, NumberPoints, NumberPointsLeft, NumberPointsRight, 0, PolynomialOrder);
  savgol(CoefsDeriv, NumberPoints, NumberPointsLeft, NumberPointsRight, 1, PolynomialOrder);
}

void SavitzkyGolayFilter::apply(vector<double> &LineDensity) {
  vector<double> temp(LineDensity.size(), 0.0);
  convlv(LineDensity, Coefs, 1, temp);
  LineDensity.assign(temp.begin(), temp.end());
}

void SavitzkyGolayFilter::calc_derivative(vector<double> &LineDensity, const double &h) {
  vector<double> temp(LineDensity.size(), 0.0);
  convlv(LineDensity, CoefsDeriv, 1, temp);
  LineDensity.assign(temp.begin(), temp.end());
}


void savgol(vector<double> &c, const int &np, const int &nl, const int &nr, const int &ld, const int &m)
{
  int j, k, imj, ipj, kk, mm;
  double d, fac, sum;
  
  if (np < nl + nr + 1 || nl < 0 || nr < 0 || ld > m || nl + nr < m)
    {
      cerr << "bad args in savgol" << endl;
      return;
    }
  vector<int> indx(m+1,0);
  vector<double> a((m+1)*(m+1),0.0);
  vector<double> b(m+1,0.0);

  for (ipj = 0; ipj <= (m << 1); ++ipj)
    {
      sum = (ipj ? 0.0 : 1.0);
      for (k = 1; k <= nr; ++k)
        sum += (int)pow((double)k, (double)ipj);
      for (k = 1; k <= nl; ++k) 
        sum += (int)pow((double)-k, (double)ipj);
      mm = (ipj < 2 * m - ipj ? ipj: 2 * m - ipj);
      for (imj = -mm; imj <= mm; imj += 2)
        a[(ipj + imj)/2 * (m + 1) + (ipj - imj)/2] = sum;
    }
  ludcmp(a, indx, d);

  for (j = 0; j < m + 1; ++j)
    b[j] = 0.0;
  b[ld] = 1.0;

  lubksb(a, indx, b);
  for (kk = 0; kk < np; ++kk)
    c[kk] = 0.0;
  for (k = -nl; k <= nr; ++k)
    {
      sum = b[0];
      fac = 1.0;
      for (mm = 1; mm <= m; ++mm)
        sum += b[mm] * (fac *= k);
      kk = (np - k) % np;
      c[kk] = sum;
    }

}

void convlv(const vector<double> &data, const vector<double> &respns, const int &isign, vector<double> &ans)
{
  int i;//, no2;
  double mag2, tmp;
  
  int n = data.size();
  int m = respns.size();

  double *temprd = (double*) fftw_malloc(sizeof(double) * n);
  fftw_complex *tempfd1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (1 + (int)floor(n/2)));
  fftw_complex *tempfd2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (1 + (int)floor(n/2)));
  fftw_plan p = fftw_plan_dft_r2c_1d(n, temprd, tempfd1, FFTW_ESTIMATE);
  fftw_plan q = fftw_plan_dft_r2c_1d(n, temprd, tempfd2, FFTW_ESTIMATE);
  fftw_plan r = fftw_plan_dft_c2r_1d(n, tempfd1, temprd, FFTW_ESTIMATE);

  temprd[0] = respns[0];
  for (i = 1; i < (m + 1) / 2; ++i) 
    {
      temprd[i] = respns[i];
      temprd[n - i] = respns[m - i];
    }
  for (i = (m + 1) / 2; i < n - (m - 1) / 2; ++i)
    temprd[i] = 0.0;
  
  fftw_execute(q);
  
  for (i = 0; i < n; ++i)
    temprd[i] = data[i];

  fftw_execute(p);

  //  no2 = n>>1;
  if (isign == 1)
    {
      for (i = 1; i < n/2; ++i)
        {
          tmp = tempfd1[i][0];
          tempfd1[i][0] = (tempfd1[i][0] * tempfd2[i][0] - tempfd1[i][1] * tempfd2[i][1]) / n;//o2;
          tempfd1[i][1] = (tempfd1[i][1] * tempfd2[i][0] + tmp           * tempfd2[i][1]) / n;//o2;
        }
      tempfd1[0][0] = tempfd1[0][0] * tempfd2[0][0] / n;//o2;
      tempfd1[0][1] = 0.0; // due to the way fftw works
      tempfd1[n/2][0] = tempfd1[n/2][0] * tempfd2[n/2][0] / n;//o2;
      tempfd1[n/2][1] = 0.0; // due to the way fftw works
    }
//   else if (isign == -1)
//     {
//       for (i = 1; i < n/2; ++i)
//         {
//           if ((mag2 = tempfd2[i][0] * tempfd2[i][0] + tempfd2[i][1] * tempfd2[i][1]) == 0.)
//             {
//               cerr << "Deconvolving at response zero in convlv" << endl;
//               return;
//             }
//           tmp = tempfd1[i][0];
//           tempfd1[i][0] = (tempfd1[i][0] * tempfd2[i][0] + tempfd1[i][1] * tempfd2[i][1]) / mag2 / n;//o2;
//           tempfd1[i][1] = (tempfd1[i][1] * tempfd2[i][0] - tmp           * tempfd2[i][1]) / mag2 / n;//o2;
//         }
//       if (tempfd2[0][0] == 0.0 || tempfd2[0][1] == 0.0)
//         {
//           cerr << "Deconvolving at response zero in convlv" << endl;
//           return;
//         }
//       tempfd1[0][0] = tempfd1[0][0] / tempfd2[0][0] / n;//o2;
//       tempfd1[0][1] = tempfd1[0][1] / tempfd2[0][1] / n;//o2;
//     }
//   else
//     {
//       cerr << "No meaning for isign in convlv" << endl;
//       return;
//     }
  
  fftw_execute(r);

  for (i = 0; i < n; ++i)
    ans[i] = temprd[i];

  fftw_destroy_plan(p);
  fftw_destroy_plan(q);
  fftw_destroy_plan(r);
  fftw_free(tempfd1);
  fftw_free(tempfd2);
  fftw_free(temprd);
}

void ludcmp(vector<double> &a, vector<int> &indx,  double &d)
{
  const double TINY = 1.0e-20;
  int i, imax, j, k;
  double big, dum, sum, temp;

  int n = indx.size();
  vector<double> vv(n,0.0);

  d = 1.0;
  for (i = 0; i < n; ++i)
    {
      big = 0.0;
      for (j = 0; j < n; ++j)
        if ((temp = fabs(a[i * n + j])) > big) big = temp;
      
      if (big == 0.0)
        {
          cerr << "Singular matrix in routine ludcmp" << endl;
          return;
        }
      vv[i] = 1./big;
    }
  
  for (j = 0; j < n; ++j)
    {
      for (i = 0; i < j; ++i)
        {
          sum = a[i * n + j];
          for (k = 0; k < i; ++k)
            sum -= a[i * n + k] * a[k * n + j];
          a[i * n + j] = sum;
        }
      big = 0.0;
      for (i = j; i < n; ++i)
        {
          sum = a[i * n + j];
          for (k = 0; k <j; ++k)
            sum -= a[i * n + k] * a[k * n + j];
          a[i * n + j] = sum;
          if ((dum = vv[i] * fabs(sum)) >= big)
            {
              big = dum;
              imax = i;
            }
        }
      
      if (j != imax)
        {
          for (k = 0; k < n; ++k)
            {
              dum = a[imax * n + k];
              a[imax * n + k] = a[j * n + k];
              a[j * n + k] = dum;
            }
          d = -d;
          vv[imax] = vv[j];
        }
      indx[j] = imax;
      if (a[j * n + j] == 0.0) a[j * n + j] = TINY;
      if (j != n -1)
        {
          dum = 1./a[j * n + j];
          for (i = j + 1; i < n; ++i)
            a[i * n + j] *= dum;
        }
    }
}

void lubksb(vector<double> &a, vector<int> &indx, vector<double> &b)
{
  int i, ii = 0, ip, j;
  double sum;
  int n = indx.size();

  for (i = 0; i < n; ++i)
    {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii != 0)
        for (j = ii - 1; j < i; ++j)
          sum -= a[i * n + j] * b[j];
      else if (sum != 0.0)
        ii = i + 1;
      b[i] = sum;
    }
  for (i = n - 1; i >= 0; --i)
    {
      sum = b[i];
      for (j = i + 1; j < n; ++j)
        sum -= a[i * n + j] * b[j];
      b[i] = sum/a[i * n + i];
    }

}
