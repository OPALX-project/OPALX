#include "Solvers/CSRWakeFunction.hh"
#include "Physics/Physics.h" 
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/SBend.h"


CSRWakeFunction::CSRWakeFunction(PartData &ref, ElementBase* element, const S_G_FilterOptions &options1, const S_G_FilterOptions &options2):
  WakeFunction(ref),
  smoother1_m(options1.np_m, options1.nl_m, options1.nr_m, options1.m_m),
  smoother2_m(options2.np_m, options2.nl_m, options2.nr_m, options2.m_m),
  lineDensity_m(),
  dlineDensitydz_m(),
  d2lineDensitydz2_m()
{
  if (dynamic_cast<RBend*>(element))
    {
      RBend* bend = dynamic_cast<RBend*>(element);
//       R_m = bend->getR();
//       Begin_m = bend->getBegin();
//       Length_m = bend->getLength();
    }
  else if (dynamic_cast<SBend*>(element))
    {
      SBend* bend = dynamic_cast<SBend*>(element);
//       R_m = bend->getR();
//       Begin_m = bend->getBegin();
//       Length_m = bend->getLength();
    }
}

void CSRWakeFunction::apply(PartBunch &bunch)
{
  //#ifdef CALCCSR
  double Phi, Phi_m, Psi, SlippageLength;
  double min_relative_s, Ds_max, Ds;
  double frac, x;
  const double R = 1.0; //itsBend_m.getR();
  const double hz = bunch.getMesh().get_meshSpacing(2);
  bunch.calcLineDensity();
  bunch.getLineDensity(lineDensity_m);
  const unsigned int N = lineDensity_m.size();
  Vector_t smin, smax;

  smoother1_m.apply(lineDensity_m);
  lineDensity_m.getFirstDerivative(dlineDensitydz_m, hz);
  
  smoother2_m.apply(dlineDensitydz_m);
  if (d2lineDensitydz2_m.size() != N)
    {
      d2lineDensitydz2_m.resize(N, 0.0);
      Ez_m.resize(N, 0.0);
    }
  IlyaPogorelovFilter::apply(dlineDensitydz_m, d2lineDensitydz2_m, hz);

  Ez_m[0] = 0.0;
  bunch.get_bounds(smin, smax);

  //min_relative_s = smin - Begin_m;
  for (int i = 1; i < N; ++i)
    {
      Phi = (min_relative_s + (i-1)*hz) / R_m;
      SlippageLength = Phi * Phi * Phi * R_m / 24.;
      Ez_m[i] = 0.0;
      if (Phi > 0.0)
        {
          if (Phi < Phi_m)
            {
              if (SlippageLength > hz)
                {
                  Ez_m[i] += 1.5 * pow(hz, 2./3.) * (dlineDensitydz_m[i-1] + 0.5 * d2lineDensitydz2_m[i-1] * hz);
                  if (i > 1)
                    {
                      int j;
                      for (j = i-1; j > 0; ++j)
                        {
                          if ((i - j) * hz > SlippageLength) break;
                          Ez_m[i] += pow((i - j) * hz, -1./3.) * dlineDensitydz_m[j];
                        }
                      Ez_m[i] -= 0.5 * (pow((i - j) * hz, -1./3.) * dlineDensitydz_m[j] * hz 
                                        + pow(hz, 2./3.) * dlineDensitydz_m[i - 1]);
                    }
                }
              if (SlippageLength < (i - 1) * hz)
                {
                  int j = i - static_cast<int>(SlippageLength / hz);
                  frac = (i - j) - SlippageLength / hz;
                  Ez_m[i] += (frac * lineDensity_m[j] + (1. - frac) * lineDensity_m[j + 1]) / pow(SlippageLength, 1./3.);
                }
              if (4. * SlippageLength < (i - 1) * hz)
                {
                  int j = i - 4 * static_cast<int>(SlippageLength / hz);
                  frac = (i - j) - 4. * SlippageLength / hz;
                  Ez_m[i] -= (frac * lineDensity_m[j + 1] + (1. - frac) * lineDensity_m[j]) / pow(SlippageLength, 1./3.);
                }
              Ez_m[i] *= -2. / (4. * Physics::pi * Physics::epsilon_0 * pow(3. * R_m*R_m, 1./3.));
            }
          else
            {
              x = Phi - Phi_m;
              Ds_max = R_m * Phi_m*Phi_m*Phi_m / 24. * (4. - 3.* Phi_m / Phi);
              
              if (Ds_max > hz)
                {
                  Psi = calcPsi(x, hz);
                  Ez_m[i] += (Psi + 2. * x) * Psi*Psi*Psi / ((Psi + x)*(Psi + x)) * dlineDensitydz_m[i-1] * R_m / 16.;

                  if (i > 1)
                    {
                      int j;
                      for (j = i - 1; j > 0; ++j)
                        {
                          Ds = (i - j) * hz;
                          if (Ds > Ds_max) break;
                          Psi = calcPsi(x, Ds);
                          Ez_m[i] += dlineDensitydz_m[j] / (Psi + 2. * x) * hz;
                        }
                      
                      Psi = calcPsi(x, Ds);
                      Ez_m[i] -= 0.5 * dlineDensitydz_m[j] / (Psi + 2. * x) * hz;
                      
                      Psi = calcPsi(x, hz);
                      Ez_m[i] -= 0.5 * dlineDensitydz_m[i-1] / (Psi + 2. * x) * hz;
                    }
                }
              
              if (Ds_max < (i - 1) * hz)
                {
                  int j = i - static_cast<int>(Ds_max / hz);
                  frac = (i - j) - Ds_max / hz;
                  Ez_m[i] += (frac * lineDensity_m[j + 1] + (1. - frac) * lineDensity_m[j]) / (Phi_m + 2. * x);
                }

              Ds_max = R_m * Phi_m*Phi_m * (Phi_m + 3. * x) / 6.;
              if (Ds_max < (i - 1) * hz)
                {
                  int j = i - static_cast<int>(Ds_max / hz);
                  frac = (i - j) - Ds_max / hz;
                  Ez_m[i] += (frac * lineDensity_m[j + 1] + (1. - frac) * lineDensity_m[j]) / (Phi_m + 2. * x);
                }
              Ez_m[i] *= 1. / (R_m * Physics::pi * Physics::epsilon_0);
            }
        }
    }
  
  //#endif
}

double CSRWakeFunction::calcPsi(double x, double Ds)
{
  /** solve the equation 
   *  \f[
   *  \Delta s = \frac{R \Psi^3}{24} \frac{\Psi + 4x}{\Psi + x}
   *  \f]
   *  for \f$\Psi\f$ using Newtons method.
  */
  
  const int Nmax = 100;
  const double eps = 1e-10;
  double residual = 0.0;
  double psi = pow(24. * Ds / R_m, 1./3.);

  for (int i = 0; i < Nmax; ++i)
    {
      residual = R_m * psi*psi*psi*psi + 4. * x * R_m * psi*psi*psi - 24. * Ds * psi - 24. * Ds * x;
      if (fabs(residual) < eps)
        return psi;
      psi -= eps/(4. * R_m * psi*psi*psi + 12. * x * psi*psi - 24. * Ds);
    }
  cerr << "In CSRWakeFunction::calcPsi(): exceed maximum number of iterations!" << endl;
  return psi;
}

const string CSRWakeFunction::getType() const
{
  return "CSRWakeFunction";
}
