#ifndef CLASSIC_PartPusher_H
#define CLASSIC_PartPusher_H

#include "Algorithms/PartBunch.h"
#include "Algorithms/PartData.h"
#include "Physics/Physics.h"
#include "Ippl.h"

class BorisPusher
{
 public:
  BorisPusher(const PartData &ref);
  BorisPusher();
  void initialise(const PartData *ref);
  void kick(const Vector_t &R, Vector_t &P, const Vector_t &Ef, const Vector_t &Bf, const double &dt) const;
  void push(Vector_t &R, const Vector_t &P, const double &dt);
 private:
  const PartData *itsReference;
};

inline BorisPusher::BorisPusher(const PartData &ref):
  itsReference(&ref)
{ }

inline BorisPusher::BorisPusher():
  itsReference(NULL)
{}

inline void BorisPusher::initialise(const PartData *ref)
{ itsReference = ref; }

inline void BorisPusher::kick(const Vector_t &R, Vector_t &P, const Vector_t &Ef, const Vector_t &Bf, const double &dt) const 
{
  using Physics::c;
  Vector_t um, s, a;
  double recpgamma, tmp;
  /** Update the momenta using the \f$D_1\f$ algorithm as described in Birdsall and Langdon's book:
   *
   *  \f[ \vec{v}_{n+1/2} = \frac{1}{2}\overline{\vec{a}_{n}} \; \Delta t + \mathbb{R} \cdot (\vec{v}_{n-1/2} + \frac{1}{2} \overline{\vec{a}_{n}} \; \Delta t) \f]
   *
   *  where the operator \f$\mathbb{R}\f$ effects a rotation through angle \f$-2\tan^{-1}(\vec{\Omega} \; \Delta t/2)\f$ where \f$\vec{\Omega} = Q \vec{B}_{n}/(m_{e} c)\f$.
   *  \f$\mathbb{R}\f$ can be written as 
   *
   *  \f[\mathbb{R} = \frac{(1 - \Theta^2) \; \mathbb{I} + 2 \; \Theta \Theta^{T} - 2 \; \Theta \times \mathbb{I}}{1+\Theta}\f]
   *
   *  where \f$\vec{\Theta} = \vec{\Omega} \; \Delta t/2\f$ and \f$\mathbb{I}\f$ is the unit tensor.
   */

  /** \f[ \vec{u}_{m} = \vec{\beta}_{n} \gamma_{n} = \vec{\beta}_{n-1/2}\; \gamma_{n-1/2} + \frac{q}{m_{e} c} \vec{E} \frac{\Delta t}{2} \f]
   * \code
   * um = P + 0.5 * Q * dt/M * c * E; 
   * \endcode
   */
  um = P + 0.5 * itsReference->getQ() * dt / itsReference->getM() * c * Ef;

  /** \f[ \vec{a}_{n} = \frac{Q \; \Delta t}{2\gamma M c^2}\vec{B} \f]
   * \code
   * recpgamma = 1.0 / sqrt(1.0 + um[0]*um[0] + um[1]*um[1] + um[2]*um[2]); 
   * 
   * tmp = 0.5 * Q * recpgamma * dt/M  * c * c;
   * a = tmp * B;
   * \endcode
   */
  recpgamma = 1.0 / sqrt(1.0 + dot(um,um));
      
  tmp = 0.5 * itsReference->getQ() * recpgamma * dt / itsReference->getM() * c * c;
  a = tmp * Bf;

  /** \f[ \vec{s} = \vec{\beta}_{n} \gamma_{n} + \frac{Q \; \Delta t}{2\gamma_{n} M c^2} \cdot \gamma_{n} \vec{\beta}_{n} \wedge \vec{B}_{n} \f]
   * \code
   * s = um + tmp * cross(um, B);
   * \endcode
   */
  s = um + tmp * cross(um,Bf);

  /** \f[\vec{u}_{m} = \frac{\mathbb{I} \vec{s} + \vec{a}_{n}\vec{a}_{n}^{T} \vec{s} - \vec{a}_{n} \wedge \vec{s}}{1 + a_{n}^{2}} \f]
   *
   * \code
   * tmp = 1.0 + dot(a,a);
   * um(0) = ((1.0 + a(0)*a(0))    * s(0) + (a(0) * a(1) + a(2)) * s(1) + (a(0) * a(2) - a(1)) * s(2)) / tmp;
   * um(1) = ((a(0) * a(1) - a(2)) * s(0) +    (1.0 + a(1)*a(1)) * s(1) + (a(1) * a(2) + a(0)) * s(2)) / tmp;
   * um(2) = ((a(0) * a(2) + a(1)) * s(0) + (a(1) * a(2) - a(0)) * s(1) +    (1.0 + a(2)*a(2)) * s(2)) / tmp;
   * \endcode
   */
  tmp = 1.0 + dot(a,a);

  // since there are no matrices and therefore the multiplication of a vector with the transpose of some other
  // vector does not exist we have to do the next step component wise...
  um(0) = ((1.0 + a(0)*a(0))    * s(0) + (a(0) * a(1) + a(2)) * s(1) + (a(0) * a(2) - a(1)) * s(2)) / tmp;
  um(1) = ((a(0) * a(1) - a(2)) * s(0) +    (1.0 + a(1)*a(1)) * s(1) + (a(1) * a(2) + a(0)) * s(2)) / tmp;
  um(2) = ((a(0) * a(2) + a(1)) * s(0) + (a(1) * a(2) - a(0)) * s(1) +    (1.0 + a(2)*a(2)) * s(2)) / tmp;

  /** \f[ \vec{\beta}_{n+1/2}\; \gamma_{n+1/2}  = \frac{\mathbb{I} \vec{s} + \vec{a}_{n}\vec{a}_{n}^{T} \vec{s} - \vec{a}_{n} \wedge \vec{s}}{1 + a_{n}^{2}} \\
   *                                            = \frac{(1 - (\frac{q}{m_{e} \gamma_{n}} \frac{\Delta t}{2} \vec{B}_{n})^{2}) \; \mathbb{I} + 2 \; \frac{q^{2}}{m_{e}^{2} \; \gamma_{n}^2}\frac{(\Delta t)^{2}}{4} \vec{B} \vec{B}^{T} - 2 \; \frac{q}{m_{e} \gamma_{n}} \frac{\Delta t}{2}\; \vec{B} \wedge \mathbb{I}}{1 + (\frac{q}{m_{e} \gamma_{n}} \frac{\Delta t}{2} \vec{B})^{2}}\; \vec{\beta}_{n} \gamma_{n} \\
   *                                            = \mathbb{R} \vec{\beta}_{n}\; \gamma_{n},\f]  
   * \code
   * P = um + 0.5 * Q * dt / M * c * E; 
   * \endcode
   */
  P = um + 0.5 * itsReference->getQ()  * dt / itsReference->getM() * c * Ef; 

}

inline void BorisPusher::push(Vector_t &R, const Vector_t &P, const double &dt)
{
  /** \f[ \vec{x}_{n+1/2} = \vec{x}_{n} + \frac{1}{2}\vec{v}_{n-1/2}\quad (= \vec{x}_{n} + \frac{\Delta t}{2} \frac{\vec{\beta}_{n-1/2}\gamma_{n-1/2}}{\gamma_{n-1/2}}) \f]
   *
   * \code
   * R[i] += 0.5 * P[i] * recpgamma; 
   * \endcode
   */
  R += 0.5 * P / sqrt(1.0 + dot(P,P));
}

#endif
