/*
Andreas Adelmann & Arno Candel
Paul Scherrer Institut PSI

$Header: /afs/psi.ch/user/a/adelmann/public/cvsroot/ChargedParticles/ListElem.h,v 1.4 2001/08/16 15:22:35 arno Exp $
$Id: ListElem.h,v 1.4 2001/08/16 15:22:35 arno Exp $
$Locker:  $

$Revision: 1.4 $
$State: Exp $
$Date: 2001/08/16 15:22:35 $

$Log: ListElem.h,v $
Revision 1.4  2001/08/16 15:22:35  arno
added ouput of Planar Projections, thus new ListElem needed, with 2 indices

Revision 1.3  2001/08/14 08:30:35  candel
List output only from 1 to N-1, assures centered R plots

Revision 1.2  2001/08/11 05:46:40  adelmann
Add header with CVS tags

add << operator


*/

#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <strstream>
using namespace std;

class ListElem;

typedef list<ListElem>::const_iterator liIt;

struct ListElem {
  double s;
  double t;
  unsigned int m,n;
  double den;
  
  ListElem(double sval, double tval, unsigned int mval, unsigned int nval, double denval) :
    s(sval),
    t(tval),
    m(mval),
    n(nval),
    den(denval)
  {}
  
  ~ListElem() 
  {}
  
  bool operator< (const ListElem& elem) const {
    return ( ( m < elem.m ) && ( n <= elem.n ) );
  }
  
};

template<class ListElem>
ofstream& operator<<(ofstream& os, const list<ListElem> &l)
{
  unsigned int mmax=0;
  unsigned int nmax=0;
  for(liIt it=l.begin();it!=l.end();++it)  {
    if (it->m > mmax) mmax=it->m;
    if (it->n > nmax) nmax=it->n;
  }

  for(liIt it=l.begin();it!=l.end();++it) 
    if ((it->m < mmax)&&(it->n < nmax)) {
      os << it->s << " " << it->t << " " << it->m << " " << it->n << " " << it->den << endl;
    }
  return os;
}
