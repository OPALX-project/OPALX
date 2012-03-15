//-*-C++-*-/////////////////////////////////////////////////////////////
// $Id: Basedata.h,v 1.2 2009/06/03 19:07:28 pope Exp $
////////////////////////////////////////////////////////////////////////

#ifndef BASEDATA_H
#define BASEDATA_H

////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////

#ifdef INITIALIZER_NAMESPACE
namespace initializer {
#endif

class Basedata {

public:

  explicit Basedata( const char *file ) :
    m_skip(-1),
    m_jinit(-1),
    m_icbc(-1),
    m_norm(-1),
    m_trans(-1),
    m_initp(-1),
    m_norder(-1),
    m_frm(-1),
    m_pbchoice(-1),
    m_entest(-1),
    m_ng(-1),
    m_ng2d(-1),
    m_ns(-1.0),
    m_iseed(-1),
    m_alpha(-1.0),
    m_zin(-1.0),
    m_zfin(-1.0),
    m_zpr(-1.0),
    m_np(-1),
    m_nsteps(-1),
    m_iprint(-1),
    m_irun(-1),
    m_hubble(-1.0),
    m_omegadm(-1.0),
    m_deut(-1.0),
    m_rL(-1.0),
    m_ss8(-1.0),
    m_qq(-1.0),
    m_nxint(-1),
    m_nxdint(-1),
    m_hdfform(-1),
    m_rprec(-1),
    m_iprec(-1),
    m_byteorder(-1),
    m_oL(-1.0),
    m_Np(-1),
    m_omegatot(-1.0),
    m_pmass(-1.0),
    m_ain(-1.0),
    m_afin(-1.0),
    m_apr(-1.0),
    m_pp(-1.0),
    m_pfin(-1.0),
    m_ppr(-1.0),
    m_adot(-1.0),
    m_prefactor(-1.0)
  {
    std::ifstream myfile;

    myfile.open(file);
    if (! myfile.is_open()) {
      std::ostringstream ost;
      ost << "basedata: cannot open '"
	  << file
	  << "'";
      throw std::runtime_error( ost.str() ); 
    }

    read_item( myfile, &m_skip );
    read_item( myfile, &m_jinit );
    read_item( myfile, &m_icbc );
    read_item( myfile, &m_norm );
    read_item( myfile, &m_trans );
    read_item( myfile, &m_initp );
    read_item( myfile, &m_norder );
    read_item( myfile, &m_frm );
    read_item( myfile, &m_pbchoice );
    read_item( myfile, &m_entest );
    read_item( myfile, &m_ng );
    read_item( myfile, &m_ng2d );
    read_item( myfile, &m_ns );
    read_item( myfile, &m_iseed );
    read_item( myfile, &m_alpha );
    read_item( myfile, &m_zin );
    read_item( myfile, &m_zfin );
    read_item( myfile, &m_zpr );
    read_item( myfile, &m_np );
    read_item( myfile, &m_nsteps );
    read_item( myfile, &m_iprint );
    read_item( myfile, &m_irun );
    read_item( myfile, &m_hubble );
    read_item( myfile, &m_omegadm );
    read_item( myfile, &m_deut );
    read_item( myfile, &m_rL );
    read_item( myfile, &m_ss8 );
    read_item( myfile, &m_qq );
    read_item( myfile, &m_nxint );
    read_item( myfile, &m_nxdint );
    //added by Adrian to match sample indat
    read_item( myfile, &m_hdfform );
    read_item( myfile, &m_rprec );
    read_item( myfile, &m_iprec );
    read_item( myfile, &m_byteorder );
    //added by Adrian for new parameters
    //read_item( myfile, &m_oL );

    m_Np = m_np * m_np * m_np;
    m_omegatot = m_omegadm + m_deut / m_hubble / m_hubble;
    m_pmass = 2.77536627e11*(m_rL*m_rL*m_rL)*(m_omegatot)/m_hubble/m_Np;

    m_ain = 1.0 / (1 + m_zin);
    m_afin = 1.0 / (1 + m_zfin);
    m_apr = 1.0 / (1 + m_zpr);
    m_pp = pow(m_ain, m_alpha);
    m_pfin = pow(m_afin, m_alpha);
    m_ppr = pow(m_apr, m_alpha);
    // m_prefactor = 1.0 / (m_alpha*m_adot*powf(m_pp,1.0+(1.0/m_alpha)));
  }

  float adot_func(float local_pp, float local_alpha) const {
    const float pp1 = powf(local_pp, (3.0 / local_alpha));
    float tmp = m_omegatot + (1 - m_omegatot) * pp1;
    tmp = tmp / (powf(local_pp, (1.0 / local_alpha)));
    return sqrtf(tmp);
  }

  int skip() const { return m_skip; }
  int jinit() const { return m_jinit; } 
  int icbc() const { return m_icbc; } 
  int norm() const { return m_norm; } 
  int trans() const { return m_trans; } 
  int initp() const { return m_initp; } 
  int norder() const { return m_norder; } 
  int frm() const { return m_frm; } 
  int pbchoice() const { return m_pbchoice; } 
  int entest() const { return m_entest; }
  int ng() const { return m_ng; } 
  int ng2d() const { return m_ng2d; } 
  float ns() const { return m_ns; } 
  int iseed() const { return m_iseed; } 
  float alpha() const { return m_alpha; } 
  float zin() const { return m_zin; } 
  float zfin() const { return m_zfin; } 
  float zpr() const { return m_zpr; } 
  int np() const { return m_np; } 
  int nsteps() const { return m_nsteps; } 
  int iprint() const { return m_iprint; } 
  int irun() const { return m_irun; } 
  float hubble() const { return m_hubble; } 
  float omegadm() const { return m_omegadm; } 
  float deut() const { return m_deut; } 
  float rL() const { return m_rL; } 
  float ss8() const { return m_ss8; } 
  float qq() const { return m_qq; } 
  int nxint() const { return m_nxint; } 
  int nxdint() const { return m_nxdint; } 
  int hdfform() const { return m_hdfform; }
  int rprec() const { return m_rprec; }
  int iprec() const { return m_iprec; }
  int byteorder() const { return m_byteorder; }
  float oL() const { return m_oL; }

  int Np() const { return m_Np; }

  float omegatot() const { return m_omegatot; } 
  float pmass() const { return m_pmass; } 

  float ain() const { return m_ain; } 
  float afin() const { return m_afin; } 
  float apr() const { return m_apr; } 
  float pp() const { return m_pp; } 
  float pfin() const { return m_pfin; } 
  float ppr() const { return m_ppr; } 
  float adot() const { return m_adot; } 
  float prefactor() const { return m_prefactor; } 

private:

  Basedata();
  Basedata( const Basedata& );
  Basedata& operator = ( const Basedata& );

  template< typename U > void 
  read_item( std::ifstream& infile, U* item )
  {
    infile >> *item;
    std::string tmp_str;
    getline( infile, tmp_str );
    //if ((! infile) && (! infile.eof())) {
    if ( infile.eof() ) {
      throw std::runtime_error( "Basedata: error reading input file" );
    }
  }

  int m_skip;
  int m_jinit;  
  int m_icbc;
  int m_norm;
  int m_trans;
  int m_initp;
  int m_norder;
  int m_frm;
  int m_pbchoice;
  int m_entest;
  int m_ng;
  int m_ng2d;
  float m_ns;
  int m_iseed;
  float m_alpha;
  float m_zin;
  float m_zfin;
  float m_zpr;
  int m_np;
  int m_nsteps;
  int m_iprint;
  int m_irun;
  float m_hubble;
  float m_omegadm;
  float m_deut;
  float m_rL;
  float m_ss8;
  float m_qq;
  int m_nxint;
  int m_nxdint;
  int m_hdfform;
  int m_rprec;
  int m_iprec;
  int m_byteorder;
  float m_oL;

  int m_Np;

  float m_omegatot;
  float m_pmass;

  float m_ain;
  float m_afin;
  float m_apr;
  float m_pp;
  float m_pfin;
  float m_ppr;
  float m_adot;
  float m_prefactor;
};

#ifdef INITIALIZER_NAMESPACE
}
#endif

////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////
