#include "Riostream.h"

const Double_t m0 = 0.511*1.0E3;

Double_t eV2bg (Double_t Eev) {

  //  Double_t tmp = 1.0 + Eev/m0;
  // return sqrt(tmp*tmp + 1.0);
  return Eev;
}


void readDist() {
  ifstream ino, ina;
  ino.open("dist.dat");
  ina.open("Q200pC_th0p249um_XYrms275um_Tfwhm9p9ps_rf0p7ps_200k.ini");
  
  Double_t x,y,t,px,py,pt,d1,d2,d3,d4;
  
  TNtuple *nto =new TNtuple("nto","OPAL data","x:y:t:px:py:pt");
  TNtuple *nta =new TNtuple("nta","Astra data","x:y:t:px:py:pt");
  TNtuple *tmp =new TNtuple("tmp","Astra data","x:y:t:px:py:pt");
  
  Int_t nlines =0;
  while (1) {
    ino >> x >> y >> t >> px >> py >> pt;
    if (!ino.good()) break;
    t += 0.5E-12;
    nto->Fill(x,y,t,px,py,pt);
    nlines++;
  }
  printf(" found %d OPAL points \n",nlines);
  
  nlines = 0;
  while (1) {
    ina >> x >> y >> d1 >> px >> py >> pt >> t >> d2  >> d3 >> d4;
    if (!ina.good()) break;
    t += 1.0E-12;
    t *= 1E-9;
    tmp->Fill(x,y,t,px,py,pt);
    nlines++;
  }
  
  Double_t mint = tmp->GetMinimum("t");
  ina.close();

  ina.open("Q200pC_th0p249um_XYrms275um_Tfwhm9p9ps_rf0p7ps_200k.ini");
  nlines = 0;
  while (1) {
    ina >> x >> y >> d1 >> px >> py >> pt >> t >> d2  >> d3 >> d4;
    if (!ina.good()) break;
    t *= 1E-9;
    t  = t + (-1.0*mint);
    Double_t pxbg = eV2bg(px);
    Double_t pybg = eV2bg(py);
    Double_t ptbg = eV2bg(pt);
    nta->Fill(x,y,t,pxbg,pybg,ptbg);
    nlines++;
  }
  mint = nta->GetMinimum("t");
  printf(" found %d Astra points \n",nlines);
  ina.close();

 
  ino.close(); 

}


