// Post processing script for OPAL-cycl in the local frame
// Run it in an interactive root session like this (make sure you have files
// called Ring-1.h5 Ring-2.h5 and Ring-3.h, N is the dumping step ):
//    
//   .x AnalyzeL2D.C(Nstep, NB1, NB2)
//
// To get the shared library, you must compile H5PartROOT with 'gmake shlib'.
// Be sure that the underlying H5Part has been compiled with the -fPIC option.
// If not, reconfigure your H5Part like this (bash shell):
//   cd $H5Part
//   ./autogen.sh
//   CFLAGS = -fPIC ./config
//   make
//
// default loading step is 0
void C9B_7B_SB_2D_1mA(int Nstep=0,int NB1=1, int NB2=1){
  // load the shared library:
  gSystem->Load("/home5/schietinger/bin/libh5root-yang.so");

  // let's do it in style:
  TH5Style::SetStyle();
  // get rid of statistic box
  gStyle->SetOptStat(0);
  
  // load the data file:
  TH5Dataset *data1 = new TH5Dataset("Scenario-3/Ring-2.h5"); //1 bunch
  TH5Dataset *data2 = new TH5Dataset("Scenario-7/Ring-2.h5"); //7 bunches
  TH5Dataset *data3 = new TH5Dataset("Scenario-6/Ring-2.h5"); //9 bunches 

  int Nstep_MB1=Nstep+(NB1-1)/2;
  int Nstep_MB2=Nstep+(NB2-1)/2;
  int Nstep_SB=Nstep;

  cout<<"loaded step of single bunch file = "<<Nstep_SB<<endl;
  cout<<"loaded step of multi-bunch  file = "<<Nstep_MB1<<endl;
  cout<<"Bunches Num. of multi-bunch file = "<<NB1<<endl;
  cout<<"loaded step of multi-bunch  file = "<<Nstep_MB2<<endl;
  cout<<"Bunches Num. of multi-bunch file = "<<NB2<<endl;

  
  // for convenience, get a dummy canvas for plotting
  TCanvas* dummy = new TCanvas("dummy","dummy",10, 10);      

  // get a 2D x-y histogram for step 100:
  //  TH2D* histoZPZ1 = data1.Histo2d("z","pz",100,-0.008,0.008,-0.002, 0.002);

  double XminL=-0.007*1000.0;
  double XmaxL=0.007*1000.0;
  double YminL=-0.007*1000.0;
  double YmaxL=0.007*1000.0;

  ///////  
  TGraph* tcx = data1->GetGraphXY(" ","centroid","x");
  TGraph* tcy = data1->GetGraphXY(" ","centroid","y");
  double* cx = tcx->GetY();
  double* cy = tcy->GetY();
  double angle = calculateAngle( cx[Nstep_SB], cy[Nstep_SB]) ;
  double cosAngle[1];
  double sinAngle[1];
   cosAngle[0] = cos(angle);
  sinAngle[0] = sin(angle);

  TNtupleD *ntd = data1->GetNtuple(Nstep_SB);

  TH2D* histoXY1 = new TH2D("histoXY1","top view in local frame",200,XminL,XmaxL, 200,YminL,YmaxL); 
  ntd->Draw(Form("1000.0*((x-%g)*%g+(y-%g)*%g):1000.0*((-x+%g)*%g+(y-%g)*%g)>>histoXY1",cx[Nstep_SB],cosAngle[0],cy[Nstep_SB],sinAngle[0],cx[Nstep_SB],sinAngle[0],cy[Nstep_SB],cosAngle[0]));
  delete ntd;
  
  ////////////
  tcx = data2->GetGraphXY(" ","centroid","x");
  tcy = data2->GetGraphXY(" ","centroid","y");
  cx = tcx->GetY();
  cy = tcy->GetY();
  angle = calculateAngle( cx[Nstep_MB1], cy[Nstep_MB1]) ;
  cosAngle[0] = cos(angle);
  sinAngle[0] = sin(angle);

  ntd = data2->GetNtuple(Nstep_MB1);

  TH2D* histoXY2 = new TH2D("histoXY2","top view in local frame",200,XminL,XmaxL,200,YminL,YmaxL); 
  ntd->Draw(Form("1000.0*((x-%g)*%g+(y-%g)*%g):1000.0*((-x+%g)*%g+(y-%g)*%g)>>histoXY2",cx[Nstep_MB1],cosAngle[0],cy[Nstep_MB1],sinAngle[0],cx[Nstep_MB1],sinAngle[0],cy[Nstep_MB1],cosAngle[0]));
  delete ntd;
  
  ///////////////
  tcx = data3->GetGraphXY(" ","centroid","x");
  tcy = data3->GetGraphXY(" ","centroid","y");
  cx = tcx->GetY();
  cy = tcy->GetY();
  angle = calculateAngle( cx[Nstep_MB2], cy[Nstep_MB2]) ;
  cosAngle[0] = cos(angle);
  sinAngle[0] = sin(angle);
  
  ntd = data3->GetNtuple(Nstep_MB2);

  TH2D* histoXY3 = new TH2D("histoXY3","top view in local frame",200,XminL,XmaxL,200,YminL,YmaxL); 
  ntd->Draw(Form("1000.0*((x-%g)*%g+(y-%g)*%g):1000.0*((-x+%g)*%g+(y-%g)*%g)>>histoXY3",cx[Nstep_MB2],cosAngle[0],cy[Nstep_MB2],sinAngle[0],cx[Nstep_MB2],sinAngle[0],cy[Nstep_MB2],cosAngle[0]));

  // root [46] TH1D* histo = data->Histo1d("sqrt(1+px*px+py*py+pz*pz)", 140, 1.502, 1.513, 1000,"sqrt(1+px*px+py*py+pz*pz)>1.505&&sqrt(1+px*px+py*py+pz*pz)<1.51$
  // root [47] TH1D* histo1 = data1->Histo1d("sqrt(1+px*px+py*py+pz*pz)", 136, 1.502, 1.513, 1000,"");

  //  *	Histo2d(const TString* varNameX, const TString* varNameY, const Int_t step = 0, const Double_t minX = -1., const Double_t maxX = 1., const Double_t minY = -1.,
  //     const Double_t maxY = 1., const Int_t nBinX = 50, const Int_t nBinY = 50, TCut cut = "", const Bool_t keep = kFALSE)
  
  ntd->Print();
  delete ntd;
  
  ///////////////

  delete dummy;

  // embellish the histograms at wish and draw it with the desired 
  // drawing option:
  histoXY1->SetLabelSize(0.07,"X");
  histoXY1->SetLabelSize(0.07,"Y");
  histoXY1->SetTitleSize(0.07,"X");
  histoXY1->SetTitleSize(0.07,"Y");
  histoXY1->SetXTitle("longitudinal (mm) ");
  histoXY1->SetYTitle("transverse (mm) ");
  histoXY1->Scale(1);
  
  histoXY2->SetLabelSize(0.07,"X");
  histoXY2->SetLabelSize(0.07,"Y");
  histoXY2->SetTitleSize(0.07,"X");
  histoXY2->SetTitleSize(0.07,"Y");
  histoXY2->SetXTitle("longitudinal (mm) ");
  histoXY2->SetYTitle("transverse (mm) ");
  histoXY2->Scale(1);
  
  histoXY3->SetLabelSize(0.07,"X");
  histoXY3->SetLabelSize(0.07,"Y");
  histoXY3->SetTitleSize(0.07,"X");
  histoXY3->SetTitleSize(0.07,"Y");
  histoXY3->SetXTitle("longitudinal (mm) ");
  histoXY3->SetYTitle("transverse (mm) ");
  histoXY2->Scale(1);

  // get a fresh canvas, tune it as you like
  TCanvas* myCanv = new TCanvas("myCanv","H5PartROOT",1200, 400);  
  myCanv->Divide(3);
  myCanv->cd(1);
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.2);

  TLegend *legend3=new TLegend(0.5,0.25,0.85,.35);
  legend3->SetTextSize(0.07);
  legend3->SetFillColor(0);   
  legend3->SetHeader("1 bunch");
  
  histoXY1->Draw("contz");
  legend3->Draw();
  //  delete legend3;
  
  myCanv->cd(2);
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.2);

  TLegend *legend4=new TLegend(0.5,0.25,0.85,0.35);
  legend4->SetTextSize(0.07);
  legend4->SetFillColor(0);   
  legend4->SetHeader("7 bunches");
  
  histoXY2->Draw("contz");
  legend4->Draw();
  //  delete legend4;
  
  myCanv->cd(3);
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.2);
  

  TLegend *legend5=new TLegend(0.5,0.25,0.85,0.35);
  legend5->SetTextSize(0.07);
  legend5->SetFillColor(0);   
  legend5->SetHeader("9 bunches");

  histoXY3->Draw("contz");
  legend5->Draw();
  //delete legend5;
  
  myCanv->Update();

  TString s1 = TString("C9B7BSB-2D-1mA-")+Nstep+".pdf";
  TString s2 = TString("C9B7BSB-2D-1mA-")+Nstep+".png";
  
  myCanv->Print(s1);
  myCanv->Print(s2);

}

double calculateAngle(const double x,const double y)
{

  const double pi= 3.14159265358979323846;
  double thetaXY;
  if      ((x>0) && (y>=0)) thetaXY=atan(y/x); 
  else if ((x<0) && (y>=0)) thetaXY=pi+atan(y/x);
  else if ((x<0) && (y<=0)) thetaXY=pi+atan(y/x); 
  else if ((x>0) && (y<=0)) thetaXY=2.0*pi+atan(y/x);
  else if ((x==0) && (y> 0)) thetaXY=pi/2.0;
  else if ((x==0) && (y< 0)) thetaXY=3.0/2.0*pi;
  cout<<"angle = "<<thetaXY/pi*180.0<<endl;
  
  return thetaXY;

}
