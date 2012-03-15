// Post processing script for OPAL-cycl in the local frame
// Run it in an interactive root session like this (make sure you have files
// called Ring-1.h5 Ring-2.h5 and Ring-3.h, N is the dumping step ):
//    
//   .x AnalyzeL.C(N)
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
void C2o_6o_10o_Th(int Nstep=0){
  // load the shared library:
  gSystem->Load("/home5/schietinger/extlib/H5Root/lib/libh5root.so");

  // let's do it in style:
  TH5Style::SetStyle();
  // get rid of statistic box
  gStyle->SetOptStat(0);
  
  // load the data file:
  TH5Dataset *data1 = new TH5Dataset("Ring-1/Ring-1.h5");
  TH5Dataset *data2 = new TH5Dataset("Ring-2/Ring-2.h5");
  TH5Dataset *data3 = new TH5Dataset("Ring-3/Ring-3.h5");

  // for convenience, get a dummy canvas for plotting
  TCanvas* dummy = new TCanvas("dummy","dummy",10, 10);      

  // get a 2D x-y histogram for step 100:
  //  TH2D* histoZPZ1 = data1.Histo2d("z","pz",100,-0.008,0.008,-0.002, 0.002);

  double XminL=-80;
  double XmaxL=80;
  double YminL=-80;
  double YmaxL=80;

  ///////  
  TGraph* tcx = data1->GetGraphXY(" ","centroid","x");
  TGraph* tcy = data1->GetGraphXY(" ","centroid","y");
  double* cx = tcx->GetY();
  double* cy = tcy->GetY();
  double angle = calculateAngle( cx[Nstep], cy[Nstep]) ;
  double cosAngle[1];
  double sinAngle[1];
  cosAngle[0] = cos(angle);
  sinAngle[0] = sin(angle);

  TNtupleD *ntd = data1->GetNtuple(Nstep);

  TH1D* h1dX1 = new TH1D("h1dX1","h1 title",256,XminL,XmaxL);
  ntd->Draw(Form("1000.0*((-x+%g)*%g+(y-%g)*%g)>>h1dX1",cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));

  delete ntd;
  
  ////////////
  tcx = data2->GetGraphXY(" ","centroid","x");
  tcy = data2->GetGraphXY(" ","centroid","y");
  cx = tcx->GetY();
  cy = tcy->GetY();
  angle = calculateAngle( cx[Nstep], cy[Nstep]) ;
  cosAngle[0] = cos(angle);
  sinAngle[0] = sin(angle);

  ntd = data2->GetNtuple(Nstep);

  TH1D* h1dX2 = new TH1D("h1dX2","h1 title",256,XminL,XmaxL);
  ntd->Draw(Form("1000.0*((-x+%g)*%g+(y-%g)*%g)>>h1dX2",cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));
  delete ntd;
  
  ///////////////
  tcx = data3->GetGraphXY(" ","centroid","x");
  tcy = data3->GetGraphXY(" ","centroid","y");
  cx = tcx->GetY();
  cy = tcy->GetY();
  angle = calculateAngle( cx[Nstep], cy[Nstep]) ;
  cosAngle[0] = cos(angle);
  sinAngle[0] = sin(angle);
  
  ntd = data3->GetNtuple(Nstep);

  TH1D* h1dX3 = new TH1D("h1dX3","h1 title",256,XminL,XmaxL);
  ntd->Draw(Form("1000.0*((-x+%g)*%g+(y-%g)*%g)>>h1dX3",cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));

  delete ntd;
  
  ///////////////

  delete dummy;
  // get a fresh canvas, tune it as you like
  TCanvas* myCanv = new TCanvas("myCanv","H5PartROOT",1000, 750);  
  myCanv->cd(1);

  TLegend *legend1=new TLegend(0.67,0.67,0.90,0.90);    
  legend1->SetTextSize(0.07);
  legend1->SetFillColor(0);
  // legend1->SetMargin(0.2);

  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.2);

  gPad->SetLogy();

  h1dX1->SetLabelSize(0.07,"X");
  h1dX1->SetLabelSize(0.07,"Y");
  h1dX1->SetTitleSize(0.07,"X");
  h1dX1->SetTitleSize(0.07,"Y");
  h1dX1->SetXTitle("longitudinal (mm)");
  h1dX1->SetYTitle("Density (a.u.)");
  h1dX1->SetFillColor(38);

  //  h1dX2->SetFillColor(42);
  //  h1dX3->SetFillColor(5);

  h1dX1->SetLineColor(1);
  h1dX2->SetLineColor(2);
  h1dX3->SetLineColor(4);
  h1dX2->SetLineWidth(5);
  h1dX3->SetLineWidth(5);

  legend1->AddEntry(h1dX1,"2 deg.","f");
  legend1->AddEntry(h1dX2,"6 deg.","l");
  legend1->AddEntry(h1dX3,"10 deg.","l");
    
  h1dX1->Draw();
  h1dX2->Draw("LSAME");
  h1dX3->Draw("LSAME");
  legend1->Draw();
  //  delete legend1;
  myCanv->Update();

  TString s1 = TString("Theta-Turn-")+Nstep+".pdf";
  TString s2 = TString("Theta-Turn-")+Nstep+".png";
  
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
