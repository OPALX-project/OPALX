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
void AnalyzeL2D(int Nstep=0){
  // load the shared library:
  gSystem->Load("/home5/schietinger/bin/libh5root-yang.so");

  // let's do it in style:
  TH5Style::SetStyle();
  // get rid of statistic box
  gStyle->SetOptStat(0);
  
  // load the data file:
  TH5Dataset *data1 = new TH5Dataset("Scenario-3/Ring-3.h5");
  TH5Dataset *data2 = new TH5Dataset("Scenario-7/Ring-3.h5");
  TH5Dataset *data3 = new TH5Dataset("Scenario-6/Ring-3.h5");

  // for convenience, get a dummy canvas for plotting
  TCanvas* dummy = new TCanvas("dummy","dummy",10, 10);      

  // get a 2D x-y histogram for step 100:
  //  TH2D* histoZPZ1 = data1.Histo2d("z","pz",100,-0.008,0.008,-0.002, 0.002);

  double XminL=-0.07*1000.0;
  double XmaxL=0.07*1000.0;
  double YminL=-0.07*1000.0;
  double YmaxL=0.07*1000.0;

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

  TH2D* histoXY1 = new TH2D("histoXY1","top view in local frame",50,XminL,XmaxL, 50,YminL,YmaxL); 
  ntd->Draw(Form("1000.0*((x-%g)*%g+(y-%g)*%g):1000.0*((-x+%g)*%g+(y-%g)*%g)>>histoXY1",cx[Nstep],cosAngle[0],cy[Nstep],sinAngle[0],cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));
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

  TH2D* histoXY2 = new TH2D("histoXY2","top view in local frame",50,XminL,XmaxL,50,YminL,YmaxL); 
  ntd->Draw(Form("1000.0*((x-%g)*%g+(y-%g)*%g):1000.0*((-x+%g)*%g+(y-%g)*%g)>>histoXY2",cx[Nstep],cosAngle[0],cy[Nstep],sinAngle[0],cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));
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

  TH2D* histoXY3 = new TH2D("histoXY3","top view in local frame",50,XminL,XmaxL,50,YminL,YmaxL); 
  ntd->Draw(Form("1000.0*((x-%g)*%g+(y-%g)*%g):1000.0*((-x+%g)*%g+(y-%g)*%g)>>histoXY3",cx[Nstep],cosAngle[0],cy[Nstep],sinAngle[0],cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));

  delete ntd;
  
  ///////////////

  delete dummy;

  // embellish the histograms at wish and draw it with the desired 
  // drawing option:
  histoXY1->SetLabelSize(0.07,"X");
  histoXY1->SetLabelSize(0.07,"Y");
  histoXY1->SetTitleSize(0.07,"X");
  histoXY1->SetTitleSize(0.07,"Y");
  histoXY1->SetXTitle("longitudinal [mm] ");
  histoXY1->SetYTitle("radial [mm] ");
  histoXY1->Scale(1);
  
  histoXY2->SetLabelSize(0.07,"X");
  histoXY2->SetLabelSize(0.07,"Y");
  histoXY2->SetTitleSize(0.07,"X");
  histoXY2->SetTitleSize(0.07,"Y");
  histoXY2->SetXTitle("longitudinal [mm] ");
  histoXY2->SetYTitle("radial [mm] ");
  histoXY2->Scale(1);
  
  histoXY3->SetLabelSize(0.07,"X");
  histoXY3->SetLabelSize(0.07,"Y");
  histoXY3->SetTitleSize(0.07,"X");
  histoXY3->SetTitleSize(0.07,"Y");
  histoXY3->SetXTitle("longitudinal [mm] ");
  histoXY3->SetYTitle("radial [mm] ");
  histoXY2->Scale(1);

  // get a fresh canvas, tune it as you like
  TCanvas* myCanv = new TCanvas("myCanv","H5PartROOT",1200, 400);  
  myCanv->Divide(3);
  myCanv->cd(1);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.2);

  TLegend *legend3=new TLegend(0.65,0.85,0.95,.95);
  legend3->SetTextSize(0.07);
  legend3->SetFillColor(0);   
  legend3->SetHeader("1 beam");
  
  histoXY1->Draw("cont1");
  legend3->Draw();
  //  delete legend3;
  
  myCanv->cd(2);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.2);

  TLegend *legend4=new TLegend(0.65,0.85,0.95,0.95);
  legend4->SetTextSize(0.07);
  legend4->SetFillColor(0);   
  legend4->SetHeader("7 beams");
  
  histoXY2->Draw("cont1");
  legend4->Draw();
  //  delete legend4;
  
  myCanv->cd(3);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.2);
  

  TLegend *legend5=new TLegend(0.65,0.85,0.95,0.95);
  legend5->SetTextSize(0.07);
  legend5->SetFillColor(0);   
  legend5->SetHeader("9 beams");

  histoXY3->Draw("cont1");
  legend5->Draw();
  //delete legend5;
  
  myCanv->Update();

  TString s1 = TString("Turn-")+Nstep+".pdf";
  TString s2 = TString("Turn-")+Nstep+".png";
  
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
