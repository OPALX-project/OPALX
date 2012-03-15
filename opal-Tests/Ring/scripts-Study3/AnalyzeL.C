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
void AnalyzeL(int Nstep=0){
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

  double XminL=-0.07;
  double XmaxL=0.07;
  double YminL=-0.07;
  double YmaxL=0.07;

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

  TH2D* histoXY1 = new TH2D("histoXY1","top view in local frame",256,XminL,XmaxL, 256,YminL,YmaxL); 
  ntd->Draw(Form("(x-%g)*%g+(y-%g)*%g:(-x+%g)*%g+(y-%g)*%g>>histoXY1",cx[Nstep],cosAngle[0],cy[Nstep],sinAngle[0],cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));

  TH1D* h1dX1 = new TH1D("h1dX1","h1 title",256,XminL,XmaxL);
  ntd->Draw(Form("(-x+%g)*%g+(y-%g)*%g>>h1dX1",cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));

  TH1D* h1dY1 = new TH1D("h1dY1","h2 title",256,YminL,YmaxL);
  ntd->Draw(Form("(x-%g)*%g+(y-%g)*%g>>h1dY1",cx[Nstep],cosAngle[0],cy[Nstep],sinAngle[0]));
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

  TH2D* histoXY2 = new TH2D("histoXY2","top view in local frame",256,XminL,XmaxL,256,YminL,YmaxL); 
  ntd->Draw(Form("(x-%g)*%g+(y-%g)*%g:(-x+%g)*%g+(y-%g)*%g>>histoXY2",cx[Nstep],cosAngle[0],cy[Nstep],sinAngle[0],cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));

  TH1D* h1dX2 = new TH1D("h1dX2","h1 title",256,XminL,XmaxL);
  ntd->Draw(Form("(-x+%g)*%g+(y-%g)*%g>>h1dX2",cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));

  TH1D* h1dY2 = new TH1D("h1dY2","h2 title",256,YminL,YmaxL);
  ntd->Draw(Form("(x-%g)*%g+(y-%g)*%g>>h1dY2",cx[Nstep],cosAngle[0],cy[Nstep],sinAngle[0]));
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

  TH2D* histoXY3 = new TH2D("histoXY3","top view in local frame",256,XminL,XmaxL,256,YminL,YmaxL); 
  ntd->Draw(Form("(x-%g)*%g+(y-%g)*%g:(-x+%g)*%g+(y-%g)*%g>>histoXY3",cx[Nstep],cosAngle[0],cy[Nstep],sinAngle[0],cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));

  TH1D* h1dX3 = new TH1D("h1dX3","h1 title",256,XminL,XmaxL);
  ntd->Draw(Form("(-x+%g)*%g+(y-%g)*%g>>h1dX3",cx[Nstep],sinAngle[0],cy[Nstep],cosAngle[0]));

  TH1D* h1dY3 = new TH1D("h1dY3","h2 title",256,YminL,YmaxL);
  ntd->Draw(Form("(x-%g)*%g+(y-%g)*%g>>h1dY3",cx[Nstep],cosAngle[0],cy[Nstep],sinAngle[0]));

  delete ntd;
  
  ///////////////

  delete dummy;

  // embellish the histograms at wish and draw it with the desired 
  // drawing option:
  histoXY1->SetLabelSize(0.05,"X");
  histoXY1->SetLabelSize(0.05,"Y");
  histoXY1->SetXTitle("longitudinal [m] ");
  histoXY1->SetYTitle("radial [m] ");
  
  histoXY2->SetLabelSize(0.05,"X");
  histoXY2->SetLabelSize(0.05,"Y");
  histoXY2->SetXTitle("longitudinal [m] ");
  histoXY2->SetYTitle("radial [m] ");

  histoXY3->SetLabelSize(0.05,"X");
  histoXY3->SetLabelSize(0.05,"Y");
  histoXY3->SetXTitle("longitudinal [m] ");
  histoXY3->SetYTitle("radial [m] ");

  // get a fresh canvas, tune it as you like
  TCanvas* myCanv = new TCanvas("myCanv","H5PartROOT",1000, 750);  
  myCanv->Divide(1,2);
  myCanv_1->Divide(2,1);
  myCanv_2->Divide(3,1);

  myCanv_1->cd(1);

  TLegend *legend1=new TLegend(0.7,0.75,0.95,0.9);    
  legend1->SetTextSize(0.04);
  legend1->SetFillColor(0);
  // legend1->SetMargin(0.2);

  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.1);
  gPad->SetLogy();
    
  h1dX1->SetXTitle("longitudinal [m]");
  h1dX1->SetFillColor(38);

  //  h1dX2->SetFillColor(42);
  //  h1dX3->SetFillColor(5);

  h1dX1->SetLineColor(1);
  h1dX2->SetLineColor(2);
  h1dX3->SetLineColor(3);
  h1dX1->SetLabelSize(0.05,"X");
  h1dX1->SetLabelSize(0.05,"Y");

  legend1->AddEntry(h1dX1,"2 deg.","f");
  legend1->AddEntry(h1dX2,"6 deg.","l");
  legend1->AddEntry(h1dX3,"10 deg.","l");
  legend1->SetHeader("Initial phase width");
    
  h1dX1->Draw();
  h1dX2->Draw("LSAME");
  h1dX3->Draw("LSAME");
  legend1->Draw();
  //  delete legend1;
  
  myCanv_1->cd(2);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.1);
  gPad->SetLogy();
  
  h1dY1->SetXTitle("radial [m] ");
  h1dY1->SetFillColor(38);

  // h1dY2->SetFillColor(42);
  // h1dY3->SetFillColor(5);
  
  h1dY1->SetLineColor(1);
  h1dY2->SetLineColor(2);
  h1dY3->SetLineColor(3);
  h1dY1->SetLabelSize(0.05,"X");
  h1dY1->SetLabelSize(0.05,"Y");

  TLegend *legend2=new TLegend(0.7,0.75,0.95,0.9);    
  legend2->SetTextSize(0.04);
  legend2->SetFillColor(0);   
  // legend2->SetMargin(0.2);

  legend2->AddEntry(h1dY1,"2 deg.","f");
  legend2->AddEntry(h1dY2,"6 deg.","l");
  legend2->AddEntry(h1dY3,"10 deg.","l");
  legend2->SetHeader("Initial phase width");
    
  h1dY1->Draw();
  h1dY2->Draw("LSAME");
  h1dY3->Draw("LSAME");
  legend2->Draw();
  // delete legend2;

  myCanv_2->cd(1);
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.1);

  TLegend *legend3=new TLegend(0.7,0.75,0.85,0.8);
  legend3->SetTextSize(0.04);
  legend3->SetFillColor(0);   
  legend3->SetHeader("2 deg.");
  
  histoXY1->Draw("colz");
  legend3->Draw();
  //  delete legend3;
  
  myCanv_2->cd(2);
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.1);

  TLegend *legend4=new TLegend(0.7,0.75,0.85,0.8);
  legend4->SetTextSize(0.04);
  legend4->SetFillColor(0);   
  legend4->SetHeader("6 deg.");
  
  histoXY2->Draw("colz");
  legend4->Draw();
  //  delete legend4;
  
  myCanv_2->cd(3);
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.1);

  
  TLegend *legend5=new TLegend(0.7,0.75,0.85,0.8);
  legend5->SetTextSize(0.04);
  legend5->SetFillColor(0);   
  legend5->SetHeader("10 deg.");


  histoXY3->Draw("colz");
  legend5->Draw();
  //delete legend5;
  
  myCanv->Update();

  TString s1 = TString("Turn-")+Nstep+".pdf";
  TString s2 = TString("Turn-")+Nstep+".jpg";
  
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
