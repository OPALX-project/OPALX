// NB total bunch number in the multi-bunch h5 file
// Nstep: the step of single bunch h5 file, which corresponding to  Nstep+(NB-1)/2
void C9B_7B_SB_R_1uA(int Nstep=0, int NB1=1, int NB2=1){

  // load the shared library:
  //gSystem->Load("/home5/schietinger/extlib/H5Root/lib/libh5root.so");
  gSystem->Load("/home5/schietinger/bin/libh5root-yang.so");

  
  // let's do it in style:
  TH5Style::SetStyle();
  // get rid of statistic box
  gStyle->SetOptStat(0);
  
  // load the data file:
  TH5Dataset *data1 = new TH5Dataset("Scenario-3/Ring-1.h5"); //1 bunch   
  TH5Dataset *data2 = new TH5Dataset("Scenario-7/Ring-1.h5"); //7 bunches 
  TH5Dataset *data3 = new TH5Dataset("Scenario-6/Ring-1.h5"); //9 bunches 
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

  double XminL=-10;
  double XmaxL=10;

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
  /////////////////
  TNtupleD *ntd = data1->GetNtuple(Nstep_SB);
  TH1D* h1dX1 = new TH1D("h1dX1","h1 title",256,XminL,XmaxL);
  ntd->Draw(Form("1000.0*((x-%g)*%g+(y-%g)*%g)>>h1dX1",cx[Nstep_SB],cosAngle[0],cy[Nstep_SB],sinAngle[0]));
  delete ntd;
  /////////////////
  TNtupleD *ntd = data2->GetNtuple(Nstep_MB1);
  TH1D* h1dX2 = new TH1D("h1dX2","h1 title",256,XminL,XmaxL);
  ntd->Draw(Form("1000.0*((x-%g)*%g+(y-%g)*%g)>>h1dX2",cx[Nstep_SB],cosAngle[0],cy[Nstep_SB],sinAngle[0]));
  delete ntd;
  /////////////////
  TNtupleD *ntd = data3->GetNtuple(Nstep_MB2);
  TH1D* h1dX3 = new TH1D("h1dX3","h1 title",256,XminL,XmaxL);
  ntd->Draw(Form("1000.0*((x-%g)*%g+(y-%g)*%g)>>h1dX3",cx[Nstep_SB],cosAngle[0],cy[Nstep_SB],sinAngle[0]));
  delete ntd;
  /////////////////  
  delete dummy;




  

  //--------------plot--------------------//
  TCanvas* myCanv = new TCanvas("myCanv","H5PartROOT",1000, 750);  
  myCanv->cd(1);

  TLegend *legend1=new TLegend(0.67,0.40,0.95,0.55);    
  legend1->SetTextSize(0.05);
  legend1->SetFillColor(0);
  // legend1->SetMargin(0.2);

  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.2);

  //  gPad->SetLogy();

  h1dX1->SetLabelSize(0.07,"X");
  h1dX1->SetLabelSize(0.07,"Y");
  h1dX1->SetTitleSize(0.07,"X");
  h1dX1->SetTitleSize(0.07,"Y");
  h1dX1->SetXTitle("transverse (mm)");
  h1dX1->SetYTitle("particles number");
  h1dX1->SetFillColor(38);
  h1dX1->SetAxisRange(0, 5000, "Y");
  
  h1dX1->SetLineColor(1);
  h1dX2->SetLineColor(2);
  h1dX2->SetLineWidth(5);  

  h1dX3->SetLineColor(3);
  h1dX3->SetLineWidth(5);

  
  legend1->AddEntry(h1dX1,"1   bunch","f");
  legend1->AddEntry(h1dX2,"7 bunches","l");
  legend1->AddEntry(h1dX3,"9 bunches","l");

  h1dX1->Draw();
  h1dX2->Draw("LSAME");
  h1dX3->Draw("LSAME");
  legend1->Draw();
  myCanv->Update();

  TString s1 = TString("C9B7BSB-R-1uA-")+Nstep+".pdf";
  TString s2 = TString("C9B7BSB-R-1uA-")+Nstep+".png";

  
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
