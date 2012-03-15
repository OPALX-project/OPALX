// Post processing script for OPAL-cycl in the Global frame
// Run it in an interactive root session like this (make sure you have files
// called Ring-1.h5 Ring-2.h5 and Ring-3.h5, N is the dumping step ):
//    
//   .x AnalyzeG.C(N)
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
void AnalyzeG(int Step=0){
  // load the shared library:
  gSystem->Load("/home5/schietinger/extlib/H5Root/lib/libh5root.so");
  
  
  // let's do it in style:
  TH5Style::SetStyle();

  gStyle->SetOptStat(0);
  
  // load the data file:
  
  TH5Dataset *data1 = new TH5Dataset("Ring-1/Ring-1.h5");
  TH5Dataset *data2 = new TH5Dataset("Ring-2/Ring-2.h5");
  TH5Dataset *data3 = new TH5Dataset("Ring-3/Ring-3.h5");

  // for convenience, get a dummy canvas for plotting
  TCanvas* dummy = new TCanvas("dummy","dummy",10, 10);      

  // get a 2D x-y histogram for step 100:
  //  TH2D* histoZPZ1 = data1.Histo2d("z","pz",100,-0.008,0.008,-0.002, 0.002);
  //S1 100 turn
  
  double Xmin=-1.55;
  double Xmax=-1.40;
  double Ymin= 3.49;
  double Ymax= 3.59;
  /*
  S5 100 turn
  double Xmin=-1.60;
  double Xmax=-1.40;
  double Ymin= 3.60;
  double Ymax= 3.69;
  */
  /*
  //0 turn
  
  double Xmin=-0.82;
  double Xmax=-0.72;
  double Ymin= 1.86;
  double Ymax= 1.92;
  */
  TH2D* histoXY1 = data1.Histo2d("x","y",Step,Xmin,Xmax,Ymin,Ymax,200,200);
  TH2D* histoXY2 = data2.Histo2d("x","y",Step,Xmin,Xmax,Ymin,Ymax,200,200);
  TH2D* histoXY3 = data3.Histo2d("x","y",Step,Xmin,Xmax,Ymin,Ymax,200,200);
  TH1D* h1dX1 = data1->Histo1d("x",Step,Xmin,Xmax,256);
  TH1D* h1dX2 = data2->Histo1d("x",Step,Xmin,Xmax,256);
  TH1D* h1dX3 = data3->Histo1d("x",Step,Xmin,Xmax,256);
  TH1D* h1dY1 = data1->Histo1d("y",Step,Ymin,Ymax,256);
  TH1D* h1dY2 = data2->Histo1d("y",Step,Ymin,Ymax,256);
  TH1D* h1dY3 = data3->Histo1d("y",Step,Ymin,Ymax,256);

  delete dummy;

  // embellish the histograms at wish and draw it with the desired 
  // drawing option:
  histoXY1->SetLabelSize(0.05,"X");
  histoXY1->SetLabelSize(0.05,"Y");
  histoXY1->SetXTitle("x [m]"); 	
  histoXY1->SetYTitle("y [m]"); 	

  
  histoXY2->SetLabelSize(0.05,"X");
  histoXY2->SetLabelSize(0.05,"Y");
  histoXY2->SetXTitle("x [m]"); 	
  histoXY2->SetYTitle("y [m]"); 	

  histoXY3->SetLabelSize(0.05,"X");
  histoXY3->SetLabelSize(0.05,"Y");
  histoXY3->SetXTitle("x [m]"); 	
  histoXY3->SetYTitle("y [m]"); 	
  

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
    
  h1dX1->SetXTitle("x [m] ");
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
  
  h1dY1->SetXTitle("y [m] ");
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

  
  TString s1 = TString("Turn-")+Step+"-Global.pdf";
  TString s2 = TString("Turn-")+Step+"-Global.jpg";
  
  myCanv->Print(s1);
  myCanv->Print(s2);

}
