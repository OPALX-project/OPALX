{
/* 
  load the shared library:
*/
  gSystem->Load("/home5/schietinger/extlib/H5Root/lib/libh5root.so"); 
  using namespace TH5Util;

  TH5Style::SetStyle();
  
  gStyle->SetPadTickY(0);
  gStyle->SetLineColor(2);
  gStyle->SetLineWidth(3);
  gStyle->SetMarkerColor(4);
  gStyle->SetMarkerStyle(21);
  gStyle->SetOptStat(11001111);
  gStyle->SetHistFillColor(kGray);

  Width_t lineWidth(3);
  Int_t nLines = 0;

  TString fnStr("../binomial.dist");

  c1 = new TCanvas("c1","Analyze Binomial Distribution",200,10,1200,1000);

  c1.Divide(3,3);

  TNtuple *ntuple =new TNtuple("ntuple","data from ascii =file","x:y:t:px:py:pt");

  /* 
     read in file data
  */
  ifstream in; 
  in.open(fnStr);
  in.ignore(64,'\n');

  while (1) {    
      Float_t x,y,t,px,py,pt;
      in >> x >> y >> t >> px >> py >> pt;
      if (!in.good()) 
	  break;      
      ntuple->Fill(x,y,t,px,py,pt);
      nLines++;   
  }   
  in.close();
  printf("Found %d lines \n",nLines);

  /* ---------------------------
         
     ---------------------------
  */
  Double_t maxX = ntuple->GetMaximum("x");

  c1.cd(1);
  ntuple->Draw("x");

  c1.cd(2);
  ntuple->Draw("px");

  c1.cd(3);
  ntuple->Draw("px:x","","colz");

  c1.cd(4);
  ntuple->Draw("y");

  c1.cd(5);
  ntuple->Draw("py");

  c1.cd(6);
  ntuple->Draw("py:y","","colz");


  c1.cd(7);
  ntuple->Draw("t");

  c1.cd(8);
  ntuple->Draw("pt");

  c1.cd(9);
  ntuple->Draw("pt:t","","colz");

  c1.Print("binomial.pdf");

  /* --------------------------------
     Calculate the covariance matrix     
     --------------------------------
  */
  Float_t x;
  Float_t px;
  Float_t y;
  Float_t py;
  Float_t t;
  Float_t pt;

  ntuple->SetBranchAddress("x",&x);
  ntuple->SetBranchAddress("px",&px);
  ntuple->SetBranchAddress("y",&y);
  ntuple->SetBranchAddress("py",&py);
  ntuple->SetBranchAddress("t",&t);
  ntuple->SetBranchAddress("pt",&pt);

  TPrincipal p( 6, "ND");
  Double_t data[6];
  for ( int ient = 0; ient < ntuple->GetEntries(); ient++ ) {
    ntuple->GetEntry(ient);
    data[0] = x;
    data[1] = px;
    data[2] = y;
    data[3] = py;    
    data[4] = t;
    data[5] = pt;
    p.AddRow( data );
  }

  TMatrixD* m = p.GetCovarianceMatrix();
  const TVectorD diag = TMatrixDDiag_const(*m); 
  m->NormByDiag(diag);
  m->Print("f= %6.3f");  
}



