void Fwhm(TH1F* h1D, Double_t &peak, Double_t &fwhm)
{
  Float_t yp = h1D->GetMaximum();
 
  Int_t lowBin = 0;
  Int_t hiBin  = h1D->GetNbinsX()+1;

  do {
    lowBin++;
  } while (h1D->GetBinContent(lowBin)<= yp/2);

  do {
    hiBin--;
  } while (h1D->GetBinContent(hiBin)<= yp/2);

  Double_t minX = h1D->GetXaxis()->GetBinLowEdge(lowBin);
  Double_t maxX = h1D->GetXaxis()->GetBinUpEdge(hiBin);

  if (minX < 0.0)
    minX *= -1.0;

  if (maxX < 0.0)
    maxX *= -1.0;

  fwhm = minX + maxX;
  peak = yp;
}

void compareDist()
{
   char *myVar("t");
   Int_t nbins = 300;	

   gStyle->SetOptStat("emri");
   
   TCanvas *c1 = new TCanvas("c1","OPAL Astra comparison",200,10,700,500);
   TPad *pad1 = new TPad("pad1","",0,0,1,1);
 
   pad1->SetFillStyle(4000); //will be transparent

   pad1->Draw();
   pad1->cd();
	
   Double_t maxo = nto->GetMaximum(myVar);
   Double_t maxa = nta->GetMaximum(myVar);

   Double_t xmax = 0.0;
   Double_t xmin = 0.0;

   if (maxo>maxa)
     xmax=maxo*1.1;
   else
     xmax=maxa*1.1;	

   if (myVar != "t")
	   xmin = -xmax;
   else {
	   Double_t mino = nto->GetMinimum(myVar);
	   Double_t mina = nta->GetMinimum(myVar);
	   if (mino<mina)
		   xmin = mino;
	   else
		   xmin = mina;
   }

   TH1F* h1 = new TH1F("OPAL","",nbins,xmin,xmax);
   TH1F* h2 = new TH1F("Astra","",nbins,xmin,xmax);

   nto->Project("OPAL" ,myVar);
   nta->Project("Astra",myVar);

   Double_t peako = 0;
   Double_t fwhmo = 0;
   Double_t peaka = 0;
   Double_t fwhma = 0;

   //   Double_t scaleh1 = 1/h1->Integral();
   //h1->Scale(scaleh1);
   //h1->GetXaxis ()->SetTitle(Form("%s (m)",myVar));
   //Double_t scaleh2 = 1/h2->Integral();
   //h2->Scale(scaleh2);

   Fwhm(h1,peako,fwhmo);
   Fwhm(h2,peaka,fwhma);

   cout << "peakOPAL = " << peako << " fwhmOPAL = " << fwhmo << endl;
   cout << "peakAstra= " << peaka << " fwhmAstra= " << fwhma << endl;


   h1->SetLineColor(kRed);
   h1->Draw();
   pad1->Update(); //this will force the generation of the "stats" box
   TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
   ps1->SetTextColor(kRed);
   pad1->Modified();
   

   h2->Draw("sames");
   pad1->Update();
  
   TPaveStats *ps2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
   ps2->SetX1NDC(0.15); 
   ps2->SetX2NDC(0.35);

   
   leg_hist = new TLegend(0.4,0.3,0.67,0.47);
   //  leg_hist->SetHeader(""); 

   if (myVar != "t") {
	   leg_hist->AddEntry(h1,Form("OPAL  Fwhm=%g (m)",fwhmo),"l");
	   leg_hist->AddEntry(h2,Form("Astra Fwhm=%g (m)",fwhma),"l");
   }
   else {
	   leg_hist->AddEntry(h1,Form("OPAL  Fwhm=%g (ps)",(fwhmo*1E12)),"l");
	   leg_hist->AddEntry(h2,Form("Astra Fwhm=%g (ps)",(fwhma*1E12)),"l");
   }
   leg_hist->Draw();

   c1->SaveAs(Form("../figures/opal-astra-%s.pdf",myVar));


  /* --------------------------------
    Calculate the covariance matrix     
   --------------------------------
  */
  Float_t xtmp;
  Float_t pxtmp;
  Float_t ytmp;
  Float_t pytmp;
  Float_t ttmp;
  Float_t pttmp;
  Double_t data[6];

  nto->SetBranchAddress("x",&xtmp);
  nto->SetBranchAddress("px",&pxtmp);
  nto->SetBranchAddress("y",&ytmp);
  nto->SetBranchAddress("py",&pytmp);
  nto->SetBranchAddress("t",&ttmp);
  nto->SetBranchAddress("pt",&pttmp);

  TPrincipal po( 6, "ND");

  for ( int ient = 0; ient < nto->GetEntries(); ient++ ) {
    nto->GetEntry(ient);
    data[0] = xtmp;
    data[1] = pxtmp;
    data[2] = ytmp;
    data[3] = pytmp;    
    data[4] = ttmp;
    data[5] = pttmp;
    po.AddRow( data );
  }

  TMatrixD* mo = po.GetCovarianceMatrix();
  const TVectorD diago = TMatrixDDiag_const(*mo); 
  mo->NormByDiag(diago);
  cout << "OPAL "; mo->Print("f= %6.5g");  


  nta->SetBranchAddress("x",&xtmp);
  nta->SetBranchAddress("px",&pxtmp);
  nta->SetBranchAddress("y",&ytmp);
  nta->SetBranchAddress("py",&pytmp);
  nta->SetBranchAddress("t",&ttmp);
  nta->SetBranchAddress("pt",&pttmp);

  TPrincipal pa( 6, "ND");

  for ( int ient = 0; ient < nta->GetEntries(); ient++ ) {
    nta->GetEntry(ient);
    data[0] = xtmp;
    data[1] = pxtmp;
    data[2] = ytmp;
    data[3] = pytmp;    
    data[4] = ttmp;
    data[5] = pttmp;
    pa.AddRow( data );
  }

  TMatrixD* ma = pa.GetCovarianceMatrix();
  const TVectorD diaga = TMatrixDDiag_const(*ma); 
  ma->NormByDiag(diaga);
  cout << "Astra "; ma->Print("f= %6.5g"); 














}
   
      
