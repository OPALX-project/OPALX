void Fwhm(TH1D* h1D, Double_t &peak, Double_t &fwhm)
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

void printLinegr(TNtuple *nt, TString ntStr, TCanvas *gr, Double_t maxvalX, Double_t maxvalY,
		 TString xaxisTitle, TString yaxisTitle, TString baseFn) {

  gr->cd();
  TVirtualPad* myPad = gr->GetPad(0);
  TH1F*  frame1 = myPad->DrawFrame(0.0,0.0,1.1*maxvalX,1.1*maxvalY);
  frame1->SetXTitle(xaxisTitle);
  frame1->SetYTitle(yaxisTitle);
  frame1->GetYaxis()->SetTitleOffset(1.75);
    
  nt->Draw(ntStr,"","goff");
  TGraph *gr1 = new TGraph(nt->GetSelectedRows(),nt->GetV2(),nt->GetV1());
  gr1->SetLineColor(2);
  gr1->SetLineStyle(1);
  gr1->SetLineWidth(4);
  gr1->Draw("Lsame");
  
  gr->Update();

  TString ofn1("pngs/" + baseFn + "-" + ntStr + ".png");
  gr->Print(ofn1);
    
  TString ofn2("pdfs/" + baseFn + "-" + ntStr + ".pdf");
  gr->Print(ofn2);
}


void divPrint() 
{
gSystem->Load("~schietinger/bin/H5root.so");
TH5Style::SetStyle();
gStyle->SetOptStat(0);

Int_t palette[50] = { 0, 0, 0, 0, 55, 56, 57, 58, 59, 60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100});
gStyle->SetPalette(50,palette); 

const Int_t dataSize = 4;
const Int_t Scenarios = 10;

TString cases[dataSize];
TString scenarios[Scenarios];

TH5Dataset *data = 0;

// **********************************************
bool doHistoPlots = false;
bool doLinePlots = true;

Int_t maxZ = 50;

const Int_t loopScenarioStart = 8;
const Int_t loopScenarioEnd = 10;

const Int_t loopCaseStart = 0;
const Int_t loopCaseEnd = 1;

cases[0] = TString("1");  
cases[1] = TString("2");  
cases[2] = TString("3");
cases[3] = TString("5");

scenarios[0] = TString("10");
scenarios[1] = TString("11");
scenarios[2] = TString("12");
scenarios[3] = TString("13");

scenarios[4] = TString("6");
scenarios[5] = TString("7");
scenarios[6] = TString("8");
scenarios[7] = TString("9");
scenarios[8] = TString("14");
scenarios[9] = TString("15");

// **********************************************

Int_t subsampling = 1;

// for density plots
TCanvas* graphs = new TCanvas("c1Graphs","Graphs",1000,800);
graphs->SetBottomMargin(4.);
graphs->Draw();
graphs->Divide(2,2);

// for line plots 
TCanvas* lgraphs = new TCanvas("lgraphs","lgr",1000,800);
lgraphs->SetBottomMargin(4.);
lgraphs->Draw();

TPaveLabel globTitle;
globTitle.SetBorderSize(0);
globTitle.SetFillColor(0);
globTitle.SetTextFont(62);   

graphs->cd(1);
gPad->SetRightMargin(0.25);
gPad->SetLeftMargin(0.2);
gPad->SetBottomMargin(0.2);

graphs->cd(2);
gPad->SetRightMargin(0.25);
gPad->SetLeftMargin(0.2);
gPad->SetBottomMargin(0.2);

graphs->cd(3);
gPad->SetRightMargin(0.25);
gPad->SetLeftMargin(0.2);
gPad->SetBottomMargin(0.2);

graphs->cd(4);
gPad->SetRightMargin(0.25);
gPad->SetLeftMargin(0.2);
gPad->SetBottomMargin(0.2);

gStyle->SetTickLength(0.03,"x");
gStyle->SetTickLength(0.01,"y");
gStyle->SetPadTickX(0);
gStyle->SetPadTickY(0);

const Int_t DATATOREAD=13; 
Double_t maxVal[DATATOREAD];	
Double_t minVal[DATATOREAD];	
Double_t x[DATATOREAD];

Double_t fwhmx,fwhmz,fwhmy,peak;
TNtuple *nt = 0;


TH1D *h1 = 0;
TH1D *h2 = 0;
TH1D *h3 = 0;

for (Int_t scenario=loopScenarioStart; scenario < loopScenarioEnd; scenario++) {
  TString globtit(TString("Ring-FlatTopScenario-" + scenarios[scenario]));
  for (Int_t c=loopCaseStart; c<loopCaseEnd; c++) {
    TString fn(globtit + "/case" + cases[c] + "/Ring-" + cases[c] + ".h5");

    for (Int_t i=0; i<=DATATOREAD; i++) {
      maxVal[i] = -10000.0;
      minVal[i] =  10000.0;
    }

    if (data)
      delete data;
    data = new TH5Dataset(fn,subsampling);
    if (nt) delete nt;
    nt = new TNtuple("ntuple","data from ascii =file","s:rmsx:minx:maxx:fwhmx:rmsz:minz:maxz:fwhmz:rmsy:miny:maxy:fwhmy");

    cout << "Loaded " << fn << " max number of steps " << (*data).GetNSteps() << endl;

    if (doLinePlots) {
      //      for (Int_t s = 0; s<(*data).GetNSteps(); s++) {
      for (Int_t s = 0; s<3; s++) {
	x[0] = data->GetSpos(s);
	cout << "Work on step " << s << endl;
	// x coordnate	  
	if(h1) delete h1;
	h1 = data->Histo1d("x",s,-0.04,0.04,256);         
	Fwhm(h1,peak,fwhmx);
	x[1] = h1->GetRMS();
	x[2] = h1->GetMinimum();
	x[3] = h1->GetMaximum();
	x[4] = fwhmx;
	
	// z coordnate	 
	if(h2) delete h2; 
	h2 = data->Histo1d("z",s,-0.13,0.13,256);         
	Fwhm(h2,peak,fwhmz);
	x[5] = h2->GetRMS();
	x[6] = h2->GetMinimum();
	x[7] = h2->GetMaximum();
	x[8] = fwhmz;
	
	// y coordnate	  
	if(h3) delete h3;
	h3 = data->Histo1d("y",s,-0.13,0.13,256);         
	Fwhm(h3,peak,fwhmy);
	x[9] =  h3->GetRMS();
	x[10] = h3->GetMinimum();
	x[11] = h3->GetMaximum();
	x[12] = fwhmy;
    
	nt->Fill(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12]);
      
	if (maxVal[0] < x[0])
	  maxVal[0]= x[0];

	for (Int_t i=1; i<=DATATOREAD;i++) {
	  if (maxVal[i] < x[i])
	    maxVal[i]= x[i];
	  else if (minVal[i] > x[i])
	    minVal[i]= x[i];
	}
      }
      
      TString ntStr1;
      TString baseFn(globtit + "-case" + cases[c]);
      
      ntStr1 = TString("rmsx:s");
      printLinegr(nt, ntStr1, lgraphs, maxVal[0], maxVal[4], TString("s [m]"), TString(" x rms [m]"), baseFn);
      
      ntStr1 = TString("rmsy:s");
      printLinegr(nt, ntStr1, lgraphs, maxVal[0], maxVal[12], TString("s [m]"), TString(" y rms [m]"), baseFn);
      
      ntStr1 = TString("rmsz:s");
      printLinegr(nt, ntStr1, lgraphs, maxVal[0], 0.2, TString("s [m]"), TString(" z rms [m]"), baseFn);
      
      
      ntStr1 = TString("fwhmx:s");
      printLinegr(nt, ntStr1, lgraphs, maxVal[0], maxVal[4], TString("s [m]"), TString(" x fwhm [m]"), baseFn);
      
      ntStr1 = TString("fwhmy:s");
      printLinegr(nt, ntStr1, lgraphs, maxVal[0], maxVal[12], TString("s [m]"), TString(" y fwhm [m]"), baseFn);
      
      ntStr1 = TString("fwhmz:s");
      printLinegr(nt, ntStr1, lgraphs, maxVal[0], 0.2, TString("s [m]"), TString(" z fwhm [m]"), baseFn);

      if (h1)
	delete h1;
      if (h2)
	delete h2;
      if (h3)
	delete h3;
    }
    
    if (doHistoPlots) {
      graphs->cd(1);
      Int_t s=0;
      TH2D* histo1 = (*data).Histo2d("x","z",s,-0.015,0.015,-0.1,0.1,512,512);
      histo1->SetMaximum(maxZ);
      histo1->Draw("colz");
      
      graphs->cd(2);
      s=60;
      TH2D* histo2 = (*data).Histo2d("x","z",s,-0.015,0.015,-0.1,0.1,512,512);
      histo2->SetMaximum(maxZ);
      histo2->Draw("colz");
      
      graphs->cd(3);
      s=120; 
      TH2D* histo3 = (*data).Histo2d("x","z",s,-0.015,0.015,-0.1,0.1,512,512);
      histo3->SetMaximum(maxZ);
      histo3->Draw("colz");
      
      graphs->cd(4);
      s= (*data).GetNSteps()-1;
      TH2D* histo4 = (*data).Histo2d("x","z",s,-0.015,0.015,-0.1,0.1,512,512);
      histo4->SetMaximum(maxZ);
      histo4->Draw("colz");

      graphs->cd(1);
      globTitle.DrawPaveLabel(0.1,0.95,0.9,0.99,fn,"NDC");
      
      graphs->Update();
      
      TString ofn;
      ofn = TSting("pngs/" + globtit + "-case" + cases[c] + "-1.png");
      graphs->Print(ofn);
      
      ofn = TString("pdfs/" + globtit + "-case" + cases[c] + "-1.pdf");
      graphs->Print(ofn);
    }
  }
 }
}

