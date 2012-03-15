#include "Riostream.h"
#include <fstream>
#include <iomanip>
#include <iostream>

void decorateCanvas(TString globtit)
{
  TPaveLabel* globTitle = new TPaveLabel(0.1,0.95,0.9,0.99,globtit);
  globTitle->SetBorderSize(0);
  globTitle->SetFillColor(0);
  globTitle->SetTextFont(62);
  globTitle->Draw();
    
  TDatime now;
  TPaveLabel* date = new TPaveLabel(0.7,0.01,0.9,0.05,now.AsString());
  date->SetBorderSize(0);
  date->SetFillColor(0);
  date->SetTextFont(62);
  date->SetTextSize(0.3);
  date->Draw();
  
  TPad* graphs = new TPad("c1Graphs","Graphs",0.01,0.05,0.95,0.95);
  graphs->SetBottomMargin(4.);
  graphs->Draw();
  graphs->Divide(1,1);
  
  graphs->cd(1);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  gStyle->SetTickLength(0.03,"x");
  gStyle->SetTickLength(0.01,"y");
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);

}

TPad* decorate4Canvases(TString globtit)
{
  TPaveLabel* globTitle = new TPaveLabel(0.1,0.95,0.9,0.99,globtit);
  globTitle->SetBorderSize(0);
  globTitle->SetFillColor(0);
  globTitle->SetTextFont(62);
  globTitle->Draw();
    
  TDatime now;
  TPaveLabel* date = new TPaveLabel(0.7,0.01,0.9,0.05,now.AsString());
  date->SetBorderSize(0);
  date->SetFillColor(0);
  date->SetTextFont(62);
  date->SetTextSize(0.3);
  date->Draw();
  
  TPad* graphs = new TPad("c1Graphs","Graphs",0.01,0.05,0.95,0.95);
  graphs->SetBottomMargin(4.);
  graphs->Draw();
  graphs->Divide(2,2);
  
  graphs->cd(1);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  graphs->cd(2);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  graphs->cd(3);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  graphs->cd(4);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  gStyle->SetTickLength(0.03,"x");
  gStyle->SetTickLength(0.01,"y");
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);

  return graphs;
}



TPad* decorate6Canvases(TString globtit)
{
  TPaveLabel* globTitle = new TPaveLabel(0.1,0.95,0.9,0.99,globtit);
  globTitle->SetBorderSize(0);
  globTitle->SetFillColor(0);
  globTitle->SetTextFont(62);
  globTitle->Draw();
    
  TDatime now;
  TPaveLabel* date = new TPaveLabel(0.7,0.01,0.9,0.05,now.AsString());
  date->SetBorderSize(0);
  date->SetFillColor(0);
  date->SetTextFont(62);
  date->SetTextSize(0.3);
  date->Draw();
  
  TPad* graphs = new TPad("c1Graphs","Graphs",0.01,0.05,0.95,0.95);
  graphs->SetBottomMargin(4.);
  graphs->Draw();
  graphs->Divide(2,3);
  
  graphs->cd(1);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  graphs->cd(2);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  graphs->cd(3);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  graphs->cd(4);
  gPad->SetRightMargin(0.25);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  
  graphs->cd(5);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  graphs->cd(6);
  gPad->SetRightMargin(0.1);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);

  gStyle->SetTickLength(0.03,"x");
  gStyle->SetTickLength(0.01,"y");
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);

  return graphs;
}




// Change the line style
//
//   linestyle = 2 dashed
//             = 3  dotted
//             = 4  dash-dotted
//              else = solid


void drawFwhm(TString baseDir, Double_t* maxVal, TString globtit, TNtuple* nt[10], Int_t startCase, Int_t endCase) 
{
  TGraph*  gr1[10];
  TGraph*  gr2[10];
  TGraph*  gr3[10];
  TGraph*  gr4[10];
  TString num1("");
  TString num2("");
  
  num1 += startCase;
  num2 += endCase;

  c1 = new TCanvas("c1",baseDir,800,1000);
  decorateCanvas(globtit);
  TH1F*  frame = c1->DrawFrame(0,0,1.1*maxVal[0],1.1*maxVal[1]);
  frame->SetXTitle("s [m] ");
  frame->SetYTitle("FWHM_{x}, RMS_{x} [m] ");
  frame->GetYaxis()->SetTitleOffset(1.75);

  TLegend *legend=new TLegend(0.7,0.65,0.88,0.85);    
  legend->SetTextSize(0.03);
  legend->SetFillColor(0);
  legend->SetMargin(0.2);
  
  Int_t pstype = 111;
  TString psfn = baseDir+TString("-fwhm")+"case"+num1+"-"+num2+TString(".ps");
  TPostScript *ps = new TPostScript(psfn,pstype);
  ps->NewPage();               
  for (Int_t d=startCase; d<endCase; d++) {
    nt[d]->Draw("fwhmx:s","","goff");
    gr1[d] = new TGraph(nt[d]->GetSelectedRows(),nt[d]->GetV2(),nt[d]->GetV1());
    gr1[d]->SetLineColor(d);
    gr1[d]->SetLineStyle(1);
    gr1[d]->Draw("Lsame");
    legend->AddEntry(gr1[d],Form("FWHM case %d",d),"l");

    nt[d]->Draw("rmsx:s","","goff");
    gr2[d] = new TGraph(nt[d]->GetSelectedRows(),nt[d]->GetV2(),nt[d]->GetV1());
    gr2[d]->SetLineColor(d);
    gr2[d]->SetLineStyle(2);
    gr2[d]->Draw("Lsame");
    legend->AddEntry(gr2[d],Form("RMS case %d",d),"l");
    /*
    nt[d]->Draw("maxx:s","","goff");
    gr3[d] = new TGraph(nt[d]->GetSelectedRows(),nt[d]->GetV2(),nt[d]->GetV1());
    gr3[d]->SetLineColor(d+2);
    gr3[d]->SetLineStyle(3);
    gr3[d]->Draw("Lsame");
    legend->AddEntry(gr3[d],Form("maxx case %d",d),"l");

    nt[d]->Draw("minx:s","","goff");
    gr4[d] = new TGraph(nt[d]->GetSelectedRows(),nt[d]->GetV2(),nt[d]->GetV1());
    gr4[d]->SetLineColor(d+2);
    gr4[d]->SetLineStyle(3);
    gr4[d]->Draw("Lsame");
    legend->AddEntry(gr4[d],Form("minx case %d",d),"l");
    */
  }
  legend->Draw();
  c1->Update();
  ps->Close();
}

void Fwhm(TH1D* h1D,Double_t &peak,Double_t &fwhm)
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

void Fwhm1(TH1D* h1D,Double_t &peak,Double_t &fwhm)
{
  Int_t npeaks = 2;

  TSpectrum *s = new TSpectrum(2*npeaks);
  Int_t nfound = s->Search(h1D,2,"goff",0.001);
  Float_t *xloc = s->GetPositionX();
  Int_t bin = h1D->GetXaxis()->FindBin(*xloc);
  Float_t yp = h1D->GetBinContent(bin);

  Int_t lowBin;
  Int_t hiBin;

  lowBin = 0;
  do {
    lowBin++;
  } while (h1D->GetBinContent(lowBin)<= yp/2);
  
  hiBin = h1D->GetNbinsX()+1;
  do {
    hiBin--;
  } while (h1D->GetBinContent(hiBin)<= yp/2);


  //  cout << "Foud " << nfound <<  " candidate peak at " << *xloc << " value " << yp << " bin " 
  //     << bin << " nbins " << h1D->GetNbinsX() << " lowBin " << lowBin << " hiBin " << hiBin 
  //     << " val low " <<h1D->GetBinContent(lowBin) << " val hi " << h1D->GetBinContent(hiBin) << endl;

  Double_t minX = h1D->GetXaxis()->GetBinLowEdge(lowBin);
  Double_t maxX = h1D->GetXaxis()->GetBinUpEdge(hiBin);

  if (minX < 0.0)
    minX *= -1.0;

  if (maxX < 0.0)
    maxX *= -1.0;

  fwhm = minX + maxX;
  peak = yp;
}


void writeLeftWWWWindow(ofstream  &ostr, TString scenarioStr, TString casenumStr, Int_t scenario) {

  TString wwwId("Scenario"+scenarioStr+"Case"+casenumStr);

  ostr << "<table width='450' cellpadding='0' cellspacing='0' border='0' style='background-color: rgb(255,255,255); padding-top: 30px;'>" << endl;
  ostr << "<TR><TD style='background-color: rgb(255,255,255);'>" << endl;
  ostr << "<h2>Ring Flat-Top <a href=\"Ring-FlatTopScenario-" << scenario << "-summary.png\"/a> Scenario " << scenario  << " case " << casenumStr  << "</h2>" << endl;
  ostr << "<TD id='"+wwwId+"_controls' style='vertical-align: top; background-color: rgb(255,255,255);'>" << endl;
  ostr << "<div><a href=\"javascript: hideDIV('"+wwwId+"');\" style='font-size:20px; text-decoration: none; float: left;'>-</a></div></TD></TR></table>" << endl;

  ostr << "<div id='"+wwwId+"' style='height: 200px; width:450px; overflow: auto;'>" << endl;
  ostr << "<table border='3' frame='void'>" 
       << "<tr>" 
       << "<td><b>Turn</b></TD>"
       << "<td><b>spos [m] </b></td> "
       << "<td><b>V Flat-Top [%] </b></td> "
       << "<td><b>V Accel [MV] </b></td> "
       << "<td><b>Results </b></td> </tr>" << endl;
}


int main() 
{ 
  //  gSystem->Load("~schietinger/bin/H5root.so");           
  Int_t mStyleArr[10] = {2,24,4,5,6,24,25,26,27,28};
  Int_t mColorArr[10] = {1,2,3,4,5,6,7,8,9,10};
 
  Double_t Q,I0,N;
  Int_t nP;
 
  
  TNtuple* nt[10];
  TGraph*  gr[10];

  TH5Dataset* data;
 
  Double_t maxVal[12];	
  Double_t minVal[12];	

  TH5Style::SetStyle();            
  
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  
  Int_t pstype = 111;

  /*
    Hold extreme values to scale output
  */
  Int_t DATATOREAD=12; 

  /*
    Define what to do:
  */

  Int_t startScenario   = 8;
  Int_t endScenario   = 12;

  // cases i.e Ring-1 Ring-2 etc.
  Int_t startCase = 1;
  Int_t endCase  = 5;
  
  // step range if endStep < 0 process all steps in the file 
  Int_t startStep = 1;
  Int_t endStep   = 3;  
  
  Int_t subsampling = 1;
  
  bool doStat     = false;
  bool doStatTest = true;    // only read 500 points 
  bool doPart     = false;
  bool doFWHM     = false;
  bool doCombined = true;
  
  if(doCombined) {
    
    gSystem->mkdir("www");
    
    /*
      Prepare some common data 
    */

    TDatime now;
    TPaveLabel* date = new TPaveLabel(0.7,0.01,0.9,0.05,now.AsString());
    date->SetBorderSize(0);
    date->SetFillColor(0);
    date->SetTextFont(62);
    date->SetTextSize(0.3);

    /*
      Open ScenarioDesc.txt
    */
    ifstream descInF;
    TString dummyS;
    descInF.open("ScenarioDesc.txt");
    
    for (Int_t scenario=startScenario; scenario<=endScenario; scenario++) {
      TString baseDir;
      TString scenarioStr("");	
      scenarioStr+=scenario;	   
      baseDir = TString("Ring-FlatTopScenario-"+scenarioStr);
      TString globtit = baseDir;
      TString myname = "Dude";

      Double_t vFlatTop = 0.0;
      Double_t vAccel = 0.0;

      for (Int_t i=0; i<DATATOREAD; i++) {
	maxVal[i] = -10000.0;
	minVal[i] =  10000.0;
      }

      descInF >> dummyS >> dummyS >> dummyS >> vFlatTop >> dummyS >> vAccel >> dummyS;
 
      c1 = new TCanvas("c1",globtit,800,1000);
      if (c1)
	delete c1;
      c1 = new TCanvas("c1",globtit,1,1);

      cout << " ********************** make www style analysis ****************" << endl;
      cout << "   ********* Analysis of " << globtit << " ********" << endl;
      cout << "   ********* Data from scenario" << startScenario << " to " << endScenario <<   " ********" << endl;
      cout << "   ********* Case " << startCase << " to " << endCase <<   " ********" << endl;
      cout << "   ********* VflatTop " << vFlatTop << " Vaccel " << vAccel <<   " ********" << endl;
      
      cout << " *********************************************************" << endl;
      
      for (Int_t d=startCase; d<=endCase; d++) {
	
	TString fn  = TString(baseDir+"/case"+d+"/Ring-"+d+".h5");
	cout << "Work on " << fn << endl;	
	
	data = new TH5Dataset(fn,subsampling);    
	
	if (endStep<0)
	  endStep = data->GetNSteps();       
	
	cout << "Read in " << endStep << " steps process data from " << startStep << " to " << endStep << endl;
	
	TString casenumStr("");
	casenumStr += d;
	
	if (c1)
	  delete c1;
	c1 = new TCanvas("c1",baseDir,800,1000);
	
	TPad* graphs = NULL;
	TPostScript *ps = NULL;
	
	TString wwwfn("www/index-scenario"+scenarioStr+"-case"+casenumStr+".html");
	ofstream ostr(wwwfn);
	writeLeftWWWWindow(ostr,scenarioStr,casenumStr,scenario);
	
	/*
	  Calculte rms and fwhm grapths
	  
	*/
	Double_t x[9];
	Double_t peak,fwhmx,fwhmz;
	nt[d] = new TNtuple("ntuple","data from ascii =file","s:rmsx:minx:maxx:fwhmx:rmsz:minz:maxz:fwhmz");
	for (Int_t s = startStep; s<endStep; s++) {
	  
	  x[0] = data->GetSpos(s);
	  
	  TH1D *h1 = data->Histo1d("x",s,-0.04,0.04,256);         
	  Fwhm(h1,peak,fwhmx);
	  x[1] = h1->GetRMS();
	  x[2] = h1->GetMinimum();
	  x[3] = h1->GetMaximum();
	  x[4] = fwhmx;
	  
	  // z coordnate	  
	  TH1D* h2 = data->Histo1d("z",s,-0.1,0.1,256);         
	  Fwhm(h2,peak,fwhmz);
	  x[5] = h2->GetRMS();
	  x[6] = h2->GetMinimum();
	  x[7] = h2->GetMaximum();
	  x[8] = fwhmz;

	  nt[d]->Fill(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]);
	
	  if (maxVal[0] < x[0])
	    maxVal[0]= x[0];

	  for (Int_t i=1; i<=8;i++) {
	    if (maxVal[i] < x[i])
	      maxVal[i]= x[i];
	    else if (minVal[i] < x[i])
	      minVal[i]= x[i];
	  }
	  cout << "calculate rms/fwhm at step " << s << " out of " << endStep << " FWHMz " << fwhmz << " FWHMx " << fwhmx << endl;	
	  
	  if (h1)
	    delete h1;
	  if (h2)
	    delete h2;
	}   
	
	for (Int_t s = startStep; s<endStep; s++) {
	  TString stepnum(""); stepnum += s;
	  TString psfn = TString("www/")+baseDir+TString("-combined-case-")+casenumStr+TString("-step")+stepnum+TString(".ps");
	  TString pngfn = baseDir+TString("-combined-case-")+casenumStr+TString("-step")+stepnum+TString(".ps.png");
	  TString globtit = baseDir+TString(" case ")+casenumStr+TString(" turn")+stepnum;
	  Double_t spos = data->GetSpos(s);
	  TString spos_str("");
	  spos_str += spos;
	  
	  if (graphs)
	    delete graphs;
	  if(ps)
	    delete ps;
	  
	  graphs = decorate6Canvases(globtit);
	  ps = new TPostScript(psfn,pstype);
	  
	  cout << "Work on step: " << globtit << endl;	
	  
	  TLatex l; 
	  l.SetNDC();
	  l.SetTextSize(0.1);
	  l.SetTextColor(2);
	  // AAA 	
	  graphs->cd(1);
	  TH1D* h1D = data->Histo1d("x",s,-0.04,0.04,256);         
	  h1D->SetFillColor(45);
	  graphs->cd(2);
	  TH1D* h2D = data->Histo1d("y",s,-0.04,0.04,256);         
	  h2D->SetFillColor(38);
	  graphs->cd(3);
	  TH1D* h3D = data->Histo1d("z",s,-0.1,0.1,256);         
	  h3D->SetFillColor(5);
	  graphs->cd(4);
	  gPad->SetLogz();
	  TH2F* h4D = data->Histo2d("x","z",s,-0.04,0.04,-0.1,0.1,256,256);         
	  h4D->SetStats(kFALSE);
	  h4D->Draw("COLZ");

	  // graphs->cd(5);
	  // l.DrawLatex(0.5,0.5,"Latex test");

	  /*
	    Add line graphs
	
	  */
	  graphs->cd(5);
	  TH1F*  frame1 = c1->DrawFrame(0,0,1.1*maxVal[0],1.1*maxVal[4]);
	  frame1->SetXTitle("s [m] ");
	  frame1->SetYTitle("x [m] ");
	  frame1->GetYaxis()->SetTitleOffset(1.75);

	  TLegend *legend1=new TLegend(0.7,0.65,0.88,0.85);    
	  legend1->SetTextSize(0.03);
	  legend1->SetFillColor(0);
	  legend1->SetMargin(0.2);

	  nt[d]->Draw("fwhmx:s","","goff");
	  TGraph *gr1 = new TGraph(nt[d]->GetSelectedRows(),nt[d]->GetV2(),nt[d]->GetV1());
	  gr1->SetLineColor(d);
	  gr1->SetLineStyle(1);
	  gr1->Draw("Lsame");
	  legend1->AddEntry(gr1,"FWHM","l");
	  legend1->Draw();

	  nt[d]->Draw("rmsx:s","","goff");
	  TGraph *gr2 = new TGraph(nt[d]->GetSelectedRows(),nt[d]->GetV2(),nt[d]->GetV1());
	  gr2->SetLineColor(d+1);
	  gr2->SetLineStyle(1);
	  gr2->Draw("Lsame");
	  legend1->AddEntry(gr2,"RMS","l");
	  legend1->Draw();

	  Double_t myS = data->GetSpos(s);
	  Double_t myY = 1.1*maxVal[4];
	  TLine *line1 = new TLine(myS,0.0,myS,myY);
	  line1->Draw();

	  graphs->cd(6);
	  TH1F*  frame2 = c1->DrawFrame(0,0,1.1*maxVal[0],1.1*maxVal[8]);
	  frame2->SetXTitle("s [m] ");
	  frame2->SetYTitle("z [m] ");
	  frame2->GetYaxis()->SetTitleOffset(1.75);

	  TLegend *legend2=new TLegend(0.7,0.65,0.88,0.85);    
	  legend2->SetTextSize(0.03);
	  legend2->SetFillColor(0);
	  legend2->SetMargin(0.2);

	  nt[d]->Draw("fwhmz:s","","goff");
	  TGraph *gr3 = new TGraph(nt[d]->GetSelectedRows(),nt[d]->GetV2(),nt[d]->GetV1());
	  gr3->SetLineColor(d);
	  gr3->SetLineStyle(1);
	  gr3->Draw("Lsame");
	  legend2->AddEntry(gr3,"FWHM","l");

	  nt[d]->Draw("rmsz:s","","goff");
	  TGraph *gr4 = new TGraph(nt[d]->GetSelectedRows(),nt[d]->GetV2(),nt[d]->GetV1());
	  gr4->SetLineColor(d+1);
	  gr4->SetLineStyle(1);
	  gr4->Draw("Lsame");
	  legend2->AddEntry(gr4,"RMS","l");
	  legend2->Draw();

	  TLine *line2 = new TLine(myS,0.0,myS,1.1*maxVal[8]);
	  line2->Draw();

	  ps->Close();
	  delete h1D;
	  delete h2D;
	  delete h3D;
	  delete h4D;
	  delete gr1;
	  delete gr2;
	  delete gr3;
	  delete gr4;

	  ostr << "<tr>   <td> " << s << " </td> <td>" << spos << "</td> <td> " << vFlatTop*100.0 << " </td>  <td> " << vAccel << " </td>  <td> <a href=\"javascript: showResult('" 
	       << pngfn << "');\"> <img border=0  width=400 src='ImpactT.mpeg' ALT='link' TITLE='link'></a> </td> </tr>" << endl;
	}   
	ostr << "</table> </div>" << endl;
      }
      ostr.close();
    }
  }
  
  
  if(doFWHM) {

    if (c1)
      delete c1;
    
    c1 = new TCanvas("c1",baseDir,1,1);

    gStyle->SetMarkerSize(0.1);

    TString globtit = baseDir;
    TString myname = "Dude";

    c1->Update();
    
    cout << endl << " *********************** doFWHM *****************************" << endl;
    cout << "   ********* Analysis of " << baseDir << " ********" << endl;
    cout << "   ********* Data samples from " << startCase << " to " << endCase <<   " ********" << endl;
    cout << " *********************************************************" << endl;

   
    for (Int_t d=startCase; d<endCase; d++) {
      TString mynum("");	
      mynum+=d;	
      fn[d]    = new TString(baseDir+"/Ring-"+d+"/Ring-"+d+".h5");
      cout << *fn[d] << endl;	
    }
    
    for (Int_t d=startCase; d<endCase; d++) {
      TString mynum("");
      if (data)
	delete data;
      data = new TH5Dataset(*fn[d],subsampling);    
     
      // negative endSteps means all data
      if (endStep<0)
	endStep = data->GetNSteps();       //      endStep = 20;

      cout << "Read in " << endStep << " steps process data from " << startStep << " to " << endStep << endl;
      Double_t x[6];
      nt[d] = new TNtuple("ntuple","data from ascii =file","s:fwhmx:rmsx:meanx:minx:maxx");
      for (Int_t s = startStep; s<endStep; s++) {

	// x coordnate

 	TH1D* h1 = data->Histo1d("x",s,-0.04,0.04,256);         
	x[2] = h1->GetRMS();
	x[3] = h1->GetMean();

	x[4] = h1->GetMinimum();
	x[5] = h1->GetMaximum();

	Fwhm(h1,peak,fwhm);
	x[0] = data->GetSpos(s);
	x[1] = fwhm;

	// z coordnate

 	TH1D* h2 = data->Histo1d("z",s,-0.1,0.1,256);         
	x[2] = h2->GetRMS();
	x[3] = h2->GetMean();

	x[4] = h2->GetMinimum();
	x[5] = h2->GetMaximum();

	Fwhm(h2,peak,fwhm);
	x[1] = fwhm;


	nt[d]->Fill(x[0],x[1],x[2],x[3],x[4],x[5]);
	
	if (maxVal[0] < x[0])
	  maxVal[0]= x[0];

	for (Int_t i=1; i<=5;i++) {
	  if (maxVal[i] < x[i])
	    maxVal[i]= x[i];
	  else if (minVal[i] < x[i])
	    minVal[i]= x[i];
	}
	cout << "Work on step " << s << " out of " << endStep << " FWHM " << fwhm << endl;	

	if (h1D)
	  delete h1D;
      }   
    }
    
    if (c1)
      delete c1;
    
    drawFwhm(baseDir,maxVal,globtit,nt,startCase,endCase);
    

  }

    //TLatex  t;	
    //t.SetTextSize(0.03);
    //t.DrawLatex(0.0+0.1*maxS,0.1*maxZ,Form("Data set %d Q= %g [Cb] I_0= %g [A] N= %g",d+1,Q,I0,N)); 



 
  if(doStat) {

    gStyle->SetMarkerSize(0.1);

    TString globtit = baseDir;
    TString myname = "Dude";

    TPaveLabel* globTitle = new TPaveLabel(0.1,0.95,0.9,0.99,globtit);
    globTitle->SetBorderSize(0);
    globTitle->SetFillColor(0);
    globTitle->SetTextFont(62);
    globTitle->Draw();
  
TDatime now;
TPaveLabel* date = new TPaveLabel(0.7,0.01,0.9,0.05,now.AsString());
date->SetBorderSize(0);
date->SetFillColor(0);
date->SetTextFont(62);
date->SetTextSize(0.3);
date->Draw();

TPad* graphs = new TPad("c1Graphs","Graphs",0.01,0.05,0.95,0.95);
graphs->SetBottomMargin(4.);
graphs->Draw();
graphs->Divide(1,2);

graphs->cd(1);gPad->SetRightMargin(0.1);gPad->SetBottomMargin(0.2);
graphs->cd(2);gPad->SetRightMargin(0.1);gPad->SetBottomMargin(0.2);

gStyle->SetTickLength(0.03,"x");
gStyle->SetTickLength(0.01,"y");
gStyle->SetPadTickX(0);
gStyle->SetPadTickY(0);
   
c1->Update();


  /* 
     Make up the filenames we need	

  */

  for (Int_t d=startCase; d<=endCase; d++) {
    TString mynum("");	
    mynum+=d;	
    fn[d] = new TString(baseDir+"/Ring-"+mynum+"/Ring-"+mynum);
    cout << *fn[d] << endl;	
  }

  cout << " *********************************************************" << endl;
  cout << "   ********* Analysis of " << baseDir << " ********" << endl;
  cout << "   ********* Data samples from " << startCase << " to " << endCase <<   " ********" << endl;
  cout << " *********************************************************" << endl;

  for (Int_t d=startCase; d<=endCase; d++) {
    TString infn = *fn[d] + TString(".stat");
    cout << infn << endl;
    nt[d] = new TNtuple("ntuple","data from ascii =file","s:phi:en:x:y:z:px:py:pz:ex:ey:ez");

    ifstream in(infn);

    Int_t nL = 34;
    Char_t ch = in.get();
    for (int i=0;i<nL;i++) {
      while (ch != '\n') 
	ch = in.get();
      ch = in.get();
    }
    in >> Q ;
    in >> I0;
    in >> N;
    in >> nP;
    cout << "Processing " << *fn[d] << " Q= " << Q << " N= " << N << " I0= " << I0 << " NP= " << nP << endl;

    Double_t x[26+nP];	
    Int_t counter=0;  
    if (doStatTest) {
      while (counter<1000) {  
	counter++;
	for (Int_t l=0;l<26+nP;l++) {
	  in >> x[l];
	}
	if (!in.good()) break;
	nt[d]->Fill(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11]);
	for (Int_t i=0; i<DATATOREAD; i++) {
	  if (maxVal[i] < x[i])
	    maxVal[i]= x[i];
	  else 
	    if (minVal[i] > x[i])
	      minVal[i]= x[i];
	}
      }
      
    } else {
      while (1) { 
	for (Int_t l=0;l<26+nP;l++) {
	  in >> x[l];
	}
	if (!in.good()) break;
	nt[d]->Fill(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11]);
	for (Int_t i=0; i<DATATOREAD; i++) {
	  if (maxVal[i] < x[i])
	    maxVal[i]= x[i];
	  else 
	    if (minVal[i] > x[i])
	      minVal[i]= x[i];
	}
      }
    }
    in.close();
  }

  graphs->cd(1);

  TH1F* frame;

  TLegend *legend=new TLegend(0.7,0.65,0.88,0.85);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->SetMargin(0.2);
  
  for (Int_t d=startCase; d<=endCase; d++) {
    if (d==(startIndex)) {
      frame = c1->DrawFrame(0,0,1.1*maxVal[0],1.1*maxVal[5]);
      frame->SetXTitle("s [m] ");
      frame->SetYTitle("#sigma_{z} [m] ");
      //TLatex  t;	
      //t.SetTextSize(0.03);
      //t.DrawLatex(0.0+0.1*maxS,0.1*maxZ,Form("Data set %d Q= %g [Cb] I_0= %g [A] N= %g",d+1,Q,I0,N)); 
    }
    nt[d]->Draw("z:s","","goff");
    gr[d] = new TGraph(nt[d]->GetSelectedRows(),nt[d]->GetV2(),nt[d]->GetV1());
    gr[d]->SetLineColor(d);
    gr[d]->Draw("Lsame"); //draw graph in current pad
    legend->AddEntry(gr[d],Form("case %d",d),"l");
  }
  legend->Draw();

  graphs->cd(2);

  for (Int_t d=startCase; d<=endCase; d++) {
    if (d==(startIndex)) {
      frame = c1->DrawFrame(0,0,1.1*maxVal[0],1.1*maxVal[3]);
      frame->SetXTitle("s [m] ");
      frame->SetYTitle("#sigma_{x} [m] ");
    }
    nt[d]->SetMarkerColor(d);
    nt[d]->Draw("x:s","","same");
  }
  //  legend->Draw();

  TString ofn = baseDir+TString("-plot1")+TString(".eps");
  cout << "Save " << ofn << endl; 

  c1->Update();
  c1->Print(ofn);
  c1->Update();
}

  if (doPart) {
    Int_t picperpage = 16;
    TH2F* xzh[16];
    /* 
       Make up the filenames we need	
    */
    for (Int_t d=startCase; d<endCase; d++) {
      TString mynum("");	
      mynum+=d;	
      fn[d]    = new TString(baseDir+"/Ring-"+d+"/Ring-"+d+".h5");
      fnwww[d] = new TString(baseDir+"/Ring-"+d+"/Ring-"+d+".html");
      cout << *fn[d] << endl;	
    }

    TString globtit = baseDir;
    TString myname = "Dude";

    TPaveLabel* globTitle = new TPaveLabel(0.1,0.95,0.9,0.99,globtit);
    globTitle->SetBorderSize(0);
    globTitle->SetFillColor(0);
    globTitle->SetTextFont(62);

    TDatime now;
    TPaveLabel* date = new TPaveLabel(0.7,0.01,0.9,0.05,now.AsString());
    date->SetBorderSize(0);
    date->SetFillColor(0);
    date->SetTextFont(62);
    date->SetTextSize(0.3);

    Int_t pstype = 111;
    TString psfn = baseDir+TString("-xz")+TString(".ps");
    TPostScript *ps = new TPostScript(psfn,pstype);

    TLatex l; 
    l.SetNDC();
    l.SetTextSize(0.04);
    l.SetTextColor(2);

    date->Draw();
    globTitle->Draw(); 

    TPad* graphs = new TPad("c1Graphs","Graphs",0.01,0.05,0.95,0.95);
    graphs->SetBottomMargin(4.);
    graphs->Draw();
    graphs->Divide(4,4);
    for (Int_t p = 1; p<=picperpage; p++) {    
      graphs->cd(p);
      gPad->SetRightMargin(0.1);
      gPad->SetBottomMargin(0.2);
    }
      
    gStyle->SetTickLength(0.03,"x");
    gStyle->SetTickLength(0.01,"y");
    gStyle->SetPadTickX(0);
    gStyle->SetPadTickY(0);
    c1->Update();

    gStyle->SetMarkerSize(0.22);
    gStyle->SetMarkerColor(4);

    for (Int_t d=startCase; d<endCase; d++) {

      data = new TH5Dataset(*fn[d],subsampling);    
      Int_t nstep = data->GetNSteps();

      TString casenum(""); casenum += d;
      
      // negative endSteps means all data
      if (endStep<0)
	endStep = nstep;

      cout << "Read in " << nstep << " steps process data from " << startStep << " to " << endStep << endl;

      for (Int_t s = startStep; s<endStep; ) {
	ps->NewPage();           
	for (Int_t p = 1; p<=picperpage; p++) {
	  graphs->cd(p);
	  xzh[p-1] = data->Histo2d("x","z",s,-0.04,0.04,-0.2,0.2);         
	  TString mynum(""); mynum += s;
	  TString s1 = TString("Turn ")+mynum+" case"+casenum;
	  l.DrawLatex(0.15,0.85,s1); 
	  s++;
	}       
	c1->Update();
      }
    } 
    ps->Close();
  }
  return 0;
}

/*

	c1->cd(2);
	TH2F* xpxh = data.Histo2d("x","px",s,-0.04,0.04,-0.02,0.02);  
	TString s2("s= ");
	s2 += data.GetSpos(s);
	s2 += " [m]"; 
	l.DrawLatex(0.15,0.015,s2);
       
	c1->cd(3);
	TH2F* ypyh = data.Histo2d("y","py",s,-0.04,0.04,-0.02,0.02);  
	c1->cd(4);
	TH2F* zpzh = data.Histo2d("z","pz",s,-0.2,0.2,-0.001,0.001);  
	c1->cd(5);
	gPad->SetLogy();
	TH1D* xh = data.Histo1d("x",s,-0.04,0.04);  
	xh->SetFillColor(46);
	c1->cd(6);
	gPad->SetLogy();
	TH1D* zh = data.Histo1d("z",s,-0.2,0.2);  
	zh->SetFillColor(36);
       


*/
