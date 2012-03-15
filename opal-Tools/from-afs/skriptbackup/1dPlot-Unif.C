{
gSystem->Load("~schietinger/bin/H5root.so");
TH5Style::SetStyle();
TCanvas* mydummy = new TCanvas("mydummy","  ",1,1);
mydummy->cd();

const Int_t dataSize = 4;
const Int_t Scenarios = 4;
TString cases[dataSize];
TString scenarios[Scenarios];
TString scenariosDesc[Scenarios];

TH5Dataset *data[Scenarios][dataSize];

// **********************************************

cases[0] = TString("1");
cases[1] = TString("2");
cases[2] = TString("3");
cases[3] = TString("5");

scenarios[0] = TString("6");
scenarios[1] = TString("7");
scenarios[2] = TString("8");
scenarios[3] = TString("9");

scenariosDesc[0] = TString("V_{accel}=0.9 [MV], V_{flatTop}=100 %");
scenariosDesc[1] = TString("V_{accel}=0.9 [MV], V_{flatTop}=80 %");
scenariosDesc[2] = TString("V_{accel}=1.0 [MV], V_{flatTop}=100 %");
scenariosDesc[3] = TString("V_{accel}=1.0 [MV], V_{flatTop}=80 %");

// **********************************************

Double_t xd[Scenarios][dataSize];
Double_t yd[Scenarios][dataSize];

TGraph  *gr[Scenarios];

Int_t subsampling = 1;

for (Int_t scenario=0; scenario < Scenarios; scenario++) {
  TString globtit(TString("Ring-FlatTopScenario-" + scenarios[scenario]));
  for (Int_t c=0; c<dataSize; c++) {
    TString fn(globtit + "/case" + cases[c] + "/Ring-" + cases[c] + ".h5");
    data[scenario][c] = new TH5Dataset(fn,subsampling);
  }
}


// ===============================================================================
// #sigma_{z}/#sigma_{z-initial}
// ===============================================================================

TCanvas* myc = new TCanvas("myc","  ",800,600);
myc->cd();
TH1F*  myFrame = myc->DrawFrame(0.001,1,0.03,4);
myFrame->SetXTitle("#sigma_{z-initial}");
myFrame->SetYTitle("#sigma_{z}/#sigma_{z-initial}");
TLegend *legend=new TLegend(0.3,0.73,0.68,0.88);    
legend->SetTextSize(0.03);
legend->SetFillColor(0);
legend->SetMargin(0.2);
mydummy->cd();  

for (Int_t scenario=0; scenario < Scenarios; scenario++) {
  TString globtit(TString("Ring-FlatTopScenario-" + scenarios[scenario]));
  for (Int_t c=0; c<dataSize; c++) {
    TString fn(globtit + "/case" + cases[c] + "/Ring-" + cases[c] + ".h5");
    cout <<  globtit << " working on case " << cases[c] << endl;

    Int_t s = 0;
    TH1D* hzi = (*data[scenario][c]).Histo1d("z",s,-0.13,0.13,256);         
  
    s = (*data[scenario][c]).GetNSteps()-1;
    TH1D* hzf = (*data[scenario][c]).Histo1d("z",s,-0.13,0.13,256);         

    xd[scenario][c] = hzi->GetRMS();     
    yd[scenario][c] = hzf->GetRMS()/hzi->GetRMS(); 
    gr[scenario] = new TGraph(dataSize,xd[scenario],yd[scenario]);
  }
  myc->cd();
  legend->AddEntry(gr[scenario],scenariosDesc[scenario],"p");
  gr[scenario]->SetMarkerColor(scenario+1);
  gr[scenario]->SetMarkerSize(1);
  gr[scenario]->SetMarkerStyle(21);
  gr[scenario]->Draw("psame");
  mydummy->cd();  
}
myc->cd();
legend->Draw();
TString ofn = TString("sigma-z-ratios")+TString(".unif.pdf");
myc->Print(ofn);
myc->Update();


// ===============================================================================
// #max_{z}/#max_{z-initial}
// ===============================================================================

TCanvas* myc = new TCanvas("myc","  ",800,600);
myc->cd();
TH1F*  myFrame = myc->DrawFrame(0.001,1,0.03,4);
myFrame->SetXTitle("#sigma_{z-initial}");
myFrame->SetYTitle("z_{max}/z-initial_{max}");
TLegend *legend=new TLegend(0.3,0.73,0.68,0.88);    
legend->SetTextSize(0.03);
legend->SetFillColor(0);
legend->SetMargin(0.2);
mydummy->cd();  

for (Int_t scenario=0; scenario < Scenarios; scenario++) {
  TString globtit(TString("Ring-FlatTopScenario-" + scenarios[scenario]));
  for (Int_t c=0; c<dataSize; c++) {
    TString fn(globtit + "/case" + cases[c] + "/Ring-" + cases[c] + ".h5");
    cout <<  globtit << " working on case " << cases[c] << endl;
    
    Int_t s = 0;
    TH1D* hzi = (*data[scenario][c]).Histo1d("z",s,-0.13,0.13,256);         
  
    s = (*data[scenario][c]).GetNSteps()-1;
    TH1D* hzf = (*data[scenario][c]).Histo1d("z",s,-0.13,0.13,256);         

    xd[scenario][c] = hzi->GetRMS();     
    yd[scenario][c] = hzf->GetMaximum()/hzi->GetMaximum(); 
    gr[scenario] = new TGraph(dataSize,xd[scenario],yd[scenario]);
  }
  myc->cd();
  legend->AddEntry(gr[scenario],scenariosDesc[scenario],"p");
  gr[scenario]->SetMarkerColor(scenario+1);
  gr[scenario]->SetMarkerSize(1);
  gr[scenario]->SetMarkerStyle(21);
  gr[scenario]->Draw("psame");
  mydummy->cd();  
}
myc->cd();
legend->Draw();
TString ofn = TString("max-z-ratios")+TString(".unif.pdf");
myc->Print(ofn);
myc->Update();



// ===============================================================================
// #sigma_{x}/#sigma_{x-initial}
// ===============================================================================

TCanvas* myc = new TCanvas("myc","  ",800,600);
myc->cd();
TH1F*  myFrame = myc->DrawFrame(0.001,0,0.03,4);
myFrame->SetXTitle("#sigma_{z-initial}");
myFrame->SetYTitle("#sigma_{x}/#sigma_{x-initial}");
TLegend *legend=new TLegend(0.3,0.73,0.68,0.88);    
legend->SetTextSize(0.03);
legend->SetFillColor(0);
legend->SetMargin(0.2);
mydummy->cd();  

for (Int_t scenario=0; scenario < Scenarios; scenario++) {
  TString globtit(TString("Ring-FlatTopScenario-" + scenarios[scenario]));
  for (Int_t c=0; c<dataSize; c++) {
    TString fn(globtit + "/case" + cases[c] + "/Ring-" + cases[c] + ".h5");
    cout <<  globtit << " working on case " << cases[c] << endl;

    Int_t s = 0;
    TH1D* hxi = (*data[scenario][c]).Histo1d("x",s,-0.03,0.03,256);         
    TH1D* hzi = (*data[scenario][c]).Histo1d("z",s,-0.13,0.13,256);         
  
    s = (*data[scenario][c]).GetNSteps()-1;
    TH1D* hxf = (*data[scenario][c]).Histo1d("x",s,-0.03,0.03,256);         

    xd[scenario][c] = hzi->GetRMS();     
    yd[scenario][c] = hxf->GetRMS()/hxi->GetRMS();
    cout << hxf->GetRMS() << "  " << hxi->GetRMS() << endl;
    gr[scenario] = new TGraph(dataSize,xd[scenario],yd[scenario]);
  }
  myc->cd();
  legend->AddEntry(gr[scenario],scenariosDesc[scenario],"p");
  gr[scenario]->SetMarkerColor(scenario+1);
  gr[scenario]->SetMarkerSize(1);
  gr[scenario]->SetMarkerStyle(21);
  gr[scenario]->Draw("psame");
  mydummy->cd();  
}
myc->cd();
legend->Draw();
TString ofn = TString("sigma-x-ratios")+TString(".unif.pdf");
myc->Print(ofn);
myc->Update();


// ===============================================================================
// #max_{x}/#max_{x-initial}
// ===============================================================================

TCanvas* myc = new TCanvas("myc","  ",800,600);
myc->cd();
TH1F*  myFrame = myc->DrawFrame(0.001,1,0.03,4);
myFrame->SetXTitle("#sigma_{z-initial}");
myFrame->SetYTitle("x_{max}/x-initial_{max}");
TLegend *legend=new TLegend(0.3,0.73,0.68,0.88);    
legend->SetTextSize(0.03);
legend->SetFillColor(0);
legend->SetMargin(0.2);
mydummy->cd();  

for (Int_t scenario=0; scenario < Scenarios; scenario++) {
  TString globtit(TString("Ring-FlatTopScenario-" + scenarios[scenario]));
  for (Int_t c=0; c<dataSize; c++) {
    TString fn(globtit + "/case" + cases[c] + "/Ring-" + cases[c] + ".h5");
    cout <<  globtit << " working on case " << cases[c] << endl;
    
    Int_t s = 0;
    TH1D* hxi = (*data[scenario][c]).Histo1d("x",s,-0.03,0.03,256);         
    TH1D* hzi = (*data[scenario][c]).Histo1d("z",s,-0.13,0.13,256);         
  
    s = (*data[scenario][c]).GetNSteps()-1;
    TH1D* hxf = (*data[scenario][c]).Histo1d("x",s,-0.03,0.03,256);         

    xd[scenario][c] = hzi->GetRMS();     
    yd[scenario][c] = hxf->GetMaximum()/hxi->GetMaximum(); 
    gr[scenario] = new TGraph(dataSize,xd[scenario],yd[scenario]);
  }
  myc->cd();
  legend->AddEntry(gr[scenario],scenariosDesc[scenario],"p");
  gr[scenario]->SetMarkerColor(scenario+1);
  gr[scenario]->SetMarkerSize(1);
  gr[scenario]->SetMarkerStyle(21);
  gr[scenario]->Draw("psame");
  mydummy->cd();  
}
myc->cd();
legend->Draw();
TString ofn = TString("max-x-ratios")+TString(".unif.pdf");
myc->Print(ofn);
myc->Update();
}




// draw a frame to define the range
// TMultiGraph *mg = new TMultiGraph();
//mg->Add(gr[0]); 
//mg->Draw("p");

