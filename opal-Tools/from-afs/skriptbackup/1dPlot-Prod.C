{
gSystem->Load("~schietinger/bin/H5root.so");
TH5Style::SetStyle();
TCanvas* mydummy = new TCanvas("mydummy","  ",1,1);
mydummy->cd();

const Int_t dataSize = 2;
const Int_t Scenarios = 1;
TString cases[dataSize];
TString scenarios[Scenarios];
TString scenariosDesc[dataSize];

TH5Dataset *data[Scenarios][dataSize];

// **********************************************

cases[0] = TString("1");  
cases[1] = TString("2");

scenarios[0] = TString("14");

scenariosDesc[0] = TString("V_{accel}=0.7 [MV], V_{flatTop}=100 % gauss");
scenariosDesc[1] = TString("V_{accel}=0.7 [MV], V_{flatTop}=100 % parabolic");


// **********************************************

Double_t xd[dataSize];
Double_t yd[dataSize];

TGraph  *gr[dataSize];

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
TH1F*  myFrame = myc->DrawFrame(0.005,0.1,0.1,2);
myFrame->SetXTitle("#sigma_{z-initial}");
myFrame->SetYTitle("#sigma_{z}/#sigma_{z-initial}");
TLegend *legend=new TLegend(0.2,0.73,0.68,0.88);    
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

    cout << "zini-rms= " << hzi->GetRMS() << " zfinal-rms= " << hzf->GetRMS() << endl;

    xd[c] = hzi->GetRMS();     
    yd[c] = hzf->GetRMS()/hzi->GetRMS(); 
    myc->cd();
    gr[c] = new TGraph(dataSize,xd,yd);
    //    gr[c]->SetMarkerColor(c+1);
    gr[c]->SetMarkerSize(1);
    gr[c]->SetMarkerStyle(21+c);
    gr[c]->Draw("psame");
    legend->AddEntry(gr[c],scenariosDesc[c],"p");
    mydummy->cd();  
  }
}
myc->cd();
legend->Draw();
TString ofn = TString("sigma-z-ratios")+TString(".prod.png");
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
TLegend *legend=new TLegend(0.2,0.73,0.68,0.88);    
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

    xd[c] = hzi->GetRMS();     
    yd[c] = hzf->GetMaximum()/hzi->GetMaximum(); 
 
    cout << "zini= " << hzi->GetRMS() << " zmax= " << hzf->GetMaximum() << " zmin= " << hzi->GetMaximum() << endl;
    myc->cd();
    gr[c] = new TGraph(dataSize,xd,yd);
    //    gr[c]->SetMarkerColor(c+1);
    gr[c]->SetMarkerSize(1);
    gr[c]->SetMarkerStyle(21+c);
    gr[c]->Draw("psame");
    legend->AddEntry(gr[c],scenariosDesc[c],"p");
    mydummy->cd();  
  }
}
myc->cd();
legend->Draw();
TString ofn = TString("max-z-ratios")+TString(".prod.png");
myc->Print(ofn);
myc->Update();



// ===============================================================================
// #sigma_{x}/#sigma_{z-initial}
// ===============================================================================

TCanvas* myc = new TCanvas("myc","  ",800,600);

myc->cd();
TH1F*  myFrame = myc->DrawFrame(0.005,0.1,0.1,2);
myFrame->SetXTitle("#sigma_{z-initial}");
myFrame->SetYTitle("#sigma_{x}/#sigma_{z-initial}");
TLegend *legend=new TLegend(0.2,0.73,0.68,0.88);    
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
    TH1D* hzf = (*data[scenario][c]).Histo1d("x",s,-0.13,0.13,256);         

    xd[c] = hzi->GetRMS();     
    yd[c] = hzf->GetRMS()/hzi->GetRMS(); 
    myc->cd();
    gr[c] = new TGraph(dataSize,xd,yd);
    //    gr[c]->SetMarkerColor(c+1);
    gr[c]->SetMarkerSize(1);
    gr[c]->SetMarkerStyle(21+c);
    gr[c]->Draw("psame");
    legend->AddEntry(gr[c],scenariosDesc[c],"p");
    mydummy->cd();  
  }
}
myc->cd();
legend->Draw();
TString ofn = TString("sigma-x-ratios")+TString(".prod.png");
myc->Print(ofn);
myc->Update();



// ===============================================================================
// #max_{x}/#max_{z-initial}
// ===============================================================================

TCanvas* myc = new TCanvas("myc","  ",800,600);
myc->cd();
TH1F*  myFrame = myc->DrawFrame(0.001,1,0.03,4);
myFrame->SetXTitle("#sigma_{z-initial}");
myFrame->SetYTitle("x_{max}/z-initial");
TLegend *legend=new TLegend(0.2,0.73,0.68,0.88);    
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
    TH1D* hzf = (*data[scenario][c]).Histo1d("x",s,-0.13,0.13,256);         

    xd[c] = hzi->GetRMS();     
    yd[c] = hzf->GetMaximum()/hzi->GetMaximum(); 
    myc->cd();
    gr[c] = new TGraph(dataSize,xd,yd);
    //    gr[c]->SetMarkerColor(c+1);
    gr[c]->SetMarkerSize(1);
    gr[c]->SetMarkerStyle(21+c);
    gr[c]->Draw("psame");
    legend->AddEntry(gr[c],scenariosDesc[c],"p");
    mydummy->cd();  
  }
}
myc->cd();
legend->Draw();
TString ofn = TString("max-x-ratios")+TString(".prod.png");
myc->Print(ofn);
myc->Update();

}

