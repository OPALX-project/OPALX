TCanvas* decorate4Canvases(TString globtit)
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
  
  TCanvas* graphs = new TCanvas("c1Graphs","Graphs",800,1000);
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


