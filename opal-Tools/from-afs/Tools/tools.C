void processFortFile(TCanvas *c1, TString dir, TString fn) {

 {
   cout << " *********************************************************" << endl;
   cout << "   ********* Analysis of fort.22 ********" << endl;
   cout << " *********************************************************" << endl;

   TString baseFn=dir+"/"+fn+"/"+fn;
   TString infn = dir + "/" + fn + TString("/fort.22");
   cout << "Read from " << infn << endl; 
   TNtuple *nt = new TNtuple("ntuple","data from ascii =file","x:px:y:py:t:pz:z:en:bin");

   ifstream in(infn);

   Double_t x[9];	

   Int_t counter=0;  
    //  while (counter < 10000) {
    while (1) {                     // 
     for (Int_t l=0;l<9;l++) {
       in >> x[l];
     }
     if (!in.good()) break;
     nt->Fill(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]);
     counter++;
    }
    in.close();
    cout << "read in of " << counter << " particles" << endl; 
   }
  
 
   TH1F *h1 = new TH1F("h1","en", 100,nt->GetMinimum("en"),nt->GetMaximum("en"));
   h1->SetFillColor(46);
   h1->GetXaxis()->SetTitle("E [MeV]");
   nt->Draw("en>>h1");
   c1->Print(baseFn + "-enatscreen.png");   
   delete h1;
  
   TH2F *h2 = new TH2F("h2","t pz",100, nt->GetMinimum("t"), nt->GetMaximum("t"),
		                   100, nt->GetMinimum("pz"),nt->GetMaximum("pz"));
   h2->GetXaxis()->SetTitle("t [s]");
   h2->GetYaxis()->SetTitle("#beta_{z}#gamma");

   h2->SetMarkerSize(0.25);
   nt->Draw("pz:t>>h2");  
   c1->Print(baseFn+"-tpzatscreen.png");   
   delete h2;


   h1 = new TH1F("h1","t", 256,nt->GetMinimum("t"),nt->GetMaximum("t"));
   h1->SetFillColor(46);
   nt->Draw("t>>h1");
   c1->Print(baseFn+"-tatscreen.png");   
   delete h1;

   TH2F *h2 = new TH2F("h2","x px",100, nt->GetMinimum("x"), nt->GetMaximum("x"),
		                   100, nt->GetMinimum("px"),nt->GetMaximum("px"));
   h2->GetXaxis()->SetTitle("x [m]");
   h2->GetYaxis()->SetTitle("#beta_{x}#gamma");
   gStyle->SetStatY(0.414);
   nt->Draw("px:x>>h2");  
   h2->SetMarkerSize(0.25);
   c1->Print(baseFn+"-xpxatscreen.png");   
   delete h2;

   h1 = new TH1F("h1","x", 256,nt->GetMinimum("x"),nt->GetMaximum("x"));
   h1->SetFillColor(46);
   h1->GetXaxis()->SetTitle("x [m]");
   nt->Draw("x>>h1");
   c1->Print(baseFn+"-xatscreen.png");   
   delete h1;

   h1 = new TH1F("h1","px", 256,nt->GetMinimum("px"),nt->GetMaximum("px"));
   h1->SetFillColor(46);
   h1->GetXaxis()->SetTitle("#beta_{x}#gamma");
   nt->Draw("px>>h1");
   c1->Print(baseFn+"-pxatscreen.png");   

  }

void processMafiaFile(TCanvas *c1, TString dir, TString fn) {

 {
   cout << " *********************************************************" << endl;
   cout << "   ********* Analysis of Mafia data ********" << endl;
   cout << " *********************************************************" << endl;

   TString baseFn=dir+"/"+fn;
   TString infn = dir + "/" + fn;
   cout << "Read from " << infn << endl; 
   TNtuple *nt = new TNtuple("ntuple","data from ascii =file","x:y:z:px:py:pz:t");

   ifstream in(infn);

   Double_t x[7];	

   Int_t counter=0;  
    //  while (counter < 10000) {
    while (1) {                     // 
     for (Int_t l=0;l<7;l++) {
       in >> x[l];
     }
     if (!in.good()) break;
     nt->Fill(x[0],x[1],x[2],x[3],x[4],x[5],x[6]);
     counter++;
    }
    in.close();
    cout << "read in of " << counter << " particles" << endl; 
   }
  
 
   TH2F *h2 = new TH2F("h2","t pz",100, nt->GetMinimum("t"), nt->GetMaximum("t"),
		                   100, nt->GetMinimum("pz"),nt->GetMaximum("pz"));
   h2->GetXaxis()->SetTitle("t [s]");
   h2->GetYaxis()->SetTitle("#beta_{z}#gamma");
   h2->SetMarkerSize(0.25);
   nt->Draw("pz:t>>h2");  
   c1->Print(baseFn+"-tpzatscreen.png");   
   delete h2;

     
   TH1F *h1 = new TH1F("h1","t", 256,nt->GetMinimum("t"),nt->GetMaximum("t"));
   h1->SetFillColor(46);
   nt->Draw("t>>h1");
   c1->Print(baseFn+"-tatscreen.png");   
   delete h1;

   h2 = new TH2F("h2","x px",100, nt->GetMinimum("x"), nt->GetMaximum("x"),		                   100, nt->GetMinimum("px"),nt->GetMaximum("px"));
   h2->GetXaxis()->SetTitle("x [m]");
   h2->GetYaxis()->SetTitle("#beta_{x}#gamma");
   gStyle->SetStatY(0.414);
   nt->Draw("px:x>>h2");  
   h2->SetMarkerSize(0.25);
   c1->Print(baseFn+"-xpxatscreen.png");   
   delete h2;

   h1 = new TH1F("h1","x", 256,nt->GetMinimum("x"),nt->GetMaximum("x"));
   h1->SetFillColor(46);
   h1->GetXaxis()->SetTitle("x [m]");
   nt->Draw("x>>h1");
   c1->Print(baseFn+"-xatscreen.png");   
   delete h1;

   h1 = new TH1F("h1","px", 256,nt->GetMinimum("px"),nt->GetMaximum("px"));
   h1->SetFillColor(46);
   h1->GetXaxis()->SetTitle("#beta_{x}#gamma");
   nt->Draw("px>>h1");
   c1->Print(baseFn+"-pxatscreen.png");   

  }


TGraph * processH5File(TCanvas *c1, TString dir, TString fn, TString fext) {

   TString baseFn=dir+"/"+fn+"/"+fn;

   TH5Dataset data(dir+"/"+fn+TString("/ImpactT.h5"),50);     

   TGraph* emit = data.GetGraphXY("SPOS","#varepsilon","x");
   Double_t smin = 1.1*TMath::MinElement(emit->GetN(),emit->GetX());
   Double_t smax = 1.1*TMath::MaxElement(emit->GetN(),emit->GetX());

   Double_t emitmin = 1.1*TMath::MinElement(emit->GetN(),emit->GetY());
   Double_t emitmax = 1.1*TMath::MaxElement(emit->GetN(),emit->GetY());

   TH1F *frame = c1->DrawFrame(smin,emitmin,smax,emitmax);  
   gPad->SetRightMargin(0.1);
   frame->SetXTitle("s [m]");  
   frame->SetYTitle("#varepsilon_{nx} [m rad] projected ");  

   TLatex l; 
   l.SetTextSize(0.04);
   TString s1 = "Slice option: mean and 1 mm"; 
   TString s2 =  fn; 
   l.DrawLatex(0.1*smax,0.93*emitmax,s1);
   l.DrawLatex(0.1*smax,1.02*emitmax,s2);

   emit->Draw("LPsame");  

   data.SetSlicePosition("mean");
   data.SetSliceLength(1.0,"mm");

   TGraph* slemit = data.GetGraphXY("SPOS","emittance(x,px)","slice");


   smin = 1.1*TMath::MinElement(slemit->GetN(),slemit->GetX());
   smax = 1.1*TMath::MaxElement(slemit->GetN(),slemit->GetX());
   emitmin = 1.1*TMath::MinElement(slemit->GetN(),slemit->GetY());
   emitmax = 1.1*TMath::MaxElement(slemit->GetN(),slemit->GetY());


   // second axis to the right:
   TGaxis *axis = new TGaxis(smax,emitmin,smax,emitmax,emitmin,emitmax, 506,"SL+");
   axis->SetTitle("#varepsilon_{nx} [m rad] slice");  
   axis->SetTitleSize(frame->GetTitleSize());
   axis->SetTitleOffset(.8);
   axis->SetLabelSize(frame->GetLabelSize());
   // axis->SetLabelOffset(0.
   //   axis->SetTickSize(frame->GetTickSize());
   axis->Draw();
   slemit->SetMarkerColor(2);
   slemit->Draw("LPsame");  

   c1->Print(baseFn+"-emitx."+fext);   

  }



