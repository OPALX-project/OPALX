// SPECIAL VERSION FOR LINAC10
// Macro producing plots from OPAL simulation data
// Execute in a root session like this:
//
//   .x plotFinPhase1.C
//
// To get the shared library, you must compile H5PartROOT with 'gmake shlib'.
// Be sure that the underlying H5Part has been compiled with the -fPIC option.
// If not, reconfigure your H5Part like this (bash shell):
//   cd $H5Part
//   ./autogen.sh
//   CFLAGS = -fPIC ./config
//   make
//
{
   /////////////////////////////////////////////////////////////////////
   // USER INPUTS:

   // the simulation file to analyse:
   //  TString simPath("/scratch2/amas/schietinger/FinOPAL/FIN_Phase1_GUNSOLB=214/");

   TString simPath("FinPhase1_FIND1_MSOL10_ks=0.164_QBUNCH=100_SIGX=300_SIGY=300_MX=32_MY=32_MT=64_NPART=2000000/");
   TString simFilename("FinPhase1.h5");

   // the measurement files:
   TString measPath("/home2/adelmann/scratch3/250MeVInjector/Linac2010/");
   TString rmsFilename("rmsData23June2010.dat");
   TString emitFilename("emitData24June2010.dat");

   // switch for slice emittance (takes long time!)
   const Bool_t doSlice = kTRUE;

   // plot dimensions
   const Double_t sMin = 0.0;
   const Double_t sMax = 5.0;
   const Double_t EmittMax = 2.8;
   const Double_t RmsMax = 4.;
   const Double_t EneMax = 8.; // in MeV
   const Double_t dEMax = 60.; // in keV
   const Double_t Bmin = -0.25;
   const Double_t Bmax = 0.25;
   const Double_t Emin=-150.;
   const Double_t Emax=150.;

   // scale factors
   const Double_t emittScale = 1.e6;
   const Double_t rmsScale = 1.e3*EmittMax/RmsMax;
   const Double_t dEScale = 1.e3*EneMax/dEMax; // factor 1.e3 for keV

   // load the shared library:
   gSystem->Load("/home5/schietinger/extlib/H5Root/lib/libh5root.so");
   using namespace TH5Util;

   // let's do it in style:
   TH5Style::SetStyle();
   gStyle->SetPadTickY(0);
   Width_t lineWidth(3);

   const Double_t titSize = 0.07;
   const Double_t labSize = 0.06;
   const Double_t titOffset = 0.5;

   gStyle->SetTitleSize(titSize,"X");
   gStyle->SetTitleSize(titSize,"Y");
   gStyle->SetLabelSize(labSize,"X");
   gStyle->SetLabelSize(labSize,"Y");

   // load the data file:
   TH5Dataset data(simPath+simFilename);

   Int_t nStep = data.GetNSteps();

   // retrieve (kinetic!) energy (for normalized emittances)
   // NOTE: "#gamma" is not relativistic gamma, but Ekin/m
   TGraph* egr = data.GetGraphXY("SPOS","ENERGY");
   Double_t* ekin = egr->GetY();
   Double_t betagamma[nStep];
   const Double_t mel = 0.51099892;
   for (Int_t i = 0; i<nStep; i++) {
     betagamma[i] = sqrt(ekin[i]*ekin[i]+2.*ekin[i]*mel)/mel;
   }

   // set slice options:
   data.SetSlicePosition("mean");
   data.SetSlicePercent(90.);   // 10% slice by default

   //---------------------------------------------------------------------
   // emittances and RMS:
   //---------------------------------------------------------------------

   // emittance stored as #varepsilon is already normalized!
   // emittance in x
   TGraph* emittx = data.GetGraphXY("SPOS","#varepsilon","x",emittScale);
   emittx->SetLineWidth(lineWidth);
   emittx->SetLineStyle(kSolid);
   emittx->SetLineColor(kRed);

   // emittance in y
   TGraph* emitty = data.GetGraphXY("SPOS","#varepsilon","y",emittScale);
   emitty->SetLineWidth(lineWidth);
   emitty->SetLineStyle(kDashed);
   emitty->SetLineColor(kBlue);

   // slice emittance in x
   TGraph* semittx;
   if (doSlice) {
     //    semittx = data.GetGraphXY("SPOS","#varepsilon","x",emittScale);
     semittx = data.GetGraphXY("SPOS","emittance(x,px)","slice",emittScale);
   } else {
     semittx = data.GetGraphXY("SPOS","#varepsilon","x",emittScale);
   }
   semittx->SetLineWidth(lineWidth-2);
   semittx->SetLineStyle(kSolid);
   semittx->SetLineColor(kRed);

   // slice emittance in y
   TGraph* semitty;
   if (doSlice) {
     //    semitty = data.GetGraphXY("SPOS","#varepsilon","y",emittScale);
     semitty = data.GetGraphXY("SPOS","emittance(y,py)","slice",emittScale);
   } else {
     semitty = data.GetGraphXY("SPOS","#varepsilon","y",emittScale);
   }
   semitty->SetLineWidth(lineWidth-2);
   semitty->SetLineStyle(kDashed);
   semitty->SetLineColor(kBlue);

   // rms in x
   TGraph* rmsx = data.GetGraphXY("SPOS","RMSX","x",rmsScale);
   rmsx->SetLineWidth(lineWidth);
   rmsx->SetLineStyle(kSolid);
   rmsx->SetLineColor(kBlack);

   // rms in y
   TGraph* rmsy = data.GetGraphXY("SPOS","RMSX","y",rmsScale);
   rmsy->SetLineWidth(lineWidth);
   rmsy->SetLineStyle(kDashed);
   rmsy->SetLineColor(kGreen);

   //------------------------------------------------------------------------------

   TCanvas* finCanv = new TCanvas("finCanv","H5PartROOT",800,400);

   TString globtit = TString("SwissFEL injector test facility, phase 1");
   //  globtit += filename;
   TString myname = TString("OPAL");

   TPaveLabel* globTitle = new TPaveLabel(0.1,0.95,0.9,0.99,globtit);
   globTitle->SetBorderSize(0);
   globTitle->SetFillColor(0);
   globTitle->SetTextFont(62);
   //  globTitle->Draw();
   TDatime now;
   TPaveLabel* date = new TPaveLabel(0.7,0.01,0.9,0.05,now.AsString());
   date->SetBorderSize(0);
   date->SetFillColor(0);
   date->SetTextFont(62);
   date->SetTextSize(0.3);
   //  date->Draw();
   TPaveLabel* nameLabel = new TPaveLabel(0.05,0.01,0.25,0.05,myname);
   nameLabel->SetBorderSize(0);
   nameLabel->SetFillColor(0);
   nameLabel->SetTextFont(62);
   nameLabel->SetTextSize(0.3);
   //  nameLabel->Draw();
   TPad* graphs = new TPad("c1Graphs","Graphs",0.01,0.05,0.95,0.95);
   graphs->SetBottomMargin(4.);
   graphs->Draw();
   //  graphs->Divide(1,3);

   graphs->cd(1);gPad->SetRightMargin(0.1);gPad->SetBottomMargin(0.16);
   graphs->cd(2);gPad->SetRightMargin(0.1);gPad->SetBottomMargin(0.16);
   graphs->cd(3);gPad->SetRightMargin(0.1);gPad->SetBottomMargin(0.16);

   gStyle->SetTickLength(0.03,"x");
   gStyle->SetTickLength(0.01,"y");
   gStyle->SetPadTickX(0);
   gStyle->SetPadTickY(0);

   graphs->cd(1);
   gPad->SetGrid();

   TH1F *frame = gPad->DrawFrame(sMin,0.,sMax,EmittMax);
   frame->SetXTitle("s [m]");
   frame->SetTitleSize(titSize);
   frame->SetLabelSize(labSize);
   frame->GetYaxis()->SetTitleOffset(titOffset+0.2);
   frame->SetYTitle("emittance [mm mrad]");

   // second axis to the right:
   TGaxis *axis = new TGaxis(sMax,0.,sMax,EmittMax,0.,RmsMax,506,"SL+");
   axis->SetTitle("beam size RMS [mm]");
   axis->SetTitleSize(titSize);
   axis->SetTitleOffset(titOffset);
   axis->SetLabelSize(labSize);
   axis->SetLabelOffset(0.005);
   axis->SetTickSize(0.01);
   axis->Draw();

   emittx->Draw("Lsame");
   emitty->Draw("Lsame");
   semittx->Draw("Lsame");
   semitty->Draw("Lsame");
   rmsx->Draw("Lsame");
   rmsy->Draw("Lsame");

   // legends
   TLegend *legEmit;
   if (doSlice) {
     legEmit = new TLegend(0.7,0.2,0.88,0.5);
     legEmit->SetMargin(0.4);
   } else {
     legEmit = new TLegend(0.7,0.2,0.85,0.35);
     legEmit->SetMargin(0.5);
   }
   legEmit->SetTextFont(62);
   legEmit->SetTextSize(0.05);
   legEmit->SetFillColor(0);
   legEmit->SetShadowColor(0);
   legEmit->AddEntry(emittx, "#varepsilon_{x}","l");
   legEmit->AddEntry(emitty, "#varepsilon_{y}","l");
   if (doSlice) {
     legEmit->AddEntry(semittx, "#varepsilon_{x,core (90%)}","l");
     legEmit->AddEntry(semitty, "#varepsilon_{y,core (90%)}","l");
   }
   legEmit->Draw();

   TLegend *legRMS;
   legRMS = new TLegend(0.75,0.55,0.88,0.7);
   legRMS->SetMargin(0.5);
   legRMS->SetTextFont(62);
   legRMS->SetTextSize(0.05);
   legRMS->SetFillColor(0);
   legRMS->SetShadowColor(0);
   legRMS->AddEntry(rmsx, "#sigma_{x}","l");
   legRMS->AddEntry(rmsy, "#sigma_{y}","l");
   legRMS->Draw();

   //---------------------------------------------------------------------
   // RMS and emittance measurements
   //---------------------------------------------------------------------

   // example for how to add a few measurements as TGraphErrors:
   // (read in a file, or add measurment by hand)

   TGraphErrors *measSigmaX = new TGraphErrors(10); // <-- number of points
   measSigmaX->SetMarkerColor(kBlack);
   measSigmaX->SetLineColor(kBlack);
   measSigmaX->SetMarkerSize(0.8);
   measSigmaX->SetMarkerStyle(20);

   TGraphErrors *measSigmaY = new TGraphErrors(10);
   measSigmaY->SetMarkerColor(kGreen+3);
   measSigmaY->SetLineColor(kGreen+3);
   measSigmaY->SetMarkerSize(0.8);
   measSigmaY->SetMarkerStyle(20);

   // read in data from file:
   TString dataFile = measPath+rmsFilename;
   FILE *fileRMS = fopen(dataFile.Data(),"r");

   Int_t ncols;
   Int_t nlines = 0;
   // variables for screen number, sigma(s) from gauss(g) and rms(r) with error(e)
   Int_t n;
   Float_t z, sgx, sgy, segx, segy, srx, sry, serx, sery;

   // skip the first line (header)
   char header[256];
   fgets(header,256,fileRMS);

   Double_t offset = 0.01;

   while (1) {
     ncols = fscanf(fileRMS,
		   "%d %e %e %e %e %e %e %e %e %e",
		   &n, &z, &sgx, &sgy, &segx, &segy, &srx, &sry, &serx, &sery);
     if (ncols < 0) break;
     //    if (nlines < 5) printf("n = %d sgx=%8f, sgy=%8f\n",n,sgx,sgy);

     measSigmaX->SetPoint(nlines,z-offset,srx*rmsScale/1000.);       // values are in mm already
     measSigmaX->SetPointError(nlines,0.,serx*rmsScale/1000.);

     measSigmaY->SetPoint(nlines,z+offset,sry*rmsScale/1000.);
     measSigmaY->SetPointError(nlines,0.,sery*rmsScale/1000.);

     nlines++;
   }

   measSigmaX->Draw("Psame");
   measSigmaY->Draw("Psame");

   TGraphErrors *measEmittX = new TGraphErrors(1);
   measEmittX->SetMarkerColor(kRed+2);
   measEmittX->SetLineColor(kRed+2);
   measEmittX->SetMarkerSize(1.0);
   measEmittX->SetMarkerStyle(20);

   TGraphErrors *measEmittY = new TGraphErrors(1);
   measEmittY->SetMarkerColor(kBlue+2);
   measEmittY->SetLineColor(kBlue+2);
   measEmittY->SetMarkerSize(1.0);
   measEmittY->SetMarkerStyle(20);

   // just one mesaurement of emittance => do it by hand
   //measEmittX->SetPoint(0,2.662-offset,1.2806);
   //measEmittX->SetPointError(0,0.,0.0551);

   //measEmittY->SetPoint(0,2.662+offset,1.2118);
   //measEmittY->SetPointError(0,0.,0.0551);

   // correted value obtained from: Natalia
   measEmittX->SetPoint(0,2.662-offset,0.77);
   measEmittX->SetPointError(0,0.,0.11);

   //   measEmittY->SetPoint(0,2.662+offset,2.3);
   // measEmittY->SetPointError(0,0.,0.3);

   measEmittX->Draw("Psame");
   measEmittY->Draw("Psame");

   TText *t = new TText();
   t->SetTextFont(62);
   t->SetTextColor(1);
   t->SetTextSize(0.04);
   t->SetTextAlign(12);
   t->DrawTextNDC(0.27,0.85,"23 June 2010");
   t->DrawTextNDC(0.27,0.81,"Q = 100 pC");
   t->DrawTextNDC(0.27,0.77,"p = 7.0 MeV/c");
   t->DrawTextNDC(0.27,0.73,"laser aperture: 2.6 mm");

   finCanv->Update();

   finCanv->Print("plotFinPhase1LINAC-100pc.pdf");
   finCanv->Print("plotFinPhase1LINAC-100pc.png");
   finCanv->Print("plotFinPhase1LINAC-100pc.eps");
   finCanv->Print("plotFinPhase1LINAC-100pc.cxx");

}
