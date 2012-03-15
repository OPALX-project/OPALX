{
   #include "Riostream.h"

   Int_t mStyleArr[10] = {2,24,4,5,6,24,25,26,27,28};

   TString* fn[10];
   TNtuple* nt[10];
   TString baseDir;

   gSystem->Load("~schietinger/bin/H5root.so");          

   TH5Style::SetStyle();            

   gStyle->SetStatX(0.9);
   gStyle->SetStatY(0.9);

   c1 = new TCanvas("c1","PSI-FEL/LEG Gun Analysis" ,200,10,700,500);
   gPad->SetRightMargin(0.1);


/* 
   Make up the filenames we need	
   fn[0] = new TString("gun0.5MV-2.cfg-gun-SL380-BINS32-NP3125-DT1.0e-13");
*/
Int_t dataSets = 5;  //  fn[0] = new TString("gun0.5MV-2.cfg-gun-SL380-BINS32-NP3125-DT1.0e-13");
   TString dir("0.5MV");
   TString fext("pdf");
   fn[0] = new TString("gun0.5MV-1.cfg-gun-SL380-BINS32-NP31250-DT1.0e-13");
   fn[1] = new TString("gun0.5MV-2.cfg-gun-SL380-BINS32-NP3125-DT1.0e-13");
   fn[2] = new TString("gun0.5MV-3.cfg-gun-SL380-BINS32-NP3125-DT1.0e-14");
   fn[3] = new TString("gun0.5MV-4.cfg-gun-SL380-BINS32-NP3125-DT1.0e-14");
   fn[4] = new TString("gun0.5MV-5.cfg-gun-SL380-BINS32-NP31250-DT1.0e-14");



/*
   TString dir("1MV");
   fn[0] = new TString("gun1MV-1.cfg-gun-SL380-BINS32-NP31250-DT1.0e-13");
   fn[1] = new TString("gun1MV-2.cfg-gun-SL380-BINS32-NP31250-DT1.0e-14");
   fn[2] = new TString("gun1MV-3.cfg-gun-SL380-BINS32-NP31250-DT1.0e-13");
*/

   for (Int_t i=0; i<dataSets; i++) {
    TGraph *h5Tgr = processH5File(c1, dir, *fn[i], fext);
    //    processFortFile(c1, dir, *fn[i]);
   }

//   fn[3] = new TString("mafia_100k.data");
//   processMafiaFile(c1, TString("."), *fn[3]);




}
