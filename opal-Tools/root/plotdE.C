//
//   Execute in root: .x plotdE.C
//   Change gauss-2.h5 to your desired input filename
//
{
   TString simPath("./");
   TString simFilename("Gun_2mm_225-optimized.h5");

   // load the shared library:
   gSystem->Load("/home5/schietinger/extlib/H5Root/lib/libh5root.so");
   using namespace TH5Util;

   // let's do it in style:
   TH5Style::SetStyle();
   gStyle->SetPadTickY(0);
   Width_t lineWidth(3);

   // load the data file:
   TH5Dataset data(simPath+simFilename);



   Int_t nsteps = data.GetNSteps();

   Double_t p[nsteps];
   Double_t s[nsteps];
   
   for(Int_t i = 0; i<nsteps; i++) {
       Double_t pmax = data.GetMaximum("pz",i);
       Double_t pmin = data.GetMinimum("pz",i);
       Double_t spos = data.GetSpos(i);

       TH1D *h = data.Histo1d("pz",i,pmin,pmax,256);
       Double_t prms = h.GetRMS();
       Double_t pmea = h.GetMean();

       cout << "Step= " << i << " prms= " << prms << " pmean= " << pmea << " SPOS= " << spos << endl;
       p[i] = prms;
       s[i] = spos;
   }

   TCanvas* finCanv = new TCanvas("finCanv","OPAL-Astra",800,400);

   gr = new TGraph(nsteps,s,p);
   gr->Draw("ACP");
   finCanv->Update();
   finCanv->Print("opal-astra-de.pdf");
   finCanv->Print("opal-astra-de.png");
}
