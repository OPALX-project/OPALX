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
