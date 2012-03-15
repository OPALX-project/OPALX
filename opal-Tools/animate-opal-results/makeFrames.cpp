#include <cmath>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

#include "TString.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TRandom.h"

#include "Riostream.h"

#include "TH5Style.h"
#include "TH5Dataset.h"

using namespace std;
using namespace TH5Util;

#define epsilon 0.001
#define NDuplicates 4
#define TIMETOWAIT 2

typedef vector<TH2D*>::iterator HIterator;
typedef vector<TPad*>::iterator PIterator;

TRandom* RGenerator;
TString fileFormat;

struct PLOT {
    PLOT(Double_t x1, Double_t y1, Double_t x2, Double_t y2, TString vx, TString vy, TString op="") {
        PosX1 = x1;
        PosY1 = y1;
        PosX2 = x2;
        PosY2 = y2;
        VarX = vx;
        VarY = vy;
        Options = op;
    }
    
    Double_t PosX1;
    Double_t PosY1;
    Double_t PosX2;
    Double_t PosY2;
    Double_t MinX;
    Double_t MaxX;
    Double_t MinY;
    Double_t MaxY;
    Double_t BinMax;
    TString VarX;
    TString VarY;
    TString Options;
};

struct GRAPH {
    GRAPH(Double_t x1, Double_t y1, Double_t x2, Double_t y2, TString vx, TString vy, TString cp, TString op="") {
        PosX1 = x1;
        PosY1 = y1;
        PosX2 = x2;
        PosY2 = y2;
        VarX = vx;
        VarY = vy;
        Component = cp;
        Options = op;
    }
    
    Double_t PosX1;
    Double_t PosY1;
    Double_t PosX2;
    Double_t PosY2;
    Double_t MinX;
    Double_t MaxX;
    Double_t MinY;
    Double_t MaxY;
    TString VarX;
    TString VarY;
    TString Component;
    TString Options;
};

typedef vector<PLOT*>::iterator PLIterator;
typedef vector<GRAPH*>::iterator GIterator;

void make_movie(TString Name, Int_t first, Int_t last) {
    char printFilename[100], sysCommand[1000], tmpCommand[1000];
    char copyFrom[100], copyTo[100];
    TString copyFromTmpl = TString("Frames/") + Name + fileFormat;
    TString copyToTmpl = TString("Frames/dir/") + Name + fileFormat;
    TString videoFormat(".mpg");
    TString videoFName = TString("Frames/dir/") + Name + videoFormat;
    TString videoOptions("-b 10000k");

    mkdir("Frames/dir", 484);
    
    
    for (Int_t i = last; i >= first; -- i) {
        sprintf(copyFrom, copyFromTmpl.Data(), i);
        for (Int_t j = NDuplicates; j > 0; -- j) {
            sprintf(copyTo, copyToTmpl.Data(), NDuplicates*i + j - 1);
            sprintf(sysCommand, "cp %s %s", copyFrom, copyTo);
            system(sysCommand);
        }
    }
    sprintf(sysCommand, "ffmpeg -y -an -i %s %s %s", copyToTmpl.Data(), videoOptions.Data(), videoFName.Data());
    system(sysCommand);
    
    sprintf(sysCommand, "mv %s %s%s", videoFName.Data(), Name.Data(), videoFormat.Data());
    system(sysCommand);

}

TString GetUniqueName(Int_t Length) {
    TString Name("");
    TString ABC("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789");
    for (Int_t i = 0; i < Length; ++ i) {
        Int_t idx = static_cast<Int_t>(ceil(RGenerator->Uniform(ABC.Sizeof())));
        Name += ABC(idx,1);
    }
    return Name;
}

TString getFileNameBase(const char* fName, Int_t &first, Int_t &last) {
    TString Filename(fName);
    
    TRegexp reRange("\\[[1-9][0-9]*:[1-9][0-9]*\\]$");
    TRegexp reExtension("\\.h5$");
    TString strRange(Filename(reRange));

    if (strRange.Sizeof() >= 6) {
        Ssiz_t idxColumn = strRange.First(":");
        Ssiz_t idxBracket = strRange.Sizeof() - 2;
        stringstream strFirst;
        stringstream strLast;

        strFirst << strRange(1, idxColumn-1);
        strFirst >> first;

        strLast << strRange(idxColumn+1, idxBracket - idxColumn - 1);
        strLast >> last;

        Filename.Remove(Filename.Index(reRange));
        
    }
    Int_t idxExtension = Filename.Index(reExtension);
    if (idxExtension > 0) {
        Filename.Remove(idxExtension);
    }

    return Filename;
}

void ProcessStep(TH5Dataset & data, 
                 Int_t step, 
                 Int_t printStep, 
                 vector<PLOT*> histoplots, 
                 vector<GRAPH*> graphplots, 
                 TString fNBase) {
    cout << "processing step " << step << endl;
    Int_t LineColors[6] = {kRed, kBlue, kGreen, kMagenta, kOrange};
    vector<TH2D*> histos;
    vector<TPad*> pads;
    vector<vector<TGraph*> > graphs;
    vector<TLegend*> legends;
    char printFilename[100];

    fNBase += fileFormat;
    sprintf(printFilename, fNBase, printStep);

    TCanvas* dummy = new TCanvas("dummy","dummy",10,10);

    Double_t minX, maxX, minY, maxY, dX, dY, meanX, meanY;
    Int_t fNContour = 20;

    /*
    head = maxX = data.GetMaximum(vxname1,step);
    tail = minX = data.GetMinimum(vxname1,step);
    if (movingWindow){
      if(maxX > 0.8 * Winmax){
        Winmax = maxX + 0.1;
        Winmin = Winmax - 0.5;
      }
      if(Winmax > Winmax_orig){
        Winmax = Winmax_orig;
        Winmin = Winmax_orig - 0.5;
      }
    }
    */

    for (PLIterator it = histoplots.begin(); it != histoplots.end(); ++ it) {
        
        minX = data.GetMinimum((*it)->VarX, step);
        maxX = data.GetMaximum((*it)->VarX, step);
        minY = data.GetMinimum((*it)->VarY, step);
        maxY = data.GetMaximum((*it)->VarY, step);

        meanX = 0.5 * ((*it)->MaxX + (*it)->MinX);
        meanY = 0.5 * ((*it)->MaxY + (*it)->MinY);
        dX = (*it)->MaxX - (*it)->MinX;
        dY = (*it)->MaxY - (*it)->MinY;

        if (maxX - meanX < 0.4 * dX || maxX - meanX > 0.5 * dX)
            (*it)->MaxX = 1.05 * maxX - 0.05 * minX;
        if (meanX - minX < 0.4 * dX || meanX - minX > 0.5 * dX) 
            (*it)->MinX = 1.05 * minX - 0.05 * maxX;
	if ((*it)->MaxX - (*it)->MinX < epsilon) {
          (*it)->MaxX = (*it)->MaxX + epsilon;
          (*it)->MinX = (*it)->MinX - epsilon;
	}
        
        if (maxY - meanY < 0.4 * dY || maxY - meanY > 0.5 * dY) 
            (*it)->MaxY = 1.05 * maxY - 0.05 * minY;
        if (meanY - minY < 0.4 * dY || meanY - minY > 0.5 * dY) 
            (*it)->MinY = 1.05 * minY - 0.05 * maxY;
	if ((*it)->MaxY - (*it)->MinY < epsilon) {
          (*it)->MaxY = (*it)->MaxY + epsilon;
          (*it)->MinY = (*it)->MinY - epsilon;
	}
        
        histos.push_back(data.Histo2d((*it)->VarX.Data(), (*it)->VarY.Data(), step, 
                                      (*it)->MinX, (*it)->MaxX, (*it)->MinY, (*it)->MaxY));
        histos.back()->SetStats(kFALSE);
        histos.back()->SetLabelSize(0.04,"X");
        histos.back()->SetLabelSize(0.04,"Y");
        histos.back()->SetLabelSize(0.04,"Z");
        //histo->SetTitleOffset(1,3, "Y");
        histos.back()->GetXaxis()->SetTitle((*it)->VarX + TH5Util::AddBrackets(data.GetUnit((*it)->VarX)));
        histos.back()->GetYaxis()->SetTitle((*it)->VarY + TH5Util::AddBrackets(data.GetUnit((*it)->VarY)));
        Double_t BinLocalMax = histos.back()->GetMaximum();
        if (BinLocalMax < 0.8 * (*it)->BinMax || BinLocalMax > (*it)->BinMax) 
            (*it)->BinMax = ceil(BinLocalMax * 0.1) * 12.0;
        histos.back()->SetAxisRange(0.0, (*it)->BinMax, "Z");
        histos.back()->SetContour(fNContour);
    }

    for (GIterator it = graphplots.begin(); it != graphplots.end(); ++ it) {
        vector<TGraph*> tmp;
        Double_t max = 1.0;
        TString yvars((*it)->VarY);
        Int_t NGraphs = 0;
        
        legends.push_back(new TLegend(0.9, 0.8, 1.0, 1.0));
        legends.back()->SetTextFont(62);
        legends.back()->SetTextSize(0.08);
        legends.back()->SetFillColor(kWhite);
        legends.back()->SetBorderSize(0);
        legends.back()->SetEntrySeparation(0.5);

        if (yvars.First("+") > 0) {
            Int_t idx = yvars.First("+");
            TString var = yvars(0, idx);

            TGraph* tmp = data.GetGraphXY((*it)->VarX, var, (*it)->Component);
            TH1F* tmphst = tmp->GetHistogram();
            max = fabs(tmphst->GetMaximum());
            max = max > fabs(tmphst->GetMinimum())? max : fabs(tmphst->GetMinimum());
        }

        while (yvars.Last('+') > 0) {
            Int_t idx = yvars.Last('+');
            TString var = yvars(idx+1, yvars.Sizeof() - idx - 2);
            yvars.Remove(idx);

            TGraph* tmpgraph = data.GetGraphXY((*it)->VarX, var, (*it)->Component);
            TH1F* tmphst = tmpgraph->GetHistogram();
            Double_t ratio = fabs(tmphst->GetMaximum());
            ratio = ratio > fabs(tmphst->GetMinimum())? ratio : fabs(tmphst->GetMinimum());
            if (ratio > epsilon) {
                ratio = max / ratio;
            } else {
                ratio = 1.0;
            }

            TString legendEntry = var;
            if ((*it)->Component == TString("X") || 
                (*it)->Component == TString("Y") || 
                (*it)->Component == TString("Z")) {
                TString comp = (*it)->Component;
                comp.ToLower();
                legendEntry += TString("{") + comp +  TString("}");
            }

            TGraph *graph = data.GetGraphXY((*it)->VarX, TString(""), 1.0, var, (*it)->Component, ratio);

            graph->SetMarkerStyle(1);
            graph->SetMarkerColor(LineColors[NGraphs % 5]);
            graph->SetMarkerSize(0.5);
            graph->SetLineColor(LineColors[NGraphs % 5]);
            graph->SetLineWidth(2);
            tmp.push_back(graph);


            legends.back()->AddEntry(graph, legendEntry, "l");
            ++ NGraphs;
        }

        TString legendEntry = yvars;
        if ((*it)->Component == TString("X") || 
            (*it)->Component == TString("Y") || 
            (*it)->Component == TString("Z")) {
            TString comp = (*it)->Component;
            comp.ToLower();
            legendEntry += TString("{") + comp +  TString("}");
        }

        TGraph *graph = data.GetGraphXY((*it)->VarX, yvars, (*it)->Component);
        graph->SetMarkerStyle(1);
        graph->SetMarkerColor(LineColors[NGraphs % 5]);
        graph->SetMarkerSize(0.5);
        graph->SetLineColor(LineColors[NGraphs % 5]);
        graph->SetLineWidth(2);
        tmp.push_back(graph);

        legends.back()->AddEntry(graph, legendEntry, "l");

        graphs.push_back(tmp);
        
    }


    delete dummy;

    
    TCanvas* myCanv = new TCanvas("myCanv","H5PartROOT",1026,660);

    HIterator hit = histos.begin();
    for (PLIterator it = histoplots.begin(); it != histoplots.end(); ++ it) {
        TString padName = GetUniqueName(10);

        pads.push_back(new TPad(padName.Data(), padName.Data(), (*it)->PosX1, (*it)->PosY1, (*it)->PosX2, (*it)->PosY2));
        pads.back()->SetTopMargin(0.1);
        pads.back()->SetRightMargin(0.1);
        pads.back()->SetBottomMargin(0.1);
        pads.back()->SetLeftMargin(0.1);
        pads.back()->Draw();
        pads.back()->cd();
        (*hit)->Draw("CONT");
        myCanv->cd();
        ++ hit;
    }
    
    vector<vector<TGraph*> >::iterator gvit = graphs.begin();
    vector<TLegend*>::iterator lit = legends.begin();
    vector<TLine*> tlv;
    for (GIterator it = graphplots.begin(); it != graphplots.end(); ++ it) {
        TString padName = GetUniqueName(10);
        Double_t MinY = 1.05 * (*it)->MinY;
        Double_t MaxY = 1.05 * (*it)->MaxY;

        pads.push_back(new TPad(padName.Data(), padName.Data(), (*it)->PosX1, (*it)->PosY1, (*it)->PosX2, (*it)->PosY2));
        pads.back()->SetTopMargin(0.1);
        pads.back()->SetRightMargin(0.1);
        pads.back()->SetBottomMargin(0.1);
        pads.back()->SetLeftMargin(0.1);
        pads.back()->Draw();
        pads.back()->cd();
        
        TH1F *frame = pads.back()->DrawFrame((*it)->MinX, MinY, (*it)->MaxX, MaxY);
        frame->SetXTitle("");
        frame->SetYTitle("");
        frame->SetLabelSize(0.04, "X");
        frame->SetTickLength(0.0, "Y");
        frame->SetLabelSize(0.0, "Y");

        for (vector<TGraph*>::iterator git = (*gvit).begin(); git != (*gvit).end(); ++ git) {
            (*git)->Draw("LPsame");
        }

        if ((*it)->VarX == TString("spos-ref")) {
            Double_t tail = data.GetMinimum(TString("z"), step);
            Double_t head = data.GetMaximum(TString("z"), step);

            TLine *l1 = new TLine(tail, MinY, tail, MaxY);
            TLine *l2 = new TLine(head, MinY, head, MaxY);
        
            l1->Draw();
            l2->Draw();
            
            tlv.push_back(l1);
            tlv.push_back(l2);
        }
	if ((*it)->VarX == TString("SPOS")) {
	    TString var("SPOS");
	    Double_t spos = data.GetScalarStepVar(&var, step);

            TLine *l = new TLine(spos, MinY, spos, MaxY);
        
            l->Draw();
            
            tlv.push_back(l);
	}
        
        (*lit)->Draw();

        myCanv->cd();
        ++ gvit;        
        ++ lit;
    }

    myCanv->Update();
    myCanv->Print(TString(printFilename));

    for (vector<TLine*>::iterator it = tlv.begin(); it != tlv.end(); ++ it) {
        delete (*it);
    }
    delete myCanv;
}

int main(int argc, char **argv) {

  RGenerator = new TRandom();
  fileFormat = TString("%07d.png");
  mkdir("Frames",484);
  system("rm Frames/*.* Frames/dir/*");

  char * name;
  int every_std = 1;
  if (argc > 1) {
    name = argv[1];
  } else {
    cerr << "Usage: AnimateSimulation input.h5" << endl;
    return 0;
  }

  Int_t fStep = 0, lStep = 10000000;

  TH5Style::SetStyle();

  TString fName(getFileNameBase(name, fStep, lStep));
  // load the data file:
  TH5Dataset data(fName + TString(".h5"), fStep, lStep);

  fStep = data.GetFirstStep();
  lStep = data.GetLastStep();

  Int_t step = fStep, every = every_std;
  Double_t minX, maxX, minY, maxY, dX, dY;
  Double_t tail, head;
  vector<PLOT*> histos;
  vector<GRAPH*> graphs;

  histos.push_back(new PLOT(0.0, 0.285, 0.683, 1.0, "x", "y"));
  histos.push_back(new PLOT(0.683, 0.5, 1.0, 1.0, "z", "pz"));
  histos.push_back(new PLOT(0.683, 0.0, 1.0, 0.5, "x", "px"));

  graphs.push_back(new GRAPH(0.0, 0.0, 0.683, 0.2, "SPOS", "B-ref", "Z"));

  for (PLIterator it = histos.begin(); it != histos.end(); ++ it) {
      Double_t DataMax = data.GetMaximum((*it)->VarX, fStep);
      Double_t DataMin = data.GetMinimum((*it)->VarX, fStep);
      Double_t DataRange = DataMax - DataMin;
      if (DataRange > 2 * epsilon) {
          (*it)->MaxY = DataMin + 1.05 * DataRange;
          (*it)->MinY = DataMax - 1.05 * DataRange;
      } else {
          (*it)->MaxY = DataMin + 2 * epsilon;
          (*it)->MinY = DataMax - 2 * epsilon;
      }

      DataMax = data.GetMaximum((*it)->VarY, fStep);
      DataMin = data.GetMinimum((*it)->VarY, fStep);
      DataRange = DataMax - DataMin;
      if (DataRange > 2 * epsilon) {
          (*it)->MaxY = DataMin + 1.05 * DataRange;
          (*it)->MinY = DataMax - 1.05 * DataRange;
      } else {
          (*it)->MaxY = DataMin + 2 * epsilon;
          (*it)->MinY = DataMin - 2 * epsilon;
      }
  }

  for (GIterator it = graphs.begin(); it != graphs.end(); ++ it) {
      TString yvars((*it)->VarY);
      Double_t max = 0.0, min = 0.0;
      while (yvars.Index("+") > 0) {
          Int_t idx = yvars.Index("+");
          TString var = yvars(idx+1, yvars.Sizeof() - idx - 2);
          yvars.Remove(idx);
          TGraph* graph = data.GetGraphXY((*it)->VarX, var, (*it)->Component);
          TH1F* hist = graph->GetHistogram();
          max = max > hist->GetMaximum()? max: hist->GetMaximum();
          min = min < hist->GetMinimum()? min: hist->GetMinimum();
      }

      TGraph* graph = data.GetGraphXY((*it)->VarX, yvars, (*it)->Component);
      TH1F* hist = graph->GetHistogram();
      max = max > hist->GetMaximum()? max: hist->GetMaximum();
      min = min < hist->GetMinimum()? min: hist->GetMinimum();
      
      (*it)->MaxY = max > fabs(min)? max : fabs(min);
      (*it)->MinY = - (*it)->MaxY;

      (*it)->MinX = 9999999.;
      (*it)->MaxX = -9999999.;
      Double_t *X = graph->GetX();
      for (Int_t i = 0; i < graph->GetN(); ++ i) {
          if (X[i] > (*it)->MaxX) (*it)->MaxX = X[i];
          if (X[i] < (*it)->MinX) (*it)->MinX = X[i];
      }
  }

  for (step = fStep; step <= lStep; step += every){
      Int_t printStep = (step - fStep) / every;
      ProcessStep(data, step, printStep, histos, graphs, TString("Frames/")+fName);
  }

  make_movie(fName, 0, (lStep-fStep) / every);

  while (true) {
      if (data.DoReload()) {
	Int_t newLastStep = data.GetLastStep();
        
          if (newLastStep > lStep) {
              for (GIterator it = graphs.begin(); it != graphs.end(); ++ it) {
                  TString yvars((*it)->VarY);
                  Double_t max = 0.0, min = 0.0;
                  while (yvars.Index("+") > 0) {
                      Int_t idx = yvars.Index("+");
                      TString var = yvars(idx+1, yvars.Sizeof() - idx - 2);
                      yvars.Remove(idx);
                      TGraph* graph = data.GetGraphXY((*it)->VarX, var, (*it)->Component);
                      TH1F* hist = graph->GetHistogram();
                      max = max > hist->GetMaximum()? max: hist->GetMaximum();
                      min = min < hist->GetMinimum()? min: hist->GetMinimum();
                  }
                  TGraph* graph = data.GetGraphXY((*it)->VarX, yvars, (*it)->Component);
                  TH1F* hist = graph->GetHistogram();
                  max = max > hist->GetMaximum()? max: hist->GetMaximum();
                  min = min < hist->GetMinimum()? min: hist->GetMinimum();
      
                  (*it)->MaxY = max > fabs(min)? max : fabs(min);
                  (*it)->MinY = - (*it)->MaxY;

                  Double_t *X = graph->GetX();
                  for (Int_t i = lStep - fStep + 1; i < graph->GetN(); ++ i) {
                      if (X[i] > (*it)->MaxX) (*it)->MaxX = X[i];
                      if (X[i] < (*it)->MinX) (*it)->MinX = X[i];
                  }
              }

              for (step = lStep + 1; step <= newLastStep; step += every) {
                  Int_t printStep = (step - fStep) / every;
                  ProcessStep(data, step, printStep, histos, graphs, TString("Frames/") + fName);
              }

              make_movie(fName, (lStep-fStep) / every + 1, (newLastStep - fStep) / every);
              lStep = newLastStep;
          }

      }
        
      TH5Util::Wait(TIMETOWAIT);
  }

  delete RGenerator;
}

