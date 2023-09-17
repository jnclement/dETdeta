#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TMath.h"
#include <utility>  // for std::pair
#include <cstdio>
#include <iostream>
#include "TGraph.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include "stdlib.h"
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>
#include <TFile.h>
#include <algorithm>
#include "dlUtility.h"
#include "mbd_info.h"
#include "/home/jocl/Documents/main/physics/projects/sphenix_macros/macros/macros/sPHENIXStyle/sPhenixStyle.h"
#include "/home/jocl/Documents/main/physics/projects/sphenix_macros/macros/macros/sPHENIXStyle/sPhenixStyle.C"

void plothists(TCanvas* ca, TH1* histarr, int nhist, string* text, int ntext, int logy)
{
  ca->cd();
}

void plotsimdat(TCanvas* ca, TH1* dathist, TH1* simhist, int logy, string cal, float sc, float sub, int run, float mine, float zcut, string xlabel, int percent0, int percent1, string name, string dir)
{
  const int par = 4;
  float parval[par];
  float mult[par];
  mult[0] = 1;
  mult[1] = 1000;
  mult[2] = 1000;
  mult[3] = 1;
  parval[0] = sc;
  parval[1] = sub;
  parval[2] = mine;
  parval[3] = zcut;
  string params[par];
  stringstream streams[par];
  int precision[par];
  precision[0] = 2;
  precision[1] = 0;
  precision[2] = 0;
  precision[3] = 0;
  for(int i=0; i<par; ++i)
    {
      streams[i] << std::fixed << std::setprecision(precision[i]) << parval[i]*mult[i];
      params[i] = streams[i].str();
    }
  const int ntext = 4;
  string texts[ntext];
  texts[0] = "Red data, blue sim, sim scaled by " + params[0];
  texts[1] = "Area normed to 1, |z|<" + params[3] +", min tower E = " +params[2];
  texts[2] = params[1] + " subtracted from each tower";
  texts[3] = "Run " + to_string(run) + " " + to_string(5*(20-percent1)) +"-"+ to_string(5*(20-percent2))+"% centrality"
  ca->cd();
  if(logy) gPad->SetLogy();
  else gPad->SetLogy(0);
  gPad->SetTicks(1);
  dathist->Scale(1./dathist->Integral());
  simhist->Scale(1./simhist->Integral());
  dathist->SetLineColor(kRed);
  simhist->SetLineColor(kBlue);
  float maxval = max(dathist->GetMaximum(),simhist->GetMaximum());
  float minval = min(dathist->GetBinContent(dathist->FindLastBinAbove(0,1)),simhist->GetBinContent(simhist->FindLastBinAbove(0,1)));
  dathist->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  dathist->GetXaxis()->SetTitle((xlabel).c_str());
  dathist->GetYaxis()->SetTitle("Counts");
  dathist->GetYaxis()->SetLabelSize(0.025);
  dathist->GetXaxis()->SetLabelSize(0.025);
  dathist->Draw();
  simhist->Draw("SAME");
  sphenixtext();
  multitext(texts, ntext);
  ca->SaveAs((dir+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(5*(20-percent1)) +"-"+ to_string(5*(20-percent2))+".pdf").c_str());

  dathist->Divide(simhist);
  dathist->GetYaxis()->SetTitle("Data/Sim");
  dathist->GetYaxis()->SetRangeUser(0.01,10);
  dathist->Draw();
  sphenixtext();
  drawText(("Run " + to_string(run) + " " + to_string(5*(20-percent1)) +"-"+ to_string(5*(20-percent2))+"% centrality").c_str(),0.85,0.85,1,kBlack,0.025);
  ca->SaveAs((dir+"ratio"+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(5*(20-percent1)) +"-"+ to_string(5*(20-percent2))+".pdf").c_str());
}

int plot(string histfilename = "savedhists_subtr_0_minE_0_scale_1.30_zcut_30.root", string treename = "ttree", string datdir = "/home/jocl/datatemp/", string plotdir = "/home/jocl/datatemp/plots/")
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  SetsPhenixStyle();
  const int centbins = 20;
  TFile* histfile = TFile::Open(datdir+histfilename);
  TTree* tree = histfile->Get<TTree>("treename");
}
