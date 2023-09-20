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
#include <iomanip>
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

void plotsimdat(TCanvas* ca, TH1* dathist, TH1* simhist, int logy, string cal, float sc, float sub, int run, float mine, float zcut, string xlabel, int percent0, int percent1, string name, string dir, string subdir, int centbins)
{
  int centrange = 90/centbins;
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
  const int ntext = 5;
  string texts[ntext];
  texts[0] = "Sim scaled by " + params[0] + ", Area normed to 1";
  texts[2] = "|z|<" + params[3] +" cm, min tower E = " +params[2] + " MeV";
  texts[1] = params[1] + " MeV subtracted from each tower";
  texts[3] = "Run " + to_string(run) + " " + to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  texts[4] = "Red data, blue sim";

  const int ntext2 = 3;
  string texts2[ntext2];
  texts2[0] = "Sim scaled by " + params[0] + ", |z|<" + params[3] + " cm, min tower E " + params[2] + " MeV";
  texts2[1] = params[1] + " MeV subtracted from each tower";
  texts2[2] = "Run " + to_string(run) + " " + to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  ca->cd();
  if(logy) gPad->SetLogy();
  else gPad->SetLogy(0);
  gPad->SetTicks(1);
  dathist->Scale(1./dathist->Integral());
  if(simhist) simhist->Scale(1./simhist->Integral());
  dathist->SetLineColor(kRed);
  dathist->SetMarkerColor(kRed);
  if(simhist)
    {
      simhist->SetLineColor(kBlue);
      simhist->SetMarkerColor(kBlue);
    }
  float maxval;
  if(simhist) maxval = max(dathist->GetMaximum(),simhist->GetMaximum());
  else maxval = dathist->GetMaximum();
  float minval;
  if(simhist) minval = min(dathist->GetBinContent(dathist->FindLastBinAbove(0,1)),simhist->GetBinContent(simhist->FindLastBinAbove(0,1)));
  else minval = dathist->GetBinContent(dathist->FindLastBinAbove(0,1));
  dathist->GetYaxis()->SetRangeUser(minval/2.,maxval*10.);
  dathist->GetXaxis()->SetTitle((xlabel).c_str());
  dathist->GetYaxis()->SetTitle("Counts");
  dathist->GetYaxis()->SetLabelSize(0.025);
  dathist->GetXaxis()->SetLabelSize(0.025);
  dathist->Draw();
  if(simhist) simhist->Draw("SAME");
  sphenixtext();
  multitext(texts, ntext);
  ca->SaveAs((dir+"pdf/"+subdir+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
  ca->SaveAs((dir+"png/"+subdir+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
  if(simhist)
    {
      dathist->Divide(simhist);
      dathist->GetYaxis()->SetTitle("Data/Sim");
      maxval = dathist->GetMaximum();
      minval = dathist->GetMinimum();
      dathist->GetYaxis()->SetRangeUser(0.01,100);
      dathist->Draw();
      sphenixtext();
      multitext(texts2, ntext2);
      ca->SaveAs((dir+"pdf/"+"ratio_"+subdir+"ratio_"+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
      ca->SaveAs((dir+"png/"+"ratio_"+subdir+"ratio_"+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
    }
}

int called_plot(string histfilename = "savedhists_subtr_0_minE_0_scale_1.30_zcut_30_run_21615.root", string treename = "ttree", string datdir = "/home/jocl/datatemp/", string plotdir = "/home/jocl/datatemp/plots/")
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  SetsPhenixStyle();
  const int centbins = 9;
  TFile* histfile = TFile::Open((datdir+histfilename).c_str());
  TTree* tree = histfile->Get<TTree>(treename.c_str());
  TH1D* centtow[2][3][centbins];
  TH1D* centet[2][3][centbins];
  TH1D* ET[2][3];
  TH1D* TW[2][3];
  TH1D* sumev[2];
  TH1D* sumtw[2];
  TH1D* mbh[2];
  TH1D* zhist;
  float sub;
  float scale[2];
  int frac[2];
  float mine;
  float zcut;
  int run;
  mbh[0] = (TH1D*)histfile->Get("smbh");
  mbh[1] = (TH1D*)histfile->Get("dmbh");
  zhist = (TH1D*)histfile->Get("zhist");
  for(int i=0; i<2; ++i)
    {
      sumev[i] = (TH1D*)histfile->Get(("sumev" + to_string(i)).c_str());
      sumtw[i] = (TH1D*)histfile->Get(("sumtw" + to_string(i)).c_str());
      for(int j=0; j<3; ++j)
	{
	  ET[i][j] = (TH1D*)histfile->Get(("et"+to_string(i)+to_string(j)).c_str());
	  TW[i][j] = (TH1D*)histfile->Get(("tw"+to_string(i)+to_string(j)).c_str());
	  for(int k=0; k<centbins; ++k)
	    {
	      centtow[i][j][k] = (TH1D*)histfile->Get(("centtow"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str());
	      centet[i][j][k] = (TH1D*)histfile->Get(("centet"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str());
	    }
	}
    }
  string cal[3] = {"EMCal","IHCal","OHCal"};
  string xlabel;
  tree->SetBranchAddress("sub",&sub);
  tree->SetBranchAddress("scale",scale);
  tree->SetBranchAddress("frac",frac);
  tree->SetBranchAddress("mine",&mine);
  tree->SetBranchAddress("zcut",&zcut);
  tree->SetBranchAddress("run",&run);
  tree->GetEvent(0);
  TCanvas* c1 = new TCanvas("c1","c1",1000,1000);
  xlabel = "MBD Z vertex [cm]";
  plotsimdat(c1, zhist, NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, 0, centbins, "zvtx", plotdir,"all/", centbins);
  xlabel = "MBD charge sum [??]";
  plotsimdat(c1, mbh[1], mbh[0], 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, 0, centbins, "mboverlay", plotdir,"all/", centbins);
  //plotsimdat(c1, mbh[1], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, 0, 20, "mbdat", plotdir);
  //plotsimdat(c1, mbh[0], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, 0, 20, "mbsim", plotdir);
  xlabel = "E_{T, event calorimeter sum}";
  plotsimdat(c1, sumev[1], sumev[0], 1, "total", scale[0], sub, run, mine, zcut, xlabel, 0, centbins, "et_total", plotdir,"all/", centbins);
  xlabel = "E_{T, stacked calorimeter towers}";
  plotsimdat(c1, sumtw[1], sumtw[0], 1, "total", scale[0], sub, run, mine, zcut, xlabel, 0, centbins, "et_sumtow", plotdir,"all/", centbins);
  
  for(int j=0; j<3; ++j)
    {
      xlabel = "E_{T," + cal[j] +" event} [GeV]";
      plotsimdat(c1, ET[1][j], ET[0][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, 0, centbins, "et_event", plotdir, "all/", centbins);
      xlabel = "E_{T," + cal[j] +" tower} [GeV]";
      plotsimdat(c1, TW[1][j], TW[0][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, 0, centbins, "et_tower", plotdir, "all/", centbins);
      for(int k=0; k<centbins; ++k)
	{
	  xlabel = "E_{T," + cal[j] +" tower} [GeV]";
	  plotsimdat(c1, centtow[1][j][k], centtow[0][j][k], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, k, k+1, "centtow", plotdir, "cent/", centbins);
	  xlabel = "E_{T," + cal[j] +" event} [GeV]";
	  plotsimdat(c1, centet[1][j][k], centet[0][j][k], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, k, k+1, "centet", plotdir, "cent/", centbins);
	}
    }
  
  return 0;
}

int plot()
{
  const int nfiles = 1;
  string filenames[nfiles] = {
    "savedhists_subtr_0_minE_0_scale_1.30_zcut_30_run_21615.root"
    /*
    "savedhists_subtr_0_minE_0_scale_1.30_zcut_30_run_21615.root"
    "savedhists_subtr_0_minE_0_scale_1.50_zcut_30_run_21615.root",
    "savedhists_subtr_0_minE_0_scale_1.30_zcut_30_run_21615.root",
    "savedhists_subtr_0_minE_0_scale_1.30_zcut_10_run_21615.root",
    "savedhists_subtr_0_minE_0_scale_1.00_zcut_30_run_21615.root",
    "savedhists_subtr_0_minE_0_scale_1.00_zcut_10_run_21615.root",
    "savedhists_subtr_0_minE_5_scale_1.30_zcut_30_run_21615.root",
    "savedhists_subtr_0_minE_5_scale_1.30_zcut_10_run_21615.root",
    "savedhists_subtr_0_minE_5_scale_1.00_zcut_30_run_21615.root",
    "savedhists_subtr_0_minE_5_scale_1.00_zcut_10_run_21615.root",
    "savedhists_subtr_18_minE_0_scale_1.30_zcut_30_run_21615.root",
    "savedhists_subtr_18_minE_0_scale_1.30_zcut_10_run_21615.root",
    "savedhists_subtr_18_minE_0_scale_1.00_zcut_30_run_21615.root",
    "savedhists_subtr_18_minE_0_scale_1.00_zcut_10_run_21615.root",
    "savedhists_subtr_18_minE_5_scale_1.30_zcut_30_run_21615.root",
    "savedhists_subtr_18_minE_5_scale_1.30_zcut_10_run_21615.root",
    "savedhists_subtr_18_minE_5_scale_1.00_zcut_30_run_21615.root",
    "savedhists_subtr_18_minE_5_scale_1.00_zcut_10_run_21615.root",
    */
  };

  for(int i=0; i<nfiles; ++i)
    {
      int dummy = called_plot(filenames[i]);
    }
  return 0;
}
