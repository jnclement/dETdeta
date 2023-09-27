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
#include <TSystem.h>
#include "dlUtility.h"
#include "mbd_info.h"
#include "/home/jocl/Documents/main/physics/projects/sphenix_macros/macros/macros/sPHENIXStyle/sPhenixStyle.h"
#include "/home/jocl/Documents/main/physics/projects/sphenix_macros/macros/macros/sPHENIXStyle/sPhenixStyle.C"

void multiplot(string options, TCanvas* ca, TH1D** hists, int logy, string cal, float sc, float sub, int run, float mine, float zcut, string xlabel, string ylabel, int percent0, int percent1, string name, string dir, string subdir, int centbins, int datorsim, int calnum)
{
  string typ = "";
  if(datorsim) typ = "Data";
  else typ = "HIJING";
  int kcodes[9] = {0,-4,-7,-9,-10,-8,-5,-1,4};
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
  const int ntext = 4;
  string texts[ntext];
  string ztext = ((zcut > 100)?"No z cut,":"|z|<"+params[3]+" cm,");
  texts[0] = params[1] + " MeV subtracted from each tower";
  texts[2] = ztext +" min tower E = " +params[2] + " MeV";
  texts[1] = "Run " + to_string(run) + " " + to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  if(sc > 1.) texts[3] = "HIJING scaled by " + params[0];
  else texts[3] = "";
  auto leg = new TLegend(0.,0.,0.15,0.13);
  leg->SetTextSize(0.015);
  for(int i=0; i<centbins; ++i)
    {
      hists[i]->SetMarkerColor(kRed+kcodes[i]);
      leg->AddEntry(hists[i],(typ+" " + to_string((centbins-i-1)*10)+"-"+to_string((centbins-i)*10) + "% centrality").c_str(),"P");
    }
  const int ntext2 = 3;
  string texts2[ntext2];
  texts2[0] = "|z|<" + params[3] + " cm, min tower E " + params[2] + " MeV";
  texts2[1] = params[1] + " MeV subtracted from each tower";
  texts2[2] = "Run " + to_string(run) + " " + to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  ca->cd();
  if(logy) gPad->SetLogy();
  else gPad->SetLogy(0);
  gPad->SetTicks(1);
  float maxval = -999999;
  for(int i=0; i<centbins; ++i)
    {
      maxval = max(hists[i]->GetMaximum(),maxval);
    }
  float minval = 999999;
  
  for(int i=0; i<centbins; ++i)
    {
      for(int j=0; j<hists[i]->GetNbinsX(); ++j)
	{
	  minval = min(minval, hists[i]->GetBinContent(j+1));
	}
    }
  if(minval < 0 && logy) minval = 1E-10;
  if(calnum == 0) maxval = 60;
  else if(calnum == 1) maxval = 6;
  else if(calnum == 2) maxval = 20;
  if(logy) hists[0]->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  else hists[0]->GetYaxis()->SetRangeUser(min(0,minval)-abs(min(0,minval))/10.,maxval+abs(maxval)/10.);
  hists[0]->GetXaxis()->SetTitle(xlabel.c_str());
  hists[0]->GetYaxis()->SetTitle(ylabel.c_str());
  hists[0]->GetYaxis()->SetLabelSize(0.025);
  hists[0]->GetXaxis()->SetLabelSize(0.025);
  hists[0]->Draw(options.c_str());
  for(int i=0; i<centbins; ++i) hists[i]->Draw(("SAME "+options).c_str());
  sphenixtext();
  multitext(texts, ntext, 0.2);
  drawText(typ.c_str(), 0.15, 0.96, 0, kBlack, 0.04);
  leg->Draw();
  ca->SaveAs((dir+"pdf/"+subdir+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
  ca->SaveAs((dir+"png/"+subdir+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
}

void plotsimdat(string options, TCanvas* ca, TH1* dathist, TH1* simhist, int logy, string cal, float sc, float sub, int run, float mine, float zcut, string xlabel, string ylabel, int percent0, int percent1, string name, string dir, string subdir, int centbins, int norm)
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
  string ztext = ((zcut > 100)?"No z cut,":"|z|<"+params[3]+" cm,");
  if(norm) texts[0] =  "Histogram areas normed to 1";
  else texts[0] = "";
  texts[2] = ztext +" min tower E = " +params[2] + " MeV";
  texts[1] = params[1] + " MeV subtracted from each tower";
  texts[3] = "Run " + to_string(run) + " " + to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  texts[4] = "HIJING scaled by " + params[0];
  auto leg = new TLegend(0.15,0.96,0.25,0.995);
  leg->SetTextSize(0.02);
  if(simhist) leg->AddEntry(simhist,"HIJING","P");
  if(simhist) leg->AddEntry(dathist,"Data","P");
  const int ntext2 = 3;
  string texts2[ntext2];
  texts2[0] = "|z|<" + params[3] + " cm, min tower E " + params[2] + " MeV";
  texts2[1] = params[1] + " MeV subtracted from each tower";
  texts2[2] = "Run " + to_string(run) + " " + to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  ca->cd();
  if(logy) gPad->SetLogy();
  else gPad->SetLogy(0);
  gPad->SetTicks(1);
  if(norm) dathist->Scale(1./dathist->Integral());
  if(simhist && norm) simhist->Scale(1./simhist->Integral());
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
  float mintemp[2];

  if(logy)
    {
      if(simhist)
	{
	  mintemp[0] = min(dathist->GetBinContent(dathist->FindLastBinAbove(0,1)),simhist->GetBinContent(simhist->FindLastBinAbove(0,1)));
	  mintemp[1] = min(dathist->GetBinContent(dathist->FindFirstBinAbove(0,1)),simhist->GetBinContent(simhist->FindFirstBinAbove(0,1)));
	  minval = min(mintemp[0],mintemp[1]);
	}
      else minval = min(dathist->GetBinContent(dathist->FindFirstBinAbove(0,1)),dathist->GetBinContent(dathist->FindLastBinAbove(0,1)));
    }
  else
    {
      if(simhist) minval = min(dathist->GetMinimum(),simhist->GetMaximum());
      else minval = dathist->GetMinimum();
    }
  if(logy) dathist->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  else dathist->GetYaxis()->SetRangeUser(min(0.,minval)-abs(min(0,minval))/10.,maxval+abs(maxval)/10.);
  dathist->GetXaxis()->SetTitle(xlabel.c_str());
  dathist->GetYaxis()->SetTitle(ylabel.c_str());
  dathist->GetYaxis()->SetLabelSize(0.025);
  dathist->GetXaxis()->SetLabelSize(0.025);
  dathist->Draw(options.c_str());
  if(simhist) simhist->Draw(("SAME "+options).c_str());
  sphenixtext();
  multitext(texts, ntext, 0.03, 0.125);
  if(simhist) leg->Draw();
  ca->SaveAs((dir+"pdf/"+subdir+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
  ca->SaveAs((dir+"png/"+subdir+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
  if(simhist)
    {
      dathist->Divide(simhist);
      dathist->GetYaxis()->SetTitle("Data/HIJING");
      maxval = dathist->GetMaximum();
      minval = dathist->GetMinimum();
      gPad->SetLogy();
      dathist->GetYaxis()->SetRangeUser(0.01,100);
      dathist->Draw(options.c_str());
      sphenixtext();
      multitext(texts2, ntext2);
      ca->SaveAs((dir+"pdf/"+"ratio_"+subdir+"ratio_"+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
      ca->SaveAs((dir+"png/"+"ratio_"+subdir+"ratio_"+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
    }
}

int called_plot(string histfilename = "savedhists_fracsim_1_fracdat_1_subtr_0_minE_0_scale_1.30_zcut_30_run_21615.root", string treename = "ttree", string datdir = "/home/jocl/datatemp/", string plotdir = "/home/jocl/datatemp/plots/")
{   
  //gSystem->RedirectOutput("test.txt");
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  SetsPhenixStyle();
  const int centbins = 9;
  TH1D* means[2][3];
  TH1D* sigs[2][3];
  TH1D* meandiffnoavg[3];
  for(int i=0; i<3; ++i)
    {
      meandiffnoavg[i] = new TH1D(("meandiffnoavg"+to_string(i)).c_str(),"",centbins,0,90);
    }
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  means[i][j] = new TH1D(("mean"+to_string(i)+to_string(j)).c_str(),"",centbins,0,90);
	  sigs[i][j] = new TH1D(("sig"+to_string(i)+to_string(j)).c_str(),"",centbins,0,90);
	}
    }
  TFile* histfile = TFile::Open((datdir+histfilename).c_str());
  TTree* tree = histfile->Get<TTree>(treename.c_str());
  TH1D* centtow[2][3][centbins];
  TH1D* centet[2][3][centbins];
  TH1D* ET[2][3];
  TH1D* dET[2][3];
  TH1D* dETcent[2][3][centbins];
  TH1D* TW[2][3];
  TH1D* sumev[2];
  TH1D* sumtw[2];
  TH1D* mbh[2];
  TH1D* zhist;
  TH1D* meandiff[3];
  TH1D* sigmu[2][3];
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
	  if(i==0) meandiff[j] = (TH1D*)histfile->Get(("md"+to_string(j)).c_str());
	  sigmu[i][j] = (TH1D*)histfile->Get(("sigmu"+to_string(i)+to_string(j)).c_str());
	  ET[i][j] = (TH1D*)histfile->Get(("et"+to_string(i)+to_string(j)).c_str());
	  TW[i][j] = (TH1D*)histfile->Get(("tw"+to_string(i)+to_string(j)).c_str());
	  dET[i][j] = (TH1D*)histfile->Get(("dET"+to_string(i)+to_string(j)).c_str());
	  for(int k=0; k<centbins; ++k)
	    {
	      centtow[i][j][k] = (TH1D*)histfile->Get(("centtow"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str());
	      centet[i][j][k] = (TH1D*)histfile->Get(("centet"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str());
	      dETcent[i][j][k] = (TH1D*)histfile->Get(("dETcent"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str());
	    }
	}
    }
  cout << "Finished getting hists" << endl;
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  for(int k=0; k<centbins; ++k)
	    {
	      means[i][j]->SetBinContent(centbins-k,centet[i][j][k]->GetMean());
	      means[i][j]->SetBinError(centbins-k,0);
	      sigs[i][j]->SetBinContent(centbins-k,centet[i][j][k]->GetStdDev());
	      sigs[i][j]->SetBinError(centbins-k,0);
	    }
	}
    }
  for(int i=0; i<3; ++i)
    {
      for(int j=0; j<centbins; ++j)
	{
	  meandiffnoavg[i]->SetBinContent(centbins-j,centet[1][i][j]->GetMean()-centet[0][i][j]->GetMean());
	  meandiffnoavg[i]->SetBinError(centbins-j,0);
	}
    }
  cout << "Finished filling new hists" << endl;
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
  string ylabel = "Counts";
  string options = "";
  cout << "Starting to plot" << endl;
  plotsimdat(options, c1, zhist, NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "zvtx", plotdir,"all/", centbins, 1);
  xlabel = "MBD charge sum [??]";
  plotsimdat(options, c1, mbh[1], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "mboverlay", plotdir,"all/", centbins, 1);
  //plotsimdat(options, c1, mbh[1], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, 0, 20, "mbdat", plotdir);
  //plotsimdat(options, c1, mbh[0], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, 0, 20, "mbsim", plotdir);
  xlabel = "E_{T, event calorimeter sum}";
  plotsimdat(options, c1, sumev[1], sumev[0], 1, "total", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_total", plotdir,"all/", centbins, 1);
  xlabel = "E_{T, stacked calorimeter towers}";
  plotsimdat(options, c1, sumtw[1], sumtw[0], 1, "total", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_sumtow", plotdir,"all/", centbins, 1);
  cout << "Plotted non-loops" << endl;
  for(int j=0; j<3; ++j)
    {
      options = "p";
      xlabel = "MBD Centrality [%]";
     ylabel = "#mu_{E_{T} event data}^{"+cal[j]+"}-#mu_{E_{T} event HIJING}^{"+cal[j]+"} [GeV]";
      plotsimdat(options, c1, meandiffnoavg[j], NULL, 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "meandiffnoavg", plotdir, "all/", centbins, 0);
	
      options = "p";
      ylabel = "#mu_{E_{T} event}^{"+cal[j]+"} [GeV]";
      xlabel = "Centrality MBD [%]";
      plotsimdat(options, c1, means[1][j], means[0][j], 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "meandist", plotdir, "all/", centbins, 0);
      ylabel = "#sigma_{E_{T} event}^{"+cal[j]+"} [GeV]";
      plotsimdat(options, c1, sigs[1][j], sigs[0][j], 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "sigdist", plotdir, "all/", centbins, 0);
      
      options = "hist p";
      ylabel = "2(#mu_{data}^{" + cal[j] + "}-#mu_{HIJING}^{" + cal[j] + "})/(#mu_{data}^{" + cal[j] + "}+#mu_{HIJING}^{" + cal[j] + "})";
      xlabel = "MBD Centrality [%]";
      meandiff[j]->Scale(2.);
      plotsimdat(options, c1, meandiff[j], NULL, 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "meandiff", plotdir,"all/", centbins, 0);
      ylabel = "#sigma_{"+cal[j]+"}/#mu_{"+cal[j]+"}";
      plotsimdat(options, c1, sigmu[1][j], sigmu[0][j], 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "sigmu", plotdir, "all/", centbins, 0);
      ylabel = "Counts";
      options = "";
      xlabel = "E_{T," + cal[j] +" event} [GeV]";
      plotsimdat(options, c1, ET[1][j], ET[0][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_event", plotdir, "all/", centbins, 1);
      xlabel = "E_{T," + cal[j] +" tower} [GeV]";
      plotsimdat(options, c1, TW[1][j], TW[0][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_tower", plotdir, "all/", centbins, 1);
      if(j==0) xlabel = "floor{(#eta bin " +cal[j]+")/4}";
      else xlabel = "#eta bin " + cal[j];
      ylabel = "dE_{T}/d#eta [GeV]";
      options = "p";
      plotsimdat(options, c1, dET[1][j],dET[0][j],0,cal[j], scale[0],sub,run,mine,zcut,xlabel,ylabel,0,centbins,"det",plotdir,"all/",centbins, 0);

      
      multiplot(options, c1, dETcent[0][j], 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "detcent_sim", plotdir, "all/", centbins, 0, j);
      multiplot(options, c1, dETcent[1][j], 0, cal[j], scale[1], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "detcent_data", plotdir, "all/", centbins, 1, j);
      for(int k=0; k<centbins; ++k)
	{
	  //outputon(gSystem);
	  cout << centet[1][j][k]->GetMean() << " " << dETcent[1][j][k]->Integral() << " " <<centet[0][j][k]->GetMean() << " " << dETcent[0][j][k]->Integral() << endl;
	  //outputoff(gSystem,"test.txt");
	  options = "";
	  ylabel = "Counts";
	  xlabel = "E_{T," + cal[j] +" tower} [GeV]";
	  plotsimdat(options, c1, centtow[1][j][k], centtow[0][j][k], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, k, k+1, "centtow", plotdir, "cent/", centbins, 1);
	  xlabel = "E_{T," + cal[j] +" event} [GeV]";
	  plotsimdat(options, c1, centet[1][j][k], centet[0][j][k], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, k, k+1, "centet", plotdir, "cent/", centbins, 1);
	  options = "hist p";
	  xlabel = "#eta bin";
	  ylabel = "dE_{T}/d#eta [GeV]";
	  plotsimdat(options, c1, dETcent[1][j][k], dETcent[0][j][k], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel,k,k+1,"detcent",plotdir,"cent/",centbins, 0);
	  options = "";
	}
    }
  
  return 0;
}

int plot()
{
  const int nfiles = 1;
  string filenames[nfiles] = {
    "savedhists_fracsim_12_fracdat_50_subtr_0_minE_0_scale_1.30_zcut_30_run_21615_ntc.root"
  };

  for(int i=0; i<nfiles; ++i)
    {
      int dummy = called_plot(filenames[i]);
    }
  return 0;
}
