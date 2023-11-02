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
#include <TFitResult.h>
#include "TLatex.h"
#include "stdlib.h"
#include <fstream>
#include <cstdlib>
#include <TF1.h>
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

void centoverlayplot(string options, TCanvas* ca, TH1D* mainhist, TH1D** hists, int logy, string cal, float sc, float sub, int run, float mine, float zcut, string xlabel, string ylabel, int percent0, int percent1, string name, string dir, string subdir, int centbins, int datorsim, int calnum)
{
  string typ = "";
  int kcodes[3] = {0,1,2};
  int knames[6] = {kBlue,kMagenta,kRed,kYellow,kGreen,kCyan};
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
  if(datorsim) typ = "Run-23 Au+Au";
  else typ = "HIJING scaled by " + params[0];
  if(calnum==4) typ = "";
  const int ntext = 4;
  string texts[ntext];
  string ztext = ((zcut > 100)?"No z cut,":"|z|<"+params[3]+" cm");
  texts[0] = (calnum==4?"":params[1] + " MeV subtracted from each tower");
  texts[2] = ztext +(calnum==4?", min particle E = " +params[2] + " MeV":", min tower E = " +params[2] + " MeV");
  texts[1] = "Run " + to_string(run) + " " + to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  //if(sc > 1.) texts[3] = "HIJING scaled by " + params[0];
  texts[3] = "";
  TLegend* leg;
  if(calnum!=4) leg = new TLegend(0.75,0.7,0.9,0.9);
  else leg = new TLegend(0.85,0.65,0.9,0.9);
  leg->SetTextSize(0.015);
  mainhist->SetMarkerColor(kBlack);
  leg->AddEntry(mainhist,("0-90% centrality"),"P");
  if(calnum==0 || calnum == 3 || calnum==4) mainhist->Rebin(5);
  for(int i=0; i<centbins; ++i)
    {
      if(calnum == 0 || calnum == 3 || calnum==4) hists[i]->Rebin(5);
      hists[i]->SetMarkerColor(knames[i%6]+kcodes[i/6]);
      hists[i]->SetLineColor(knames[i%6]+kcodes[i/6]);
      hists[i]->SetFillColorAlpha(knames[i%6]+kcodes[i/6],0.2);
      leg->AddEntry(hists[i],(to_string((centbins-i-1)*centrange)+"-"+to_string((centbins-i)*centrange) + "% centrality").c_str(),"P");
    }
  ca->cd();
  if(logy) gPad->SetLogy();
  else gPad->SetLogy(0);
  gPad->SetTicks(1);
  float maxval = -999999;
  maxval = mainhist->GetMaximum();
  for(int i=0; i<centbins; ++i)
    {
      maxval = max(hists[i]->GetMaximum(),maxval);
    }
  float minval = 999999;
  for(int i=0; i<centbins; ++i)
    {
      for(int j=0; j<hists[i]->GetNbinsX(); ++j)
	{
	  if(hists[i]->GetBinContent(j+1) > 0) minval = min(minval, hists[i]->GetBinContent(j+1));
	}
    }
  if(minval <= 0 && logy) minval = 1E-10;
  cout << minval << " " << maxval << endl;
  //if(calnum == 0) maxval = 10;
  //else if(calnum == 1) maxval = 6;
  //else if(calnum == 2) maxval = 20;
  if(logy) mainhist->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  else mainhist->GetYaxis()->SetRangeUser(min(0,minval)-abs(min(0,minval))/10.,maxval+abs(maxval)/10.);
  float ranges[5] = {1500,150,400,2000,5000};
  mainhist->GetXaxis()->SetRangeUser(0,ranges[calnum]);
  cout << mainhist->GetName() << endl;
  mainhist->GetXaxis()->SetTitle(xlabel.c_str());
  mainhist->GetYaxis()->SetTitle(ylabel.c_str());
  mainhist->GetYaxis()->SetLabelSize(0.025);
  mainhist->GetXaxis()->SetLabelSize(0.025);
  mainhist->Draw(options.c_str());
  mainhist->SetMarkerColor(kBlack);
  mainhist->Draw(options.c_str());
  for(int i=0; i<centbins; ++i) hists[i]->Draw(("SAME "+options).c_str());
  sphenixtext();
  //multitext(texts, ntext, 0.1);
  drawText(typ.c_str(), 0.15, 0.96, 0, kBlack, 0.04);
  leg->Draw();
  ca->SaveAs((dir+"pdf/"+subdir+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
  ca->SaveAs((dir+"png/"+subdir+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
}


void multiplot(string options, TCanvas* ca, TH1D** hists, int logy, string cal, float sc, float sub, int run, float mine, float zcut, string xlabel, string ylabel, int percent0, int percent1, string name, string dir, string subdir, int centbins, int datorsim, int calnum)
{
  string typ = "";
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
  if(datorsim) typ = "Run-23 Au+Au";
  else typ = "HIJING scaled by " + params[0];
  const int ntext = 4;
  string texts[ntext];
  string ztext = ((zcut > 100)?"No z cut,":"|z|<"+params[3]+" cm,");
  texts[0] = params[1] + " MeV subtracted from each tower";
  texts[2] = ztext +" min tower E = " +params[2] + " MeV";
  texts[1] = "Run " + to_string(run) + " " + to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  //if(sc > 1.) texts[3] = "HIJING scaled by " + params[0];
  texts[3] = "";
  auto leg = new TLegend(0.,0.,0.15,0.13);
  leg->SetTextSize(0.015);
  for(int i=0; i<centbins; ++i)
    {
      hists[i]->SetMarkerColor(kRed+kcodes[i]);
      leg->AddEntry(hists[i],(typ+" " + to_string((centbins-i-1)*centrange)+"-"+to_string((centbins-i)*centrange) + "% centrality").c_str(),"P");
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
  if(calnum == 0) maxval = 10;
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
  //multitext(texts, ntext, 0.2);
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
  const int ntext = 4;
  string texts[ntext];
  string ztext = ((zcut > 100)?"No z cut,":"|z|<"+params[3]+" cm,");
  if(norm) texts[0] =  "Histogram areas normed to 1";
  else texts[0] = "";
  texts[2] = ztext +" min tower E = " +params[2] + " MeV";
  texts[1] = params[1] + " MeV subtracted from each tower";
  if(name!="det_no_acc_cor" && name != "truthpar_et" && name != "truth_reco_et") texts[3] = "Run " + to_string(run) + " " + to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  else texts[3] = to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  //texts[4] = "HIJING scaled by " + params[0];
  //if(cal=="MBD") texts[4] = "";
  auto leg = new TLegend(0.6,0.23,0.9,0.3);
  leg->SetTextSize(0.02);
  if(simhist) leg->AddEntry(simhist,"HIJING","P");
  if(name != "fullcor" && name != "det_no_acc_cor" && name != "truthpar_et") leg->AddEntry(dathist,"Data","P");
  const int ntext2 = 4;
  string texts2[ntext2];
  texts2[0] = "|z|<" + params[3] + " cm, min tower E " + params[2] + " MeV";
  texts2[1] = params[1] + " MeV subtracted from each tower";
  if(name!="det_no_acc_cor" && name != "truthpar_et" && name != "truth_reco_et") texts2[2] = "Run " + to_string(run) + " " + to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  else texts2[2] = to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+"% centrality";
  if(name!="detcent") texts2[3] = "No acceptance cut applied";
  ca->cd();
  if(logy) gPad->SetLogy();
  else gPad->SetLogy(0);
  gPad->SetTicks(1);
  if(norm) dathist->Scale(1./dathist->Integral());
  if(simhist && norm) simhist->Scale(1./simhist->Integral());
  int dcolor;
  int scolor;
  if(dir=="/home/jocl/datatemp/plots_unc/")
    {
      dcolor = kGreen+2;
      scolor = kMagenta+2;
    }
  else
    {
      dcolor = kRed;
      scolor = kBlue;
    }
  dathist->SetLineColor(dcolor);
  dathist->SetMarkerColor(dcolor);
  if(simhist)
    {
      simhist->SetLineColor(scolor);
      simhist->SetMarkerColor(scolor);
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
  //else dathist->GetYaxis()->SetRangeUser(min(0.,minval)-abs(min(0,minval))/10.,maxval+abs(maxval)/10.);
  else dathist->GetYaxis()->SetRangeUser(minval*0.9,maxval*1.1);
  if(name=="sigmu")
    {
      dathist->GetYaxis()->SetRangeUser(0,2);
    }
  if(cal == "MBD") dathist->GetXaxis()->SetRangeUser(-zcut,zcut);
  dathist->GetXaxis()->SetTitle(xlabel.c_str());
  dathist->GetYaxis()->SetTitle(ylabel.c_str());
  dathist->GetYaxis()->SetLabelSize(0.025);
  dathist->GetXaxis()->SetLabelSize(0.025);
  TFitResultPtr hp[2];
  if((simhist && (name=="zcent")) || name=="zvtx")
    {
      hp[1] = dathist->Fit("gaus","S");
      if(simhist) hp[0] = simhist->Fit("gaus","S");
      dathist->GetFunction("gaus")->SetLineColor(dcolor);
      if(simhist) simhist->GetFunction("gaus")->SetLineColor(scolor);
      dathist->Draw(options.c_str());
      if(simhist) simhist->Draw(("SAME "+options).c_str());
      dathist->GetFunction("gaus")->Draw("SAME");
      if(simhist) simhist->GetFunction("gaus")->Draw("SAME");
    }
  TH1D* simflip;
  if(simhist) simflip = (TH1D*)simhist->Clone();
  TH1D* datflip = (TH1D*)dathist->Clone();
  int etabins = 20;
  for(int i=0; i<etabins; ++i)
    {
      if(simhist) simflip->SetBinContent(i+1,simhist->GetBinContent(etabins-i));
      if(simhist) simflip->SetBinError(i+1,simhist->GetBinError(etabins-i));
      datflip->SetBinContent(i+1,dathist->GetBinContent(etabins-i));
      datflip->SetBinError(i+1,dathist->GetBinError(etabins-i));
    }
  if(simhist)
    {
      dathist->Draw(options.c_str());
      simhist->Draw(("SAME "+options).c_str());
      datflip->SetMarkerStyle(4);
      leg->AddEntry(datflip,"Data reversed over #eta=0","P");
      datflip->Draw(("SAME "+options).c_str());
      simflip->SetMarkerStyle(4);
      leg->AddEntry(simflip,"HIJING reversed over #eta=0","P");
      simflip->Draw(("SAME "+options).c_str());
    }
  else
    {
      dathist->Draw(options.c_str());
      datflip->SetMarkerStyle(4);
      //leg->AddEntry(dathist,"Data","P");
      if(name == "det_no_acc_cor" || name == "truthpar_et") leg->AddEntry(dathist,"HIJING","P");
      if(name != "fullcor" && name != "det_no_acc_cor" && name!="truthpar_et") leg->AddEntry(datflip,"Data reversed over #eta=0","P");
      else if(name == "det_no_acc_cor" || name == "truthpar_et") leg->AddEntry(datflip,"HIJING reversed over #eta=0","P");
      if(name != "fullcor") datflip->Draw(("SAME "+options).c_str());
    }
  sphenixtext();
  multitext(texts, ntext, 0.03, 0.11);
  if(name != "fullcor") leg->Draw();
  if((name == "zcent" && simhist) || name=="zvtx")
    {
      stringstream stmean, stsig;
      stmean << std::fixed << std::setprecision(2) << hp[1]->Parameter(1);
      stsig << std::fixed << std::setprecision(2) << hp[1]->Parameter(2);
      drawText(("#mu"+(string)(simhist?"_{data}":"")+"="+stmean.str()+", #sigma"+(simhist?"_{data}":"")+"="+stsig.str()).c_str(),0.9,0.91,1,kBlack,0.025);
      stmean.str("");
      stsig.str("");
      if(simhist)
	{
	  stmean << std::fixed << std::setprecision(2) << hp[0]->Parameter(1);
	  stsig << std::fixed << std::setprecision(2) << hp[0]->Parameter(2);
	  drawText(("#mu_{HIJING}="+stmean.str()+", #sigma_{HIJING}="+stsig.str()).c_str(),0.9,0.88,1,kBlack,0.025);
	}
    }
  if(name=="det_no_acc_cor")
    {
      drawText("No acceptance cut applied",0.7,0.085,1,kBlack,0.025);
    }
  ca->SaveAs((dir+"pdf/"+subdir+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
  ca->SaveAs((dir+"png/"+subdir+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
  if(simhist)
    {
      dathist->Divide(simhist);
      if(name == "detcent") dathist->GetYaxis()->SetTitle((cal+" dE_{T}/d#eta Data/HIJING").c_str());
      maxval = dathist->GetMaximum();
      minval = dathist->GetMinimum();
      gPad->SetLogy(0);
      dathist->GetYaxis()->SetRangeUser(0,2);
      if(name == "detcent") dathist->GetYaxis()->SetRangeUser(0.8,1.2);
      if(name == "truth_reco_et" && cal=="EMCal") dathist->GetYaxis()->SetRangeUser(0.65,0.75);
      else if(name == "truth_reco_et") dathist->GetYaxis()->SetRangeUser(0,0.2);
      dathist->Draw(options.c_str());
      sphenixtext();
      multitext(texts2, ntext2);
      ca->SaveAs((dir+"pdf/"+"ratio_"+subdir+"ratio_"+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
      ca->SaveAs((dir+"png/"+"ratio_"+subdir+"ratio_"+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
    }
}

int called_plot(string histfilename = "savedhists_fracsim_1_fracdat_1_subtr_0_minE_0_scale_1.00_zcut_30_run_21615_20231018_nopileup_cor.root", string treename = "ttree", string datdir = "/home/jocl/datatemp/", string plotdir = "/home/jocl/datatemp/plots/")
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
  TH1D* fullcor[3][centbins];
  TH1D* TW[2][3];
  TH1D* sumev[2];
  TH1D* sumtw[2];
  TH1D* mbh[2];
  TH1D* zhist[2];
  TH1D* meandiff[3];
  TH1D* sigmu[2][3];
  TH1D* zcent[2][centbins];
  TH2D* deadmap[2][3][centbins];
  TH1D* ettotcent[2][centbins];
  TH1D* truthparehist;
  TH1D* truthparnhist;
  TH1D* truthparecent[centbins];
  TH1D* truthparncent[centbins];
  TH1D* truthpar_et[centbins];
  TH1D* dETcentsimunc[3][centbins];
  float sub;
  float scale[2];
  int frac[2];
  float mine;
  float zcut;
  int run;
  truthparehist = (TH1D*)histfile->Get("truthparehist");
  truthparnhist = (TH1D*)histfile->Get("truthparnhist");
  mbh[0] = (TH1D*)histfile->Get("smbh");
  mbh[1] = (TH1D*)histfile->Get("dmbh");
  zhist[0] = (TH1D*)histfile->Get("zhist_0");
  zhist[1] = (TH1D*)histfile->Get("zhist_1");
  for(int i=0; i<2; ++i)
    {
      sumev[i] = (TH1D*)histfile->Get(("sumev" + to_string(i)).c_str());
      sumtw[i] = (TH1D*)histfile->Get(("sumtw" + to_string(i)).c_str());
      for(int j=0; j<3; ++j)
	{
	  means[i][j] = (TH1D*)histfile->Get(("meancent"+to_string(i)+to_string(j)).c_str());
	  if(i==0) meandiff[j] = (TH1D*)histfile->Get(("md"+to_string(j)).c_str());
	  sigmu[i][j] = (TH1D*)histfile->Get(("sigmu"+to_string(i)+to_string(j)).c_str());
	  ET[i][j] = (TH1D*)histfile->Get(("et"+to_string(i)+to_string(j)).c_str());
	  TW[i][j] = (TH1D*)histfile->Get(("tw"+to_string(i)+to_string(j)).c_str());
	  dET[i][j] = (TH1D*)histfile->Get(("dET"+to_string(i)+to_string(j)).c_str());
	  for(int k=0; k<centbins; ++k)
	    {
	      if(j==0)
		{
		  zcent[i][k] = (TH1D*)histfile->Get(("zcent"+to_string(i)+"_"+to_string(k)).c_str());
		  ettotcent[i][k] = (TH1D*)histfile->Get(("ettotcent"+to_string(i)+"_"+to_string(k)).c_str());
		}
	      
	      deadmap[i][j][k] = (TH2D*)histfile->Get(("deadmap"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str());
	      centtow[i][j][k] = (TH1D*)histfile->Get(("centtow"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str());
	      centet[i][j][k] = (TH1D*)histfile->Get(("centet"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str());
	      dETcent[i][j][k] = (TH1D*)histfile->Get(("dETcent"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str());
	    }
	}
    }
  TCanvas* c2 = new TCanvas("","");
  c2->cd();
  cout << centet[1][0][17]->GetEntries() << endl;
  centet[1][0][17]->Draw();
  c2->SaveAs("test.png");
  cout << "Finished getting hists" << endl;
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  for(int k=0; k<centbins; ++k)
	    {
	      sigs[i][j]->Divide(sigmu[i][j],means[i][j]);
	    }
	}
    }
  for(int i=0; i<3; ++i)
    {
      for(int j=0; j<centbins; ++j)
	{
	  if(i==0)
	    {
	      truthparecent[j] = (TH1D*)histfile->Get(("truthparecent_"+to_string(j)).c_str());
	      truthparncent[j] = (TH1D*)histfile->Get(("truthparncent_"+to_string(j)).c_str());
	      truthpar_et[j] = (TH1D*)histfile->Get(("truthpar_et_"+to_string(j)).c_str());
	    }
	  dETcentsimunc[i][j] = (TH1D*)histfile->Get(("dETcentsimunc_"+to_string(i)+"_"+to_string(j)).c_str());
	  cout << ("fullcor_"+to_string(i)+to_string(j)).c_str() << endl;
	  fullcor[i][j] = (TH1D*)histfile->Get(("fullcor_"+to_string(i)+to_string(j)).c_str());
	  cout << fullcor[i][j] << endl;
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
  cout << "Starting to plot" << endl;
  string options = "hist";
  xlabel = "E_{T, event}^{all} [GeV]";
  string ylabel = "Counts";
  //centoverlayplot(options, c1, sumev[0], ettotcent[0], 1, "all", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "cent_overlay_sim", plotdir, "all/", centbins, 0, 3);
  //centoverlayplot(options, c1, sumev[1], ettotcent[1], 1, "all", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "cent_overlay_dat", plotdir, "all/", centbins, 1, 3);
  xlabel = "Truth Particle N";
  //centoverlayplot(options, c1, truthparnhist, truthparncent, 1, "truth particles", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "n_cent_overlay",plotdir, "all/", centbins, 0, 4);
  xlabel = "Truth Z Vertex [cm]";
  ylabel = "Normalized Counts";
  options = "";
  plotsimdat(options, c1, zhist[0], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "zvtx", plotdir,"all/", centbins, 1);
  xlabel = "MBD charge sum [??]";
  //plotsimdat(options, c1, mbh[1], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "mboverlay_dat", plotdir,"all/", centbins, 1);
  //plotsimdat(options, c1, mbh[0], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "mboverlay_sim", plotdir,"all/", centbins, 1);
  //plotsimdat(options, c1, mbh[1], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, 0, 20, "mbdat", plotdir);
  //plotsimdat(options, c1, mbh[0], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, 0, 20, "mbsim", plotdir);
  xlabel = "E_{T, event calorimeter sum}";
  //plotsimdat(options, c1, sumev[1], sumev[0], 1, "total", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_total", plotdir,"all/", centbins, 1);
  xlabel = "E_{T, stacked calorimeter towers}";
  //plotsimdat(options, c1, sumtw[1], sumtw[0], 1, "total", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_sumtow", plotdir,"all/", centbins, 1);
  cout << "Plotted non-loops" << endl;
		  
  for(int j=0; j<3; ++j)
    {

      xlabel = "#eta";
      ylabel = cal[j]+" Corrected dE_{T}/d#eta [GeV]";
      options="p";
      plotsimdat(options, c1, fullcor[j][centbins-1], NULL, 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, centbins-1, centbins, "fullcor", plotdir, "cent/", centbins, 0);

      ylabel = "HIJING "+cal[j]+" dE_{T}/d#eta [GeV]";
      if(j==0) plotsimdat(options, c1, dETcentsimunc[j][centbins-1],NULL, 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, centbins-1, centbins, "det_no_acc_cor",plotdir,"cent/",centbins,0);

      ylabel = "Truth Particle dE_{T}/d#eta [GeV]";
      if(j==0) plotsimdat(options, c1, truthpar_et[centbins-1], NULL, 0, "", scale[0], sub, run, mine, zcut, xlabel, ylabel, centbins-1,centbins, "truthpar_et", plotdir, "cent/", centbins, 0);

      ylabel = "(Reco dE_{T}/d#eta)/(Truth dE_{T}/d#eta)";
      plotsimdat(options, c1, dETcentsimunc[j][centbins-1],truthpar_et[centbins-1],0,cal[j],scale[0],sub,run,mine,zcut,xlabel,ylabel,centbins-1,centbins,"truth_reco_et",plotdir,"cent/",centbins,0);
      
      options = "hist";
      xlabel = "#Sigma #it{E}_{T}^{"+cal[j]+"} [GeV]";
      ylabel = "Counts";
      centoverlayplot(options, c1, ET[1][j], centet[1][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "cent_overlay_dat", plotdir, "all/", centbins, 1, j);
      //centoverlayplot(options, c1, ET[0][j], centet[0][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "cent_overlay_sim", plotdir, "all/", centbins, 0, j);
      options = "p";
      xlabel = "MBD Centrality [%]";
      ylabel = "#mu_{E_{T} event data}^{"+cal[j]+"}-#mu_{E_{T} event HIJING}^{"+cal[j]+"} [GeV]";
      //plotsimdat(options, c1, meandiffnoavg[j], NULL, 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "meandiffnoavg", plotdir, "all/", centbins, 0);
	
      options = "p";
      ylabel = "#mu_{E_{T} event}^{"+cal[j]+"} [GeV]";
      xlabel = "Centrality MBD [%]";
      //plotsimdat(options, c1, means[1][j], means[0][j], 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "meandist", plotdir, "all/", centbins, 0);
      ylabel = "#sigma_{E_{T} event}^{"+cal[j]+"} [GeV]";
      //plotsimdat(options, c1, sigs[1][j], sigs[0][j], 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "sigdist", plotdir, "all/", centbins, 0);
      
      options = "p";
      ylabel = "2(#mu_{data}^{" + cal[j] + "}-#mu_{HIJING}^{" + cal[j] + "})/(#mu_{data}^{" + cal[j] + "}+#mu_{HIJING}^{" + cal[j] + "})";
      xlabel = "MBD Centrality [%]";
      meandiff[j]->Scale(2.);
      //plotsimdat(options, c1, meandiff[j], NULL, 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "meandiff", plotdir,"all/", centbins, 0);
      ylabel = cal[j]+" #sigma_{E_{T}^{event}}";
      //plotsimdat(options, c1, sigmu[1][j], sigmu[0][j], 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "sigdist", plotdir, "all/", centbins, 0); //NOTE THAT THE HIST WITH VARIABLE NAME SIGMU NOW POINTS TO THE SIGMA DISTRIBUTION, AND VICE VERSA
      ylabel = "Normalized Counts";
      options = "";
      xlabel = "E_{T," + cal[j] +" event} [GeV]";
      //plotsimdat(options, c1, ET[1][j], ET[0][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_event", plotdir, "all/", centbins, 1);
      xlabel = "E_{T," + cal[j] +" tower} [GeV]";
      //plotsimdat(options, c1, TW[1][j], TW[0][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_tower", plotdir, "all/", centbins, 1);
      if(j==0) xlabel = "floor{(#eta bin " +cal[j]+")/4}";
      else xlabel = "#eta bin " + cal[j];
      ylabel = cal[j]+" dE_{T}/d#eta [GeV]";
      options = "p";
      //plotsimdat(options, c1, dET[1][j],dET[0][j],0,cal[j], scale[0],sub,run,mine,zcut,xlabel,ylabel,0,centbins,"det",plotdir,"all/",centbins, 0);

      
      //multiplot(options, c1, dETcent[0][j], 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "detcent_sim", plotdir, "all/", centbins, 0, j);
      //multiplot(options, c1, dETcent[1][j], 0, cal[j], scale[1], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "detcent_data", plotdir, "all/", centbins, 1, j);
      cout << "test" << endl;
      for(int k=0; k<centbins; ++k)
	{
	  if(j==0)
	    {
	      xlabel = "Z-vertex [cm]";
	      ylabel = "Normalized Counts";
	      options = "p";
	      plotsimdat(options, c1, zcent[1][k], zcent[0][k], 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, ylabel, k, k+1, "zcent", plotdir, "cent/", centbins, 1);
	    }
	  //outputon(gSystem);
	  cout << j << " " << k << " " << centet[1][j][k]->GetMean() << " " << dETcent[1][j][k]->Integral() << " " <<centet[0][j][k]->GetMean() << " " << dETcent[0][j][k]->Integral() << endl;
	  //outputoff(gSystem,"test.txt");
	  options = "P";
	  ylabel = "Normalized Counts";
	  xlabel = "E_{T," + cal[j] +" tower} [GeV]";
	  //plotsimdat(options, c1, centtow[1][j][k], centtow[0][j][k], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, k, k+1, "centtow", plotdir, "cent/", centbins, 1);
	  xlabel = "E_{T," + cal[j] +" event} [GeV]";
	  //plotsimdat(options, c1, centet[1][j][k], centet[0][j][k], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, k, k+1, "centet", plotdir, "cent/", centbins, 1);
	  options = "p";
	  xlabel = "#eta";
	  ylabel = cal[j]+" dE_{T}/d#eta [GeV]";
	  plotsimdat(options, c1, dETcent[1][j][k], dETcent[0][j][k], 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel,k,k+1,"detcent",plotdir,"cent/",centbins, 0);
	}
    }
  c1->cd();
  gPad->SetLogy(0);
  xlabel = "#eta bin";
  ylabel = "#phi bin";
  options = "COLZ";
  const int ntext = 3;
  string texts[ntext] = {"0 MeV subtracted from each tower","|z|<30 cm, min tower E = 0 MeV"};
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.15);
  for(int i=0; i<centbins; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  deadmap[1][j][i]->GetYaxis()->SetLabelSize(0.025);
	  deadmap[1][j][i]->GetXaxis()->SetLabelSize(0.025);
	  //deadmap[1][j][i]->Divide(deadmap[0][j][i]);
	  deadmap[1][j][i]->GetYaxis()->SetTitle(ylabel.c_str());
	  deadmap[1][j][i]->GetXaxis()->SetTitle(xlabel.c_str());
	  deadmap[1][j][i]->GetZaxis()->SetTitle((cal[j] +" <E^{tower}_{data}> [GeV]").c_str());
	  deadmap[0][j][i]->GetZaxis()->SetTitle((cal[j] +" <E^{tower}_{HIJING}> [GeV]").c_str());
	  deadmap[1][j][i]->GetZaxis()->SetLabelSize(0.025);
	  //deadmap[1][j][i]->GetZaxis()->SetRangeUser(0,2);
	  texts[2] = "Run 21615 " + to_string((centbins-i-1)*(90/centbins))+"-"+to_string((centbins-i)*(90/centbins))+"% centrality";
	  deadmap[1][j][i]->Draw("COLZ");
	  sphenixtext();
	  multitext(texts,ntext);
	  c1->SaveAs(("/home/jocl/datatemp/plots/png/cent/dat_deadmap_"+cal[j]+"_" +to_string((centbins-i-1)*(90/centbins))+"-"+to_string((centbins-i)*(90/centbins))+".png").c_str());
	  c1->SaveAs(("/home/jocl/datatemp/plots/pdf/cent/dat_deadmap_"+cal[j]+"_" +to_string((centbins-i-1)*(90/centbins))+"-"+to_string((centbins-i)*(90/centbins))+".pdf").c_str());
	  deadmap[0][j][i]->GetYaxis()->SetLabelSize(0.025);
	  deadmap[0][j][i]->GetXaxis()->SetLabelSize(0.025);
	  //deadmap[1][j][i]->Divide(deadmap[0][j][i]);
	  deadmap[0][j][i]->GetYaxis()->SetTitle(ylabel.c_str());
	  deadmap[0][j][i]->GetXaxis()->SetTitle(xlabel.c_str());
	  deadmap[0][j][i]->GetZaxis()->SetLabelSize(0.025);
	  //deadmap[0][j][i]->GetZaxis()->SetRangeUser(0,2);
	  texts[2] = "Run 21615 " + to_string((centbins-i-1)*(90/centbins))+"-"+to_string((centbins-i)*(90/centbins))+"% centrality";
	  deadmap[0][j][i]->Draw("COLZ");
	  sphenixtext();
	  multitext(texts,ntext);
	  //c1->SaveAs(("/home/jocl/datatemp/plots/png/cent/sim_deadmap_"+cal[j]+"_" +to_string((centbins-i-1)*(90/centbins))+"-"+to_string((centbins-i)*(90/centbins))+".png").c_str());
	  deadmap[1][j][i]->Divide(deadmap[0][j][i]);
	  deadmap[1][j][i]->GetZaxis()->SetTitle((cal[j] +" <E^{tower}_{data}>/<E^{tower}_{HIJING}>").c_str());
	  deadmap[1][j][i]->GetZaxis()->SetRangeUser(0.5,1.5);
	  deadmap[1][j][i]->Draw("COLZ");
	  sphenixtext();
	  multitext(texts,ntext);
	  //c1->SaveAs(("/home/jocl/datatemp/plots/png/cent/div_deadmap_"+cal[j]+"_" +to_string((centbins-i-1)*(90/centbins))+"-"+to_string((centbins-i)*(90/centbins))+".png").c_str());

	}
    }

  for(int i=0; i<centbins; ++i)
    {
      cout << i << " " << ettotcent[0][i]->GetBinContent(centet[0][0][i]->GetNbinsX()+2) << endl;
    }
  return 0;
}


int plot()
{
  const int nfiles = 1;
    /*;
  string tags[nfiles] =
    {
      "unc",
      "cor"
    };
  string filenames[nfiles] =
    {
      //"savedhists_fracsim_100_fracdat_100_subtr_0_minE_0_scale_1.30_zcut_30_run_21615_ntc.root"
      //"savedhists_fracsim_1_fracdat_1_subtr_0_minE_0_scale_1.30_zcut_30_run_21615_ntc.root",
      //"savedhists_fracsim_1_fracdat_1_subtr_0_minE_5_scale_1.30_zcut_30_run_21615_ntc.root",
      "savedhists_fracsim_1_fracdat_1_subtr_0_minE_0_scale_1.30_zcut_10_run_21615_ntc.root"
      //"savedhists_fracsim_1_fracdat_1_subtr_0_minE_0_scale_1.30_zcut_10_run_21615_newest.root"
    };

  for(int i=0; i<nfiles; ++i)
    {
      int dummy = called_plot(filenames[i]);
    }
    */
  //called_plot("savedhists_fracsim_1_fracdat_1_subtr_0_minE_5_scale_1.30_zcut_10_run_21615_20231009_cor.root", "ttree","/home/jocl/datatemp/","/home/jocl/datatemp/plots/");
  //called_plot("savedhists_fracsim_1_fracdat_1_subtr_0_minE_5_scale_1.30_zcut_10_run_21615_20231009_unc.root", "ttree","/home/jocl/datatemp/","/home/jocl/datatemp/plots_unc/");
  called_plot("savedhists_fracsim_1_fracdat_1_subtr_0_minE_0_scale_1.00_zcut_10_run_21615_20231018_nopileup_cor.root");
  return 0;
}
