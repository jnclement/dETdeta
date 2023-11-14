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
  if(calnum!=4) leg = new TLegend(0.4,0.75,0.9,0.9);
  else leg = new TLegend(0.4,0.65,0.9,0.9);
  leg->SetNColumns(3);
  leg->SetTextSize(0.015);
  mainhist->SetMarkerColor(kBlack);
  if(calnum==0 || calnum == 3 || calnum==4) mainhist->Rebin(5);
  for(int i=0; i<centbins; ++i)
    {
      if(calnum == 0 || calnum == 3 || calnum==4) hists[i]->Rebin(5);
      hists[i]->SetMarkerColor(knames[i%6]+kcodes[i/6]);
      hists[i]->SetLineColor(knames[i%6]+kcodes[i/6]);
      hists[i]->SetFillColorAlpha(knames[i%6]+kcodes[i/6],0.2);
      leg->AddEntry(hists[i],(to_string((centbins-i-1)*centrange)+"-"+to_string((centbins-i)*centrange) + "% centrality").c_str(),"P");
    }
  leg->AddEntry(mainhist,("0-90% centrality"),"P");
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
  //cout << minval << " " << maxval << endl;
  //if(calnum == 0) maxval = 10;
  //else if(calnum == 1) maxval = 6;
  //else if(calnum == 2) maxval = 20;
  if(logy) mainhist->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  else mainhist->GetYaxis()->SetRangeUser(min(0,minval)-abs(min(0,minval))/10.,maxval+abs(maxval)/10.);
  float ranges[5] = {1500,150,400,2000,5000};
  mainhist->GetXaxis()->SetRangeUser(0,ranges[calnum]);
  //cout << mainhist->GetName() << endl;
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
  auto leg = new TLegend(0.1,0.025,0.3,0.1);
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
      if(simhist) minval = min(dathist->GetMinimum(),simhist->GetMinimum());
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
  else if(name=="zcent")
    {
      hp[1] = dathist->Fit("gaus","S");
      dathist->GetFunction("gaus")->SetLineColor(dcolor);
      dathist->Draw(options.c_str());
      dathist->GetFunction("gaus")->Draw("SAME");
    }
  /*
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
  */
  sphenixtext();
  //multitext(texts, ntext, 0.03, 0.11);
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
  else if(name=="zcent")
    {
      stringstream stmean, stsig;
      stmean << std::fixed << std::setprecision(2) << hp[1]->Parameter(1);
      stsig << std::fixed << std::setprecision(2) << hp[1]->Parameter(2);
      drawText(("#mu"+(string)(simhist?"_{data}":"")+"="+stmean.str()+", #sigma"+(simhist?"_{data}":"")+"="+stsig.str()).c_str(),0.9,0.91,1,kBlack,0.025);
    }
  if(name=="det_no_acc_cor")
    {
      drawText("No acceptance cut applied",0.3,0.025,0,kBlack,0.025);
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
      //multitext(texts2, ntext2);
      ca->SaveAs((dir+"pdf/"+"ratio_"+subdir+"ratio_"+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
      ca->SaveAs((dir+"png/"+"ratio_"+subdir+"ratio_"+name+"_"+cal+"_scale_"+params[0]+"_subtr_"+params[1]+"_mine_"+params[2]+"_zcut_"+params[3]+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
    }
}

void newplot(string options, TCanvas* ca, TH1* hist1, TH1* hist2, int logy, string cal, int run, float zcut, string xlabel, string ylabel,string datleg, string simleg, string raty, int percent0, int percent1, string name, string dir, string subdir, int centbins, int norm, string maintag = "", int rebin = 0)
{
  TH1D* dathist;
  TH1D* simhist;
  if(name == "et_event")
    {
      //gPad->SetLogx();
    }
  else
    {
      gPad->SetLogx(0);
    }
  //gPad->SetLogx();
  //if(name == "fullcor");
  if(rebin)
    {
      dathist = (TH1D*)hist1->Rebin(rebin,"datnewbins");
      simhist = (TH1D*)hist2->Rebin(rebin,"simnewbins");
    }
  else
    {
      dathist = (TH1D*)hist1;
      simhist = (TH1D*)hist2;
    }
  int centrange = 10;
  ca->cd();
  if(logy) gPad->SetLogy();
  else gPad->SetLogy(0);
  TH1D* datcopy;
  TH1D* simcopy;
  TH1D* rathist;
  string toprint;
  if(name != "truthpar_et")
    {
      toprint = "Run "+to_string(run)+" "+to_string(centrange*(centbins-percent1))+"-"+to_string(centrange*(centbins-percent0))+"% centrality";
    }
  else
    {
      toprint = "";
    }
  TLegend* leg = new TLegend(0.7,0.85,0.9,0.925);
  leg->SetTextSize(0.02);
  leg->SetFillColorAlpha(0,0);
  dathist->SetLineColor(kRed);
  dathist->SetMarkerColor(kRed);
  leg->AddEntry(dathist,datleg.c_str(),"P");
  dathist->GetXaxis()->SetTitle(xlabel.c_str());
  dathist->GetYaxis()->SetTitle(ylabel.c_str());
  dathist->GetYaxis()->SetLabelSize(0.025);
  dathist->GetXaxis()->SetLabelSize(0.025);
  if(norm) dathist->Scale(1./dathist->Integral());
  if(simhist)
    {
      if(norm) simhist->Scale(1./simhist->Integral());
      simhist->SetLineColor(kBlue);
      simhist->SetMarkerColor(kBlue);
      leg->AddEntry(simhist,simleg.c_str(),"P");
      float minval = 9999999999;
      if(logy)
	{
	  for(int i=0; i<dathist->GetNbinsX(); ++i)
	    {
	      if(dathist->GetBinContent(i) < minval && dathist->GetBinContent(i) > 0) minval = dathist->GetBinContent(i);
	    }
	  for(int i=0; i<simhist->GetNbinsX(); ++i)
	    {
	      if(simhist->GetBinContent(i) < minval && simhist->GetBinContent(i) > 0) minval = simhist->GetBinContent(i);
	    }
	}
      else minval = min(dathist->GetBinContent(dathist->GetMinimumBin()),simhist->GetBinContent(simhist->GetMinimumBin()));
      float maxval = max(dathist->GetBinContent(dathist->GetMaximumBin()),simhist->GetBinContent(simhist->GetMaximumBin()));
      float maxx = max(dathist->GetBinLowEdge(dathist->FindLastBinAbove(0)+1),simhist->GetBinLowEdge(simhist->FindLastBinAbove(0)+1));
      minval *= 0.9;
      maxval *= 1.1;
      maxx *= 1.1;
      dathist->GetYaxis()->SetRangeUser(minval,maxval);
      if(dathist->GetBinLowEdge(1) == 0) dathist->GetXaxis()->SetRangeUser(0,maxx);
      int lastx = min(dathist->FindLastBinAbove(0),simhist->FindLastBinAbove(0));
      if(dathist->GetBinLowEdge(1) == 0)
	{
	  datcopy = new TH1D("datcopy","",lastx,0,lastx);
	  simcopy = new TH1D("simcopy","",lastx,0,lastx);
	  rathist = new TH1D("rathist","",lastx,0,lastx);
	}
      else
	{
	  int nbin = dathist->GetNbinsX();
	  datcopy = new TH1D("datcopy","",nbin,dathist->GetBinLowEdge(1),dathist->GetBinLowEdge(nbin+1));
	  simcopy = new TH1D("simcopy","",nbin,simhist->GetBinLowEdge(1),simhist->GetBinLowEdge(nbin+1));
	  rathist = new TH1D("rathist","",nbin,dathist->GetBinLowEdge(1),dathist->GetBinLowEdge(nbin+1));
	}
      for(int i=1; i<lastx+1; ++i)
	{
	  datcopy->SetBinContent(i,dathist->GetBinContent(i));
	  datcopy->SetBinError(i,dathist->GetBinError(i));
	  simcopy->SetBinContent(i,simhist->GetBinContent(i));
	  simcopy->SetBinError(i,simhist->GetBinError(i));
	}
      rathist->Divide(datcopy,simcopy);
      rathist->SetLineColor(kRed);
      rathist->SetMarkerColor(kRed);
      rathist->GetYaxis()->SetLabelSize(0.025);
      rathist->GetXaxis()->SetLabelSize(0.025);
      rathist->GetXaxis()->SetTitle(xlabel.c_str());
      rathist->GetYaxis()->SetTitle(raty.c_str());
      cout << minval << " " << maxval << " " << dathist->GetName()  << " " << simhist->GetName() << endl << endl;
      
      dathist->Draw();
      simhist->Draw("same");
      leg->Draw();
      ca->SaveAs((dir+"png/"+subdir+name+"_"+maintag+"_"+cal+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
      ca->SaveAs((dir+"pdf/"+subdir+name+"_"+maintag+"_"+cal+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());

      if(cal=="EMCal" && name=="et_event")
	{
	  dathist->GetXaxis()->SetRangeUser(-100,100);
	  dathist->Draw();
	  simhist->Draw("same");
	  leg->Draw();
	}
      else if(cal=="OHCal" && name=="et_event")
	{
	  dathist->GetXaxis()->SetRangeUser(-40,40);
	  dathist->Draw();
	  simhist->Draw("same");
	  leg->Draw();
	}
      else if(cal=="IHCal" && name=="et_event")
	{
	  dathist->GetXaxis()->SetRangeUser(-12,12);
	  dathist->Draw();
	  simhist->Draw("same");
	  leg->Draw();
	}
      if((cal=="EMCal" || cal=="OHCal" || cal=="IHCal") && name=="et_event")
	{
      drawText(toprint.c_str(),0.2,0.96,0,kBlack,0.025);
      sphenixtext();
      
      ca->SaveAs((dir+"png/"+subdir+name+"_zoom_"+maintag+"_"+cal+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
      ca->SaveAs((dir+"pdf/"+subdir+name+"_zoom_"+maintag+"_"+cal+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
	}
      gPad->SetLogy(0);
      if(rathist->GetMaximum() > 2) rathist->GetYaxis()->SetRangeUser(0,2);
      rathist->Draw();
      drawText(toprint.c_str(),0.2,0.96,0,kBlack,0.025);
      sphenixtext();
      ca->SaveAs((dir+"png/ratio_"+subdir+"ratio_"+name+"_"+maintag+"_"+cal+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
      ca->SaveAs((dir+"pdf/ratio_"+subdir+"ratio_"+name+"_"+maintag+"_"+cal+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
      delete datcopy;
      delete simcopy;
      delete rathist;
    }
  else
    {
      dathist->Draw();
      drawText(toprint.c_str(),0.2,0.96,0,kBlack,0.025);
      sphenixtext();
      ca->SaveAs((dir+"png/"+subdir+name+"_"+maintag+"_"+cal+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".png").c_str());
      ca->SaveAs((dir+"pdf/"+subdir+name+"_"+maintag+"_"+cal+"_cent_"+ to_string(centrange*(centbins-percent1)) +"-"+ to_string(centrange*(centbins-percent0))+".pdf").c_str());
    }
}

int setrebin(int j)
{
  int rebin = 0;
  if(j==0) rebin = 10;
  else if(j==2) rebin = 5;
  return rebin;
}


int called_plot(string histfilename = "savedhists_fracsim_1_fracdat_1_subtr_0_minE_0_scale_1.00_zcut_30_run_21615_20231018_nopileup_cor.root", string maintag = "", string treename = "ttree", string datdir = "/home/jocl/datatemp/", string plotdir = "/home/jocl/datatemp/plots/")
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
  TH1D* truthpar_total_ET;
  float sub;
  float scale[2];
  int frac[2];
  float mine;
  float zcut;
  int run;
  truthpar_total_ET = (TH1D*)histfile->Get("truthpar_total_ET");
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
  //cout << centet[1][0][17]->GetEntries() << endl;
  //centet[1][0][17]->Draw();
  //c2->SaveAs("test.png");
  //cout << "Finished getting hists" << endl;
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
	  //cout << ("fullcor_"+to_string(i)+to_string(j)).c_str() << endl;
	  fullcor[i][j] = (TH1D*)histfile->Get(("fullcor_"+to_string(i)+to_string(j)).c_str());
	  //cout << fullcor[i][j] << endl;
	  meandiffnoavg[i]->SetBinContent(centbins-j,centet[1][i][j]->GetMean()-centet[0][i][j]->GetMean());
	  meandiffnoavg[i]->SetBinError(centbins-j,0);
	}
    }
  //cout << "Finished filling new hists" << endl;
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
  //cout << "Starting to plot" << endl;
  string options = "hist";
  xlabel = "E_{T, event}^{all} [GeV]";
  string ylabel = "Counts";
  //centoverlayplot(options, c1, sumev[0], ettotcent[0], 1, "all", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "cent_overlay_sim", plotdir, "all/", centbins, 0, 3);
  //centoverlayplot(options, c1, sumev[1], ettotcent[1], 1, "all", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "cent_overlay_dat", plotdir, "all/", centbins, 1, 3);
  xlabel = "Truth Particle N";
  //centoverlayplot(options, c1, truthparnhist, truthparncent, 1, "truth particles", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "n_cent_overlay",plotdir, "all/", centbins, 0, 4);
  xlabel = "MBD Z Vertex [cm]";
  ylabel = "Normalized Counts";
  options = "";
  plotsimdat(options, c1, zhist[1], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "zvtx", plotdir,"all/", centbins, 1);
  xlabel = "MBD charge sum [??]";
  //plotsimdat(options, c1, mbh[1], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "mboverlay_dat", plotdir,"all/", centbins, 1);
  //plotsimdat(options, c1, mbh[0], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "mboverlay_sim", plotdir,"all/", centbins, 1);
  //plotsimdat(options, c1, mbh[1], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, 0, 20, "mbdat", plotdir);
  //plotsimdat(options, c1, mbh[0], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, 0, 20, "mbsim", plotdir);
  xlabel = "E_{T, event calorimeter sum}";
  //plotsimdat(options, c1, sumev[1], sumev[0], 1, "total", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_total", plotdir,"all/", centbins, 1);
  xlabel = "E_{T, stacked calorimeter towers}";
  //plotsimdat(options, c1, sumtw[1], sumtw[0], 1, "total", scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_sumtow", plotdir,"all/", centbins, 1);
  //cout << "Plotted non-loops" << endl;
  string datleg = "Data";
  string simleg = "HIJING";
  string raty = "";
  int norm = 0;
  for(int j=0; j<3; ++j)
    {
      norm = 0;
      xlabel = "#eta";
      ylabel = cal[j]+" Corrected dE_{T}/d#eta [GeV]";
      options="p";
      //plotsimdat(options, c1, fullcor[j][centbins-1], NULL, 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, centbins-1, centbins, "fullcor", plotdir, "cent/", centbins, 0);
      int rebin = 0;
      newplot(options, c1, fullcor[j][centbins-1], NULL, 0, cal[j], run, zcut, xlabel, ylabel, datleg, simleg, raty, centbins-1, centbins, "fullcor", plotdir, "cent/", centbins, norm, maintag);
      ylabel = "HIJING "+cal[j]+" dE_{T}/d#eta [GeV]";
      //if(j==0) plotsimdat(options, c1, dETcent[0][j][centbins-1],NULL, 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, centbins-1, centbins, "det_no_acc_cor",plotdir,"cent/",centbins,0);
      if(j==0) newplot(options, c1, dETcent[0][j][centbins-1], NULL, 0, cal[j], run, zcut, xlabel, ylabel, datleg, simleg, raty, centbins-1, centbins, "det_cent_sim", plotdir, "cent/", centbins, norm,maintag);
      ylabel = "Truth Particle dE_{T}/d#eta [GeV]";
      if(j==0) plotsimdat(options, c1, truthpar_et[centbins-1], NULL, 0, "", scale[0], sub, run, mine, zcut, xlabel, ylabel, centbins-1,centbins, "truthpar_et", plotdir, "cent/", centbins, 0);
      if(j==0) newplot(options, c1, truthpar_et[centbins-1], NULL, 0, "truth", run, zcut, xlabel, ylabel, datleg, simleg, raty, centbins-1, centbins, "truthpar_et", plotdir, "cent/", centbins, norm,maintag);
      raty = "(Reco dE_{T}/d#eta)/(Truth dE_{T}/d#eta)";
      //plotsimdat(options, c1, dETcent[0][j][centbins-1],truthpar_et[centbins-1],0,cal[j],scale[0],sub,run,mine,zcut,xlabel,ylabel,centbins-1,centbins,"truth_reco_et",plotdir,"cent/",centbins,0);
      ylabel = cal[j]+" dE_{T}/d#eta [GeV]";
      datleg = "Reco HIJING";
      simleg = "Truth HIJING";
      newplot(options, c1, dETcent[0][j][centbins-1],truthpar_et[centbins-1], 0, cal[j], run, zcut, xlabel, ylabel, datleg, simleg, raty, centbins-1, centbins, "recotruth_detdeta", plotdir, "cent/", centbins, norm, maintag);
      datleg = "Data";
      raty = "(Data dE_{T}/d#eta)/(Truth dE_{T}/d#eta)";
      //plotsimdat(options, c1, dETcent[1][j][centbins-1],truthpar_et[centbins-1],0,cal[j],scale[0],sub,run,mine,zcut,xlabel,ylabel,centbins-1,centbins,"truth_data_et",plotdir,"cent/",centbins,0);
      newplot(options, c1, dETcent[1][j][centbins-1],truthpar_et[centbins-1], 0, cal[j], run, zcut, xlabel, ylabel, datleg, simleg, raty, centbins-1, centbins, "datatruth_detdeta", plotdir, "cent/", centbins, norm, maintag);
      ylabel = "Raw dE_{T}/d#eta in "+cal[j];
      raty = cal[j]+" (Data dE_{T}/d#eta)/(Reco dE_{T}/d#eta)";
      simleg = "Reco HIJING";
      newplot(options, c1, dETcent[1][j][centbins-1],dETcent[0][j][centbins-1], 0, cal[j], run, zcut, xlabel, ylabel, datleg, simleg, raty, centbins-1, centbins, "datareco_detdeta", plotdir, "cent/", centbins, norm, maintag);
      options = "hist";
      xlabel = "#Sigma #it{E}_{T}^{"+cal[j]+"} [GeV]";
      ylabel = "Counts";
      //centoverlayplot(options, c1, ET[1][j], centet[1][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "cent_overlay_dat", plotdir, "all/", centbins, 1, j);
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
      xlabel = "#Sigma E_{T} [GeV]";
      //plotsimdat(options, c1, ET[1][j], ET[0][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_event", plotdir, "all/", centbins, 1);
      rebin = setrebin(j);
      datleg = "Data "+cal[j];
      simleg = "HIJING "+cal[j];
      raty = cal[j]+" Ratio Data/HIJING";
      norm = 1;
      newplot(options, c1, ET[1][j],ET[0][j], 1, cal[j], run, zcut, xlabel, ylabel, datleg, simleg, raty, 0, centbins, "et_event", plotdir, "cent/", centbins, norm, maintag,rebin);
      cout << "test" << endl;
	//plotsimdat(options, c1, truthpar_total_ET, ET[0][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_event_truthpar_sim", plotdir, "all/", centbins, 1);
      datleg = "Truth";
      raty = cal[j]+" Ratio Truth/HIJING";
      if(maintag == "wzs" || maintag == "nzs") newplot(options, c1, truthpar_total_ET, ET[0][j], 1, cal[j], run, zcut, xlabel, ylabel, datleg, simleg, raty, 0, centbins, "et_event_truthpar_sim", plotdir, "cent/", centbins, norm, maintag,rebin);
	//plotsimdat(options, c1, truthpar_total_ET, ET[1][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_event_truthpar_data", plotdir, "all/", centbins, 1);
      simleg = "Data "+cal[j];
      raty = cal[j]+" Ratio Truth/Data";
      if(maintag == "wzs" || maintag == "nzs") newplot(options, c1, truthpar_total_ET, ET[1][j], 1, cal[j], run, zcut, xlabel, ylabel, datleg, simleg, raty, 0, centbins, "et_event_truthpar_data", plotdir, "cent/", centbins, norm, maintag,rebin);
      xlabel = "E_{T," + cal[j] +" tower} [GeV]";
      //plotsimdat(options, c1, TW[1][j], TW[0][j], 1, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0, centbins, "et_tower", plotdir, "all/", centbins, 1);
      if(j==0) xlabel = "floor{(#eta bin " +cal[j]+")/4}";
      else xlabel = "#eta bin " + cal[j];
      ylabel = cal[j]+" dE_{T}/d#eta [GeV]";
      options = "p";
      //plotsimdat(options, c1, dET[1][j],dET[0][j],0,cal[j], scale[0],sub,run,mine,zcut,xlabel,ylabel,0,centbins,"det",plotdir,"all/",centbins, 0);

      
      //multiplot(options, c1, dETcent[0][j], 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "detcent_sim", plotdir, "all/", centbins, 0, j);
      //multiplot(options, c1, dETcent[1][j], 0, cal[j], scale[1], sub, run, mine, zcut, xlabel, ylabel, 0,centbins, "detcent_data", plotdir, "all/", centbins, 1, j);
      //cout << "test" << endl;
      for(int k=0; k<centbins; ++k)
	{
	  if(j==0)
	    {
	      xlabel = "Z-vertex [cm]";
	      ylabel = "Normalized Counts";
	      options = "p";
	      //plotsimdat(options, c1, zcent[1][k], NULL, 1, "MBD", scale[0], sub, run, mine, zcut, xlabel, ylabel, k, k+1, "zcent", plotdir, "cent/", centbins, 1);
	    }
	  //outputon(gSystem);
	  //cout << j << " " << k << " " << centet[1][j][k]->GetMean() << " " << dETcent[1][j][k]->Integral() << " " <<centet[0][j][k]->GetMean() << " " << dETcent[0][j][k]->Integral() << endl;
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
	  //plotsimdat(options, c1, dETcent[1][j][k], dETcent[0][j][k], 0, cal[j], scale[0], sub, run, mine, zcut, xlabel, ylabel,k,k+1,"detcent",plotdir,"cent/",centbins, 0);
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
  gPad->SetLogx(0);
  for(int i=centbins-1; i<centbins; ++i)
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
  /*
  for(int i=0; i<centbins; ++i)
    {
      cout << i << " " << ettotcent[0][i]->GetBinContent(centet[0][0][i]->GetNbinsX()+2) << endl;
    }
  */
  float test = 0;
  TH1D* star05 = new TH1D("star05","",1,0.49,0.51);
  TH1D* phnx05 = new TH1D("phnx05","",1,-0.01,0.01);
  TH1D* phnx510 = new TH1D("phnx510","",1,-0.01,0.01);
  TH1D* star510 = new TH1D("star510","",1,0.49,0.51);
  TH1D* spnx010[3];
  spnx010[0] = new TH1D("spnx05em","",1,-1,1);
  spnx010[1] = new TH1D("spnx05ih","",1,-1,1);
  spnx010[2] = new TH1D("spnx05oh","",1,-1,1);
  star05->SetBinContent(1,3.51);
  star05->SetBinError(1,0.19);
  star510->SetBinContent(1,3.43);
  star510->SetBinError(1,0.2);
  phnx05->SetBinContent(1,3.41);
  phnx05->SetBinError(1,0.2);
  phnx510->SetBinContent(1,3.29);
  phnx510->SetBinError(1,0.2);
  for(int i=0; i<3; ++i)
    {
      spnx010[i]->SetBinContent(1,fullcor[i][centbins-1]->Integral("width")/(2*175));
    }
  star05->SetMarkerColor(kRed-4);
  star510->SetMarkerColor(kRed+1);
  phnx05->SetMarkerColor(kGreen-4);
  phnx510->SetMarkerColor(kGreen+2);
  spnx010[0]->SetMarkerColor(kBlack);
  spnx010[1]->SetMarkerColor(kBlue+1);
  spnx010[2]->SetMarkerColor(kBlue-7);
  star05->SetLineColor(kRed-4);
  star510->SetLineColor(kRed+1);
  phnx05->SetLineColor(kGreen-4);
  phnx510->SetLineColor(kGreen+2);
  spnx010[0]->SetLineColor(kBlack);
  spnx010[1]->SetLineColor(kBlue+1);
  spnx010[2]->SetLineColor(kBlue-7);
  star05->SetMarkerSize(2);
  star510->SetMarkerSize(2);
  phnx05->SetMarkerSize(2);
  phnx510->SetMarkerSize(2);
  spnx010[0]->SetMarkerSize(2);
  spnx010[1]->SetMarkerSize(2);
  spnx010[2]->SetMarkerSize(2);
  spnx010[0]->GetYaxis()->SetTitle("dE_{T}/d#eta/(0.5N_{part}) [GeV]");
  spnx010[0]->GetXaxis()->SetTitle("#eta Range");
  spnx010[0]->GetXaxis()->SetLabelSize(0.03);
  spnx010[0]->GetYaxis()->SetLabelSize(0.03);
  spnx010[0]->GetYaxis()->SetRangeUser(3,4.5);
  TLegend* leg = new TLegend(0,0.01,0.7,0.11);
  leg->SetNColumns(2);
  leg->AddEntry(star05,"STAR 0-5% Centrality","P");
  leg->AddEntry(phnx05,"PHENIX 0-5% Centrality","P");
  leg->AddEntry(star510,"STAR 5-10% Centrality","P");
  leg->AddEntry(phnx510,"PHENIX 5-10% Centrality","P");
  leg->AddEntry(spnx010[1],"sPHENIX IHCal 0-10% Centrality","P");
  leg->AddEntry(spnx010[0],"sPHENIX EMCal 0-10% Centrality","P");
  leg->AddEntry(spnx010[2],"sPHENIX OHCal 0-10% Centrality","P");
  TCanvas* c3 = new TCanvas("","",1000,1000);
  spnx010[0]->Draw("P");
  spnx010[1]->Draw("SAME P");
  spnx010[2]->Draw("SAME P");
  phnx510->Draw("SAME P");
  phnx05->Draw("SAME P");
  star510->Draw("SAME P");
  star05->Draw("SAME P");
  leg->SetTextSize(0.02);
  leg->Draw();
  drawText("Run 21615 0-5% centrality",0.2,0.96,0,kBlack,0.025);
  sphenixtext();
  c3->SaveAs("/home/jocl/datatemp/plots/pdf/all/detcomp.pdf");
  c3->SaveAs("/home/jocl/datatemp/plots/png/all/detcomp.png");

  TLine* starline = new TLine(0,620,1,620);
  TLine* phnxline = new TLine(-0.35,597,0.35,597);
  TBox* starbox = new TBox(0,620-33,1,620+33);
  TBox* phnxbox = new TBox(-0.35,597-35,0.35,597+35);
  starline->SetLineColor(kGreen);
  starbox->SetFillColorAlpha(kGreen,0.3);
  phnxline->SetLineColor(kBlue);
  phnxbox->SetFillColorAlpha(kBlue,0.3);
  fullcor[0][centbins-1]->GetYaxis()->SetTitle("Fully Corrected dE_{T}/d#eta [GeV]");
  fullcor[0][centbins-1]->GetYaxis()->SetRangeUser(500,900);
  fullcor[0][centbins-1]->SetMarkerColor(kRed);
  fullcor[1][centbins-1]->SetMarkerColor(kRed-9);
  fullcor[2][centbins-1]->SetMarkerColor(kRed+3);
  fullcor[0][centbins-1]->SetLineColor(kRed);
  fullcor[1][centbins-1]->SetLineColor(kRed-9);
  fullcor[2][centbins-1]->SetLineColor(kRed+3);
  TLegend* tleg = new TLegend(0.5,0.8,0.9,0.92);
  tleg->SetNColumns(2);
  tleg->SetTextSize(0.02);
  tleg->AddEntry(fullcor[0][centbins-1],"sPHENIX EMCal","P");
  tleg->AddEntry(fullcor[1][centbins-1],"sPHENIX IHCal","P");
  tleg->AddEntry(fullcor[2][centbins-1],"sPHENIX OHCal","P");
  tleg->AddEntry(phnxline,"PHENIX Measurement","L");
  tleg->AddEntry(starline,"STAR Measurement","L");
  tleg->SetFillColorAlpha(0,0);
  fullcor[0][centbins-1]->Draw();
  fullcor[1][centbins-1]->Draw("SAME P");
  fullcor[2][centbins-1]->Draw("SAME P");
  phnxline->Draw("SAME");
  starline->Draw("SAME");
  starbox->Draw("SAME");
  phnxbox->Draw("SAME");
  tleg->Draw();
  drawText("Run 21615 0-5% centrality",0.2,0.96,0,kBlack,0.025);
  sphenixtext();
  c3->SaveAs("/home/jocl/datatemp/plots/pdf/all/fullcor_all.pdf");
  c3->SaveAs("/home/jocl/datatemp/plots/png/all/fullcor_all.png");
  for(int i=0; i<3; i++)
    {
      cout << cal[i] << " mean corrected dET/deta/(0.5*npart): " << fullcor[i][centbins-1]->Integral("width")/(2*175) << endl;
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
  gErrorIgnoreLevel = kWarning;
  //called_plot("savedhists_fracsim_1_fracdat_1_subtr_0_minE_-10000_scale_1.00_zcut_30_run_21615_20231106_nopileup_wzs_cor.root","wzs");
  //called_plot("savedhists_fracsim_1_fracdat_1_subtr_0_minE_-10000_scale_1.00_zcut_30_run_21615_20231018_nopileup_cor.root","new18");
  //called_plot("savedhists_fracsim_1_fracdat_1_subtr_0_minE_-10000_scale_1.00_zcut_30_run_21615_20231106_nopileup_nzs_cor.root","nzs");
  //called_plot("savedhists_fracsim_1_fracdat_1_subtr_0_minE_-10000_scale_1.00_zcut_30_run_21615_20231026_nopileup_cor.root","new26");
  called_plot("savedhists_fracsim_1_fracdat_1_subtr_0_minE_-10000_scale_1.00_zcut_30_run_21615_20231108_nopileup_wzs_nodanvtx_cor.root","nodan");
  return 0;
}
