#include "TApplication.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TRandom2.h"
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
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <TMinuit.h>
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
const int centbins = 9;
void plotcentet(TCanvas* ca, TH1* hist1, TH1* hist2, int percent, int percent2, int logy, string cal, string sc, string sub, string run)
{
  ca->cd();
  if(logy) gPad->SetLogy();
  gPad->SetTicks(1);
  hist1->Scale(1./hist1->Integral());
  hist2->Scale(1./hist2->Integral());
  hist1->SetLineColor(kRed);
  hist2->SetLineColor(kBlue);
  float maxval = max(hist1->GetMaximum(),hist2->GetMaximum());
  float minval = min(hist1->GetBinContent(hist1->FindLastBinAbove(0,1)),hist2->GetBinContent(hist2->FindLastBinAbove(0,1)));
  hist1->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  //hist2->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  hist1->GetXaxis()->SetTitle(("E_{T,total "+cal+"} [GeV]").c_str());
  hist1->GetYaxis()->SetTitle("Counts");
  hist1->GetYaxis()->SetLabelSize(0.025);
  hist1->GetXaxis()->SetLabelSize(0.025);
  hist1->Draw();
  hist2->Draw("SAME");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText("Red: data",0.85,0.76,1,kBlack,0.025);
  drawText("Blue: sim",0.85,0.73,1,kBlack,0.025);
  drawText(("Sim scaled "+sc).c_str(),0.85,0.79,1,kBlack,0.025);
  drawText("Area normed to 1",0.85,0.82,1,kBlack,0.025);
  drawText("|z|<10cm",0.85,0.7,1,kBlack,0.025);
  drawText(("Run " + run + " " + to_string(10*(centbins-percent)) +"-"+to_string(10*(centbins-percent2))+"% centrality").c_str(),0.85,0.67,1,kBlack,0.025);
  drawText((sub+" MeV subtracted from each tower").c_str(),0.85,0.85,1,kBlack,0.025);
  ca->SaveAs(("dcentet" + cal + "_" + sc + "_cent" + to_string(10*(centbins-percent)) +"-"+to_string(10*(centbins-percent2)) + ".png").c_str());
  hist1->Divide(hist2);
  hist1->GetYaxis()->SetTitle("Data/Sim");
  hist1->GetYaxis()->SetRangeUser(0.01,10);
  hist1->Draw();
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText(("Run " + run + " " + to_string(10*(centbins-percent)) +"-"+to_string(10*(centbins-percent2))+"% centrality").c_str(),0.85,0.85,1,kBlack,0.025);
  ca->SaveAs(("ratio_dcentet" + cal + "_" + sc + "_cent" + to_string(10*(centbins-percent)) +"-"+to_string(10*(centbins-percent2)) + ".png").c_str());
  
}

void plotcenttow(TCanvas* ca, TH1* hist1, TH1* hist2, int percent, int percent2, int logy, string cal, string sc, string sub, string run)
{
  ca->cd();
  if(logy) gPad->SetLogy();
  gPad->SetTicks(1);
  hist1->Scale(1./hist1->Integral());
  hist2->Scale(1./hist2->Integral());
  hist1->SetLineColor(kRed);
  hist2->SetLineColor(kBlue);
  float maxval = max(hist1->GetMaximum(),hist2->GetMaximum());
  float minval = min(hist1->GetBinContent(hist1->FindLastBinAbove(0,1)),hist2->GetBinContent(hist2->FindLastBinAbove(0,1)));
  hist1->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  //hist2->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  hist1->GetXaxis()->SetTitle(("E_{T,tower "+cal+"} [GeV]").c_str());
  hist1->GetYaxis()->SetTitle("Counts");
  hist1->GetYaxis()->SetLabelSize(0.025);
  hist1->GetXaxis()->SetLabelSize(0.025);
  hist1->Draw();
  hist2->Draw("SAME");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText("Red: data",0.85,0.76,1,kBlack,0.025);
  drawText("Blue: sim",0.85,0.73,1,kBlack,0.025);
  drawText(("Sim scaled "+sc).c_str(),0.85,0.79,1,kBlack,0.025);
  drawText("Area normed to 1",0.85,0.82,1,kBlack,0.025);
  drawText("|z|<10cm",0.85,0.7,1,kBlack,0.025);
  drawText(("Run " + run+ " " + to_string(10*(centbins-percent)) +"-"+to_string(10*(centbins-percent2))+"% centrality").c_str(),0.85,0.67,1,kBlack,0.025);
  drawText((sub+" MeV subtracted from each tower").c_str(),0.85,0.85,1,kBlack,0.025);
  ca->SaveAs(("dcenttow" + cal + "_" + sc + "_cent" + to_string(10*(centbins-percent)) +"-"+to_string(10*(centbins-percent2)) + ".png").c_str());
  hist1->Divide(hist2);
  hist1->GetYaxis()->SetTitle("Data/Sim");
  hist1->GetYaxis()->SetRangeUser(0.01,10);
  hist1->Draw();
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText(("Run " + run + " " + to_string(10*(centbins-percent)) +"-"+to_string(10*(centbins-percent2))+"% centrality").c_str(),0.85,0.85,1,kBlack,0.025);
  ca->SaveAs(("ratio_dcenttow" + cal + "_" + sc + "_cent" + to_string(10*(centbins-percent)) +"-"+to_string(10*(centbins-percent2)) + ".png").c_str());
}
int build_events()
{
  mbd_init();
  //cout << gaincorr[0] << endl;
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  float mbenrgy[25000], ohcalen[25000], ihcalen[25000], emcalen[25000];
  int emcalet[25000], emcalph[25000];
  int ihcalet[25000], ihcalph[25000];
  int ohcalet[25000], ohcalph[25000];
  int   mbdtype[25000], mbdside[25000], mbdchan[25000];
  float smbenrgy[25000], sohcalen[25000], sihcalen[25000], semcalen[25000];
  int semcalet[25000], semcalph[25000];
  int sihcalet[25000], sihcalph[25000];
  int sohcalet[25000], sohcalph[25000];
  float cents[11] = {0};
  float semtowercomb[64][24];
  float demtowercomb[64][24];
  int sectorem, sectoroh, sectorih, sectormb;
  int ssectorem, ssectoroh, ssectorih, ssectormb;
  TFile* hottowers = TFile::Open("/home/jocl/datatemp/hot_towers_21518_1.root");
  TTree* hottree = hottowers->Get<TTree>("T_hot_tower");
  TFile* file = TFile::Open("/home/jocl/datatemp/merged_dEdeta_71.root");
  TTree* tree = file->Get<TTree>("ttree");
  TFile* simf = TFile::Open("/home/jocl/datatemp/merged_dEdeta_250.root");
  TTree* simt = simf->Get<TTree>("ttree");
  TH1D* hist = new TH1D("hist","",500,0,500);
  TH1D* dmbh = new TH1D("dmbh","",3000,0,3000);
  double dcent[centbins+1] = {0};
  dcent[centbins] = 999999;
  TH1D* centclass[centbins][3];
  TH1D* dcenttow[3][centbins];
  TH1D* scenttow[3][centbins];
  TH1D* dcentet[3][centbins];
  TH1D* scentet[3][centbins];
  float et_em_range[centbins] = {100,150,275,350,400,600,800,1200,1750};
  float et_oh_range[centbins] = {35,50,80,100,140,175,225,300,400};
  float et_ih_range[centbins] = {10,15,25,35,50,75,100,120,150};
  for(int i=0; i<centbins; ++i)
    {
      dcenttow[0][i] = new TH1D(("dcenttowem" + to_string(i)).c_str(),"",50,0,10);
      dcenttow[1][i] = new TH1D(("dcenttowih" + to_string(i)).c_str(),"",50,0,2);
      dcenttow[2][i] = new TH1D(("dcenttowoh" + to_string(i)).c_str(),"",50,0,10);

      scenttow[0][i] = new TH1D(("scenttowem" + to_string(i)).c_str(),"",50,0,10);
      scenttow[1][i] = new TH1D(("scenttowih" + to_string(i)).c_str(),"",50,0,2);
      scenttow[2][i] = new TH1D(("scenttowoh" + to_string(i)).c_str(),"",50,0,10);


      dcentet[0][i] = new TH1D(("dcentetem" + to_string(i)).c_str(),"",100,0,et_em_range[i]);
      dcentet[1][i] = new TH1D(("dcentetih" + to_string(i)).c_str(),"",100,0,et_ih_range[i]);
      dcentet[2][i] = new TH1D(("dcentetoh" + to_string(i)).c_str(),"",100,0,et_oh_range[i]);

      scentet[0][i] = new TH1D(("scentetem" + to_string(i)).c_str(),"",100,0,et_em_range[i]);
      scentet[1][i] = new TH1D(("scentetih" + to_string(i)).c_str(),"",100,0,et_ih_range[i]);
      scentet[2][i] = new TH1D(("scentetoh" + to_string(i)).c_str(),"",100,0,et_oh_range[i]);
    }

  TH1D* mc20ets[3];
  TH1D* lc30ets[3];
  TH1D* mc20tws[3];
  TH1D* lc30tws[3];
  TH1D* mc20etd[3];
  TH1D* lc30etd[3];
  TH1D* mc20twd[3];
  TH1D* lc30twd[3];
  mc20ets[0] = new TH1D("mc20ets0","",100,300,1100);
  mc20ets[1] = new TH1D("mc20ets1","",100,25,100);
  mc20ets[2] = new TH1D("mc20ets2","",100,50,300);
  mc20tws[0] = new TH1D("mc20tws0","",100,0,5);
  mc20tws[1] = new TH1D("mc20tws1","",100,0,1);
  mc20tws[2] = new TH1D("mc20tws2","",100,0,5);
  mc20etd[0] = new TH1D("mc20etd0","",100,300,1100);
  mc20etd[1] = new TH1D("mc20etd1","",100,25,100);
  mc20etd[2] = new TH1D("mc20etd2","",100,50,300);
  mc20twd[0] = new TH1D("mc20twd0","",100,0,5);
  mc20twd[1] = new TH1D("mc20twd1","",100,0,1);
  mc20twd[2] = new TH1D("mc20twd2","",100,0,5);
  lc30ets[0] = new TH1D("lc30ets0","",100,0,100);
  lc30ets[1] = new TH1D("lc30ets1","",100,0,10);
  lc30ets[2] = new TH1D("lc30ets2","",100,0,30);
  lc30tws[0] = new TH1D("lc30tws0","",100,0,3);
  lc30tws[1] = new TH1D("lc30tws1","",100,0,1);
  lc30tws[2] = new TH1D("lc30tws2","",100,0,5);
  lc30etd[0] = new TH1D("lc30etd0","",100,0,100);
  lc30etd[1] = new TH1D("lc30etd1","",100,0,10);
  lc30etd[2] = new TH1D("lc30etd2","",100,0,30);
  lc30twd[0] = new TH1D("lc30twd0","",100,0,3);
  lc30twd[1] = new TH1D("lc30twd1","",100,0,1);
  lc30twd[2] = new TH1D("lc30twd2","",100,0,5);
  
  TH1D* etdata = new TH1D("etdata","",50,0,1200);
  TH1D* etsim  = new TH1D("etsim","",50,0,1200);
  TH1D* sedist = new TH1D("sedist","",50,0,10);
  TH1D* dedist = new TH1D("dedist","",50,0,10);
  TH1D* ietdata = new TH1D("ietdata","",50,0,120);
  TH1D* ietsim  = new TH1D("ietsim","",50,0,120);
  TH1D* isedist = new TH1D("isedist","",50,0,2);
  TH1D* idedist = new TH1D("idedist","",50,0,2);
  TH1D* oetdata = new TH1D("oetdata","",50,0,350);
  TH1D* oetsim  = new TH1D("oetsim","",50,0,350);
  TH1D* osedist = new TH1D("osedist","",50,0,10);
  TH1D* odedist = new TH1D("odedist","",50,0,10);
  TH1D* mthist = new TH1D("mthist","",100,-50,50);
  TH1D* rhist = new TH1D("rhist","",50,0,1200);
  TH1D* rthist = new TH1D("rthist","",50,0,10);
  TH1D* rihist = new TH1D("rihist","",50,0,120);
  TH1D* rithist = new TH1D("rithist","",50,0,2);
  TH1D* rohist = new TH1D("rohist","",50,0,350);
  TH1D* rothist = new TH1D("rothist","",50,0,10);
  TH1D* ssumhist = new TH1D("ssumhist","",50,0,1500);
  TH1D* dsumhist = new TH1D("dsumhist","",50,0,1500);
  TH1D* sumrat = new TH1D("sumrat","",50,0,1500);
  TH1D* ssumtow = new TH1D("ssumtow","",50,0,15);
  TH1D* dsumtow = new TH1D("dsumtow","",50,0,15);
  TH1D* rsumtow = new TH1D("rsumtow","",50,0,15);
  for(int i=0; i<centbins; ++i)
    {
      centclass[i][0] = new TH1D(("centhiste" + to_string(i)).c_str(),"",96,-1.1,1.1);
      centclass[i][1] = new TH1D(("centhisti" + to_string(i)).c_str(),"",24,-1.1,1.1);
      centclass[i][2] = new TH1D(("centhisto" + to_string(i)).c_str(),"",24,-1.1,1.1);
    }
  float detacente[centbins][96] = {0};
  float detacento[centbins][24] = {0};
  float detacenti[centbins][24] = {0};
  int nem[centbins][96] = {0};
  int nih[centbins][24] = {0};
  int noh[centbins][24] = {0};
  int tpn = 0;
  float subtr = 0;//0.018;
  float scale = 1.3;
  float mine = 0;//0.005;
  float tpe[25000] = {0};
  float tpet[25000] = {0};
  float tpph[25000] = {0};
  int tpem[25000] = {0};
  int npart = 0;
  int ncoll = 0;
  float bimp = 0;
  TH2D* ectee = new TH2D("ectee","",100,0,2000,100,0,3000);
  TH2D* ectea = new TH2D("ectea","",100,0,2000,100,0,3000);
  tree->SetBranchAddress("mbenrgy",mbenrgy);
  tree->SetBranchAddress("emcalen",emcalen);
  tree->SetBranchAddress("emcalet",emcalet);
  tree->SetBranchAddress("emcalph",emcalph);
  tree->SetBranchAddress("ihcalen",ihcalen);
  tree->SetBranchAddress("ihcalet",ihcalet);
  tree->SetBranchAddress("ihcalph",ihcalph);
  tree->SetBranchAddress("ohcalen",ohcalen);
  tree->SetBranchAddress("ohcalet",ohcalet);
  tree->SetBranchAddress("ohcalph",ohcalph);
  tree->SetBranchAddress("mbdtype",mbdtype);
  tree->SetBranchAddress("mbdside",mbdside);
  tree->SetBranchAddress("mbdchan",mbdchan);
  tree->SetBranchAddress("sectorem",&sectorem);
  tree->SetBranchAddress("sectorih",&sectorih);
  tree->SetBranchAddress("sectoroh",&sectoroh);
  tree->SetBranchAddress("sectormb",&sectormb);
  simt->SetBranchAddress("truthpar_e",tpe);
  simt->SetBranchAddress("truthpar_et",tpet);
  simt->SetBranchAddress("truthpar_ph",tpph);
  simt->SetBranchAddress("truthpar_em",tpem);
  simt->SetBranchAddress("truthpar_n",&tpn);
  simt->SetBranchAddress("npart",&npart);
  simt->SetBranchAddress("ncoll",&ncoll);
  simt->SetBranchAddress("bimp",&bimp);
  simt->SetBranchAddress("mbenrgy",smbenrgy);
  simt->SetBranchAddress("emcalen",semcalen);
  simt->SetBranchAddress("emcalet",semcalet);
  simt->SetBranchAddress("emcalph",semcalph);
  simt->SetBranchAddress("ihcalen",sihcalen);
  simt->SetBranchAddress("ihcalet",sihcalet);
  simt->SetBranchAddress("ihcalph",sihcalph);
  simt->SetBranchAddress("ohcalen",sohcalen);
  simt->SetBranchAddress("ohcalet",sohcalet);
  simt->SetBranchAddress("ohcalph",sohcalph);
  simt->SetBranchAddress("sectorem",&ssectorem);
  simt->SetBranchAddress("sectorih",&ssectorih);
  simt->SetBranchAddress("sectoroh",&ssectoroh);
  simt->SetBranchAddress("sectormb",&ssectormb);
  float z_v;
  simt->SetBranchAddress("zvtx",&z_v);
  TLine* lines[centbins];
  int fracdat = 50;
  int fracsim = 5;
  int counter = 0;
  int mbdsum = 0;
  int mbdsumunc = 0;
  //MBD timing channels 54 and 56 are bad
  cout << simt->GetEntries()/fracsim << endl;
  cout << tree->GetEntries()/fracdat << endl;
  /*
  for(int i=0; i<simt->GetEntries()/fracsim; ++i)
    {
      cout << i << " " <<endl
	}
  */
  for(int i=0; i<simt->GetEntries()/fracsim; ++i)
    {
      mbdsum = 0;
      mbdsumunc = 0;
      simt->GetEntry(i);
      /*
      for(int j=0; j<ssectormb; ++j)
	{
	  mbdsum += npart//smbenrgy[j];
	}
      */
      mbdsum = npart;
      if(z_v == 0) continue;
      if(abs(z_v) > 30) continue;
      hist->Fill(mbdsum);
      //cout << mbdsum << endl;
      /*
      for(int j=0; j<sectormb; ++j)
	{
	  mbdsumunc += mbenrgy[j];
	  if(mbdtype[j] == 1)
	    {
	      if(mbdside[j] == 1)
		{
		  mbdsum += mbenrgy[j]*gaincorr[mbdchan[j]+64];
		}
	      else if(mbdside[j] == 0)
		{
		  mbdsum += mbenrgy[j]*gaincorr[mbdchan[j]];
		}
	    }
	}
      */
      //cout << mbdsum << endl;
      /*
      if(mbdsumunc == 0)
	{
	  //for (int j=0; j<sectormb; ++j) cout << mbenrgy[j];
	  //cout << endl;
	  //cout << "zero uncorrected mbd sum in event " << i << endl;
	  //cout << "sectormb: " << sectormb << endl;
	  continue;
	}
      */
    }
  for(int i=0; i<tree->GetEntries()/fracdat; ++i)
    {
      mbdsum = 0;
      mbdsumunc = 0;
      tree->GetEntry(i);
      int mbdts = 0;
      int mbdtn = 0;
      int mbdns = 0;
      int mbdnn = 0;
      for(int j=0; j<sectormb; ++j)
	{
	  if(mbdtype[j] == 1)
	    {
	      if(mbdside[j] == 1)
		{
		  mbdsum += mbenrgy[j]*gaincorr[mbdchan[j]+64];
		}
	      else if(mbdside[j] == 0)
		{
		  mbdsum += mbenrgy[j]*gaincorr[mbdchan[j]];
		}
	    }
	  else
	    {
	      if(mbdside[j] == 1 && mbenrgy[j-8] > 10)
		{
		  mbdnn++;
		  mbdtn += mbenrgy[j]*9./5000-tq_t0_offsets[mbdchan[j]+64];
		}
	      else if(mbdside[j] == 0 && mbenrgy[j-8] > 10)
		{
		  mbdns++;
		  mbdts += mbenrgy[j]*9./5000-tq_t0_offsets[mbdchan[j]];
		}
	    }
	}
      if(mbdns == 0 || mbdnn == 0) continue;
      mbdts /= mbdns;
      mbdtn /= mbdnn;
      if(mbdsum<=0) continue;
      if(isnan(mbdtn-mbdts)) continue;
      if(abs(mbdtn-mbdts)*15>30) continue;
      if(mbdsum > 0) dmbh->Fill(mbdsum);
      cout << i << endl;
    }
  
  
  int n = 0;
  int nsum = 0;
  cout << "test" << endl;
  while(n<centbins)
    {
      nsum = 0;
      for(int i=1; i<hist->GetNbinsX()+1; ++i)
	{
	  //cout <<n <<" "<<i << " " << nsum << endl;
	  nsum += hist->GetBinContent(i);
	  if(n*hist->GetEntries()/centbins < nsum)
	    {
	      lines[n] = new TLine(hist->GetBinCenter(i),0,hist->GetBinCenter(i), hist->GetBinContent(i));
	      lines[n]->SetLineColor(kRed);
	      cents[n] = hist->GetBinLowEdge(i+1);
	      if(n==0) cents[n] = 0;
	      //cout << cents[n] << endl;
	      ++n;
	      break;
	    }
	}
    }

  n = 0;
  nsum = 0;
  cout << "test" << endl;
  cout << hist->GetEntries() << endl;
  cout << dmbh->GetEntries() << endl;
  cout << dmbh->GetBinContent(dmbh->GetNbinsX()+1) << endl;
  while(n<centbins)
    {
      nsum = 0;
      for(int i=1; i<dmbh->GetNbinsX()+1; ++i)
	{
	  //cout <<n <<" "<<i << " " << nsum << endl;
	  nsum += dmbh->GetBinContent(i);
	  if(n*dmbh->GetEntries()/centbins < nsum)
	    {
	      cout << n << " " << nsum << " " << i << " " << dmbh->GetBinContent(i) << endl;
	      //lines[n] = new TLine(hist->GetBinCenter(i),0,hist->GetBinCenter(i), hist->GetBinContent(i));
	      //lines[n]->SetLineColor(kRed);
	      dcent[n] = dmbh->GetBinLowEdge(i+1);
	      if(n==0) dcent[n]=0;
	      //cout << cents[n] << endl;
	      ++n;
	      break;
	    }
	}
    }

  
  cents[centbins] = 99999999999;
  cout << "test4" << endl;
  int kcodes[12] = {0,1,2,3,4,-1,-5,-8,-10,+1,+2,+3};
  for(int i=0; i<centbins; ++i)
    {
      setcolorcent(centclass[i][0], kcodes, kRed, 10, i);
      setcolorcent(centclass[i][1], kcodes, kGreen, 10, i);
      setcolorcent(centclass[i][2], kcodes, kBlue, 10, i);
    }
  
  for(int i=0; i<simt->GetEntries()/fracsim; ++i)
    {
      float emee = 0;
      float emea = 0;
      float emed = 0;
      simt->GetEntry(i);
      for(int j=0; j<tpn; ++j)
	{
	  emea += tpe[j]*sin(2*atan(exp(-tpet[j])));
	  if(tpem[j])
	    {
	      emee += tpe[j]*sin(2*atan(exp(-tpet[j])));
	    }
	}
      for(int j=0; j<ssectorem; ++j)
	{
	  emed += semcalen[j]*sin(2*atan(exp(-(semcalet[j]-48)*0.024)));
	}
      ectee->Fill(emed, emee);
      ectea->Fill(emed, emea);
    }











    
  for(int i=0; i<simt->GetEntries()/fracsim; ++i)
    {
      for(int j=0; j<64; ++j)
	{
	  for(int k=0; k<24; ++k)
	    {
	      semtowercomb[j][k] = 0;
	    }
	}
      mbdsum = 0;
      float simsum = 0;
      float sime = 0;
      simt->GetEntry(i);
      if(z_v == 0) continue;
      if(abs(z_v) > 30) continue;
      mbdsum = npart;
      float eval = 0;
      //cout << mbdsum << endl;
      for(int j=0; j<centbins; ++j)
	{
	  if(mbdsum < cents[j+1])
	    {
	      //cout <<sectorem << endl;
	      //cout << j << endl;
	      int etabin;
	      for(int k=0; k<ssectorem; ++k)
		{
		  if(semcalen[k] < mine) continue;
		  
		  if (semcalet[k] < 8) { continue; }
		  if (semcalet[k] >= 9 && semcalet[k] <= 47 && semcalph[k] >= 32 && semcalph[k] <= 39) { continue; } 
		  if (semcalet[k] >= 9 && semcalet[k] <= 15 && semcalph[k] >= 40 && semcalph[k] <= 47) { continue; } 
		  if (semcalet[k] >= 32 && semcalet[k] <= 39 && semcalph[k] >= 0 && semcalph[k] <= 7) { continue; } 
		  if (semcalet[k] >= 16 && semcalet[k] <= 23 && semcalph[k] >= 16 && semcalph[k] <= 23) { continue; } 
		  if (semcalet[k] >= 88 && semcalet[k] <= 95 && semcalph[k] >= 24 && semcalph[k] <= 32) { continue; } 
		  if (semcalet[k] >= 88 && semcalet[k] <= 95 && semcalph[k] >= 40 && semcalph[k] <= 47) { continue; } 
		  if (semcalet[k] >= 88 && semcalet[k] <= 95 && semcalph[k] >= 56 && semcalph[k] <= 72) { continue; } 
		  if (semcalet[k] >= 64 && semcalet[k] <= 72 && semcalph[k] >= 88 && semcalph[k] <= 95) { continue; } 
		  if (semcalet[k] >= 88 && semcalet[k] <= 95 && semcalph[k] >= 88 && semcalph[k] <= 95) { continue; } 
		  if (semcalet[k] >= 88 && semcalet[k] <= 95 && semcalph[k] >= 104 && semcalph[k] <= 111) { continue; } 
		  if (semcalet[k] >= 9 && semcalet[k] <= 31 && semcalph[k] >= 136 && semcalph[k] <= 143) { continue; } 
		  if (semcalet[k] >= 16 && semcalet[k] <= 47 && semcalph[k] >= 144 && semcalph[k] <= 151) { continue; } 
		  if (semcalet[k] >= 9 && semcalet[k] <= 15 && semcalph[k] >= 152 && semcalph[k] <= 159) { continue; } 
		  if (semcalet[k] >= 80 && semcalet[k] <= 95 && semcalph[k] >= 152 && semcalph[k] <= 159) { continue; } 
		  if (semcalet[k] >= 88 && semcalet[k] <= 95 && semcalph[k] >= 168 && semcalph[k] <= 171) { continue; } 
		  if (semcalet[k] >= 9 && semcalet[k] <= 15 && semcalph[k] >= 208 && semcalph[k] <= 223) { continue; } 
		  if (semcalet[k] >= 9 && semcalet[k] <= 15 && semcalph[k] >= 231 && semcalph[k] <= 239) { continue; } 
		  if (semcalet[k] >= 24 && semcalet[k] <= 31 && semcalph[k] >= 208 && semcalph[k] <= 215) { continue; } 
		  if (semcalet[k] >= 88 && semcalet[k] <= 95 && semcalph[k] >= 208 && semcalph[k] <= 227) { continue; } 
		  if (semcalet[k] >= 80 && semcalet[k] <= 87 && semcalph[k] >= 224 && semcalph[k] <= 231) { continue; } 
		  if (semcalet[k] >= 88 && semcalet[k] <= 95 && semcalph[k] >= 232 && semcalph[k] <= 239) { continue; } 
		  if (semcalet[k] >= 48 && semcalet[k] <= 95 && semcalph[k] >= 240 && semcalph[k] <= 247) { continue; }
		  
		  etabin = semcalet[k];
		  eval = scale*(semcalen[k]-subtr)*sin(2*atan(exp(-(etabin-48)*0.024)));
		  detacente[j][etabin] += eval;
		  sedist->Fill(eval);
		  if(mbdsum > cents[8]) mc20tws[0]->Fill(eval);
		  if(mbdsum < cents[3]) lc30tws[0]->Fill(eval);
		  scenttow[0][j]->Fill(eval);
		  //floor(emcalet[k]/0.023)+48;
		  
		  sime += eval;
		  //cout << j << "etabin" << etabin << " detacente " << detacente[j][etabin] << endl;
		  semtowercomb[semcalph[k]/4][etabin/4] += eval;
		  //cout << emcalen[k] << endl;
		  if(semcalen[k] > 0)
		    {
		      nem[j][etabin]++;
		    }
		}
	      simsum += sime;
	      etsim->Fill(sime);
	      scentet[0][j]->Fill(sime);
	      if(mbdsum > cents[8]) mc20ets[0]->Fill(sime);
	      if(mbdsum < cents[3]) lc30ets[0]->Fill(sime);
	      sime=0;
	      for(int k=0; k<ssectorih; ++k)
		{
		  if(sihcalen[k] <= 0) continue;
		  etabin = sihcalet[k];//floor(ihcalet[k]/0.092)+12;
		  eval = scale*(sihcalen[k]-subtr)*sin(2*atan(exp(-(etabin-12)*0.096)));
		  isedist->Fill(eval);
		  scenttow[1][j]->Fill(eval);
		  if(mbdsum > cents[8]) mc20tws[1]->Fill(eval);
		  if(mbdsum < cents[3]) lc30tws[1]->Fill(eval);
		  detacenti[j][etabin] += eval;//-ihcalet[k])));
		  sime += eval;
		  semtowercomb[sihcalph[k]][etabin] += eval;
		  nih[j][etabin]++;
		}
	      simsum += sime;
	      ietsim->Fill(sime);
	      scentet[1][j]->Fill(sime);
	      if(mbdsum > cents[8]) mc20ets[1]->Fill(sime);
	      if(mbdsum < cents[3]) lc30ets[1]->Fill(sime);
	      sime=0;
	      for(int k=0; k<ssectoroh; ++k)
		{
		  if(sohcalen[k] <= 0) continue;
		  etabin = sohcalet[k];//floor(ohcalet[k]/0.092)+12;
		  eval = scale*(sohcalen[k]-subtr)*sin(2*atan(exp(-(etabin-12)*0.096)));
		  scenttow[2][j]->Fill(eval);
		  osedist->Fill(eval);
		  if(mbdsum > cents[8]) mc20tws[2]->Fill(eval);
		  if(mbdsum < cents[3]) lc30tws[2]->Fill(eval);
		  detacento[j][etabin] += eval;
		  sime += eval;
		  semtowercomb[sohcalph[k]][etabin] += eval;
		  if(sohcalen[k] > 0)
		    {
		      noh[j][etabin]++;
		    }
		}
	      simsum += sime;
	      oetsim->Fill(sime);
	      scentet[2][j]->Fill(sime);
	      if(mbdsum > cents[8]) mc20ets[2]->Fill(sime);
	      if(mbdsum < cents[3]) lc30ets[2]->Fill(sime);
	      ssumhist->Fill(simsum);
	      simsum = 0;
	      sime=0;
	      //cout << detacenti[j][12] << endl;
	      break;
	    }
	}
      for(int j=0; j<64; ++j)
	{
	  for(int k=0; k<24; ++k)
	    {
	      if(semtowercomb[j][k] < 0.001) continue;
	      ssumtow->Fill(semtowercomb[j][k]);
	    }
	}
    }
  //cout << detacenti[5][12] << endl;
  /*
  for(int i=0; i<10; ++i)
    {
      for(int j=0; j<96; ++j)
	{
	  //cout << i << " " << j << " " << detacente[i][j] << endl;
	  if(j>7)
	    {
	      detacente[i][j] /= (nem[i][j]);
	    }
	  if(j<24)
	    {
	      //cout << "ihcal " << i << " " << j << " " << detacenti[i][j] << endl;
	      detacenti[i][j] /= (nih[i][j]);
	      detacento[i][j] /= (noh[i][j]);
	    }
	}
    }
  */
  cout << "test5" << endl;
  for(int i=0; i<centbins; ++i)
    {
      for(int j=0; j<96; ++j)
	{
	  //cout << i << " " << j << " " << detacente[i][j] << endl;
	  centclass[i][0]->SetBinContent(j+1,detacente[i][j]);
	  detacente[i][j] = 0;
	  if(j<8) centclass[i][0]->SetBinContent(j+1,centclass[i][0]->Integral()/(2.2*92/centclass[i][0]->GetNbinsX()));
	  if(j<24)
	    {
	      centclass[i][1]->SetBinContent(j+1,detacenti[i][j]);
	      centclass[i][2]->SetBinContent(j+1,detacento[i][j]);
	      detacento[i][j] = 0;
	      detacenti[i][j] = 0;
	    }
	}
    }
  hist->GetXaxis()->SetTitle("MBD Charge Channel Sum [???]");
  hist->GetYaxis()->SetTitle("Counts");
  cout << "test2" << endl;
  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->SetLogy();
  c1->cd();
  c1->SetTicks(1);
  //hist->GetXaxis()->SetRangeUser(0.1,0000);
  hist->Draw();
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.675, 1, kBlack, 0.04);
  //drawText("Simulated reco MBD charge sum", 0.85,0.825,1, kBlack,0.04);
  //drawText("Red lines represent 10\%-ile bins in",0.85,0.775,1,kBlack,0.04);
  //drawText("centrality, with most central at right",0.85,0.725,1,kBlack,0.04);
  for(int i=0; i<centbins; ++i) lines[i]->Draw();
  c1->SaveAs("cent_edep.png");
  TCanvas* c2 = new TCanvas("c2","c2",1000,1000);
  c2->SetLogy();
  c2->SetTicks(1);
  for(int j=0; j<3; ++j)
    {
      float min = 999999999;
      float max = 0;
      for(int i=0; i<centbins; ++i)
	{
	  if(centclass[i][j]->GetMaximum() > max)
	    {
	      max = centclass[i][j]->GetMaximum();
	    }
	  if(centclass[i][j]->GetMinimum() < min)
	    {
	      min = centclass[i][j]->GetMinimum();
	    }
	}
      for(int i=0; i<centbins; ++i)
	{
	  for(int k=0; k<8; ++k)
	    {
	      centclass[i][0]->SetBinContent(k+1,0);
	    }
	}
      for(int i=0; i<centbins; ++i)
	{
	  c2->cd();
	  if(i==0)
	    {
	      centclass[i][j]->GetYaxis()->SetRangeUser(min/2,max*2);
	      centclass[i][j]->Draw();
	    }
	  else centclass[i][j]->Draw("SAME");
	}
      c2->SaveAs(("centdETdeta"+to_string(j)+".png").c_str());
    }

  float date = 0;
  for(int i=0; i<tree->GetEntries()/fracdat; ++i)
    {
      float zvtx = 0;
      mbdsum = 0;
      float mbdtn = 0;
      float mbdts = 0;
      int mbdnn = 0;
      int mbdns = 0;
      float eval = 0;
      date=0;
      tree->GetEntry(i);
      for(int j=0; j<sectormb; ++j)
	{
	  if(mbdtype[j] == 1)
	    {
	      if(mbdside[j] == 1)
		{
		  mbdsum += mbenrgy[j]*gaincorr[mbdchan[j]+64];
		}
	      else if(mbdside[j] == 0)
		{
		  mbdsum += mbenrgy[j]*gaincorr[mbdchan[j]];
		}
	    }
	  else
	    {
	      if(mbdside[j] == 1 && mbenrgy[j-8] > 10)
		{
		  mbdnn++;
		  mbdtn += mbenrgy[j]*9./5000-tq_t0_offsets[mbdchan[j]+64];
		}
	      else if(mbdside[j] == 0 && mbenrgy[j-8] > 10)
		{
		  mbdns++;
		  mbdts += mbenrgy[j]*9./5000-tq_t0_offsets[mbdchan[j]];
		}
	    }
	}
      mbdts /= mbdns;
      mbdtn /= mbdnn;
      if(mbdsum<=0) continue;
      if(isnan(mbdtn-mbdts)) continue;
      if(abs(mbdtn-mbdts)*15>30) continue;
      else
      {
      	  mthist->Fill(mbdtn-mbdts);
	  //cout << mbdtn - mbdts << endl;
      }
      //cout << mbdsum << endl;
      //for(int j=0; j<10; ++j)
      //if(sectorem < 10) continue;
      for(int j=0; j<64; ++j)
	{
	  for(int k=0; k<24; ++k)
	    {
	      demtowercomb[j][k] = 0;
	    }
	}
      for(int j=0; j<centbins; ++j)
      {
	if(mbdsum < dcent[j+1])
	{
	      //cout <<sectorem << endl;
	      //cout << j << endl;
	      int etabin;
	      float datsum = 0;
	      for(int k=0; k<sectorem; ++k)
		{
		  if(emcalen[k] < mine) continue;
		  etabin = emcalet[k];
		  eval = (emcalen[k]-subtr)*sin(2*atan(exp(-(etabin-48)*0.024)));
		  dedist->Fill(eval);
		  demtowercomb[emcalph[k]/4][etabin/4] += eval;
		  //floor(emcalet[k]/0.023)+48;
		  //detacente[j][etabin] += emcalen[k]*sin(2*atan(exp(-(etabin-48)*0.024)));//-emcalet[k])));
		  dcenttow[0][j]->Fill(eval);
		  if(mbdsum > dcent[8]) mc20twd[0]->Fill(eval);
		  if(mbdsum < dcent[3]) lc30twd[0]->Fill(eval);
		  date += eval;
		  //cout << j << "etabin" << etabin << " detacente " << detacente[j][etabin] << endl;
		  //cout << emcalen[k] << endl;
		  //if(emcalen[k] > 0)
		  //{
		  //nem[j][etabin]++;
		  //}
		}
	      //cout << date/sectorem << endl;
	      etdata->Fill(date);
	      dcentet[0][j]->Fill(date);
	      datsum += date;
	      if(mbdsum > dcent[8]) mc20etd[0]->Fill(date);
	      if(mbdsum < dcent[3]) lc30etd[0]->Fill(date);
	      date=0;
	      for(int k=0; k<sectorih; ++k)
		{
		  etabin = ihcalet[k];//floor(ihcalet[k]/0.092)+12;
		  eval = (ihcalen[k]-subtr)*sin(2*atan(exp(-(etabin-12)*0.096)));
		  idedist->Fill(eval);
		  date += eval;
		  dcenttow[1][j]->Fill(eval);
		  if(mbdsum > dcent[8]) mc20twd[1]->Fill(eval);
		  if(mbdsum < dcent[3]) lc30twd[1]->Fill(eval);
		  demtowercomb[ihcalph[k]][etabin] += eval;
		  //detacenti[j][etabin] += ihcalen[k]*sin(2*atan(exp(-(etabin-12)*0.096)));//-ihcalet[k])));
		  //if(ihcalen[k] > 0)
		  //{
		  //nih[j][etabin]++;
		  //}
		}
	      ietdata->Fill(date);
	      dcentet[1][j]->Fill(date);
	      if(mbdsum > dcent[8]) mc20etd[1]->Fill(date);
	      if(mbdsum < dcent[3]) lc30etd[1]->Fill(date);
	      datsum += date;
	      date=0;
	      for(int k=0; k<sectoroh; ++k)
		{
		  etabin = ohcalet[k];//floor(ohcalet[k]/0.092)+12;
		  eval = (ohcalen[k]-subtr)*sin(2*atan(exp(-(etabin-12)*0.096)));
		  odedist->Fill(eval);
		  date += eval;
		  dcenttow[2][j]->Fill(eval);
		  if(mbdsum > dcent[8]) mc20twd[2]->Fill(eval);
		  if(mbdsum < dcent[3]) lc30twd[2]->Fill(eval);
		  demtowercomb[ohcalph[k]][etabin] += eval;
		  //detacento[j][etabin] += ohcalen[k]*sin(2*atan(exp(-(etabin-12)*0.096)));//-ohcalet[k])));
		  //if(ohcalen[k] > 0)
		  //{
		  //noh[j][etabin]++;
		  //}
		}
	      datsum += date;
	      if(mbdsum > dcent[8]) mc20etd[2]->Fill(date);
	      if(mbdsum < dcent[3]) lc30etd[2]->Fill(date);
	      oetdata->Fill(date);
	      dcentet[2][j]->Fill(date);
	      dsumhist->Fill(datsum);
	      datsum = 0;
	      date=0;
	      for(int k=0; k<64; ++k)
		{
		  for(int l=0; l<24; ++l)
		    {
		      if(demtowercomb[k][l] < 0.001) continue;
		      dsumtow->Fill(demtowercomb[k][l]);
		    }
		}
	      //if(date<0.001) continue;
	      
	      break;
	}
	else if(j>=9)
	  {
	    cout << "Somehow the MBD sum is not < an absurdly large value" << endl;
	  }
      }
	
    }

  TCanvas* c5 = new TCanvas("","",1500,1000);
  TCanvas* c7 =  new TCanvas("","",500,500);
  c5->Divide(3,2);
  c5->cd(1);
  gPad->SetLogy();
  gPad->SetTicks(1);
  etdata->Scale(1./etdata->Integral());
  etsim->Scale(1./etsim->Integral());
  etdata->SetLineColor(kRed);
  etsim->SetLineColor(kBlue);
  float maxval = max(etdata->GetMaximum(),etsim->GetMaximum());
  float minval = min(etdata->GetBinContent(etdata->FindLastBinAbove(0,1)),etsim->GetBinContent(etsim->FindLastBinAbove(0,1)));
  etdata->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  etdata->GetYaxis()->SetLabelSize(0.025);
  etdata->GetXaxis()->SetLabelSize(0.025);
  etdata->GetYaxis()->SetTitle("Counts");
  etdata->GetXaxis()->SetTitle("E_{T,total EMCal} [GeV]");
  etdata->Draw();
  etsim->Draw("SAME");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText("Red: data",0.85,0.76,1,kBlack,0.025);
  drawText("Blue: sim",0.85,0.73,1,kBlack,0.025);
  stringstream stream;
  stream << std::fixed << std::setprecision(2) << scale;
  string sc = stream.str();
  stringstream stream2;
  stream2 <<std::fixed << std::setprecision(0) << subtr*1000;
  string sub = stream2.str();
  drawText((sub+" MeV subtracted from each tower").c_str(),0.85,0.85,1,kBlack,0.025);
  drawText(("Sim scaled "+sc).c_str(),0.85,0.79,1,kBlack,0.025);
  drawText("Area normed to 1",0.85,0.82,1,kBlack,0.025);
  drawText("|z|<10cm",0.85,0.7,1,kBlack,0.025);
  doPlot(c7, etdata, etsim, sc, sub, 1, "emcal_et_event_" + sc + "_" + sub + ".png");

  for(int i=0; i<centbins; ++i)
    {
      plotcentet(c7, dcentet[0][i], scentet[0][i], i,i+1, 1, "EMCal", sc, sub, "21615");
      plotcentet(c7, dcentet[1][i], scentet[1][i], i,i+1, 1, "IHCal", sc, sub, "21615");
      plotcentet(c7, dcentet[2][i], scentet[2][i], i,i+1, 1, "OHCal", sc, sub, "21615");
      plotcenttow(c7, dcenttow[0][i], scenttow[0][i], i,i+1, 1, "EMCal", sc, sub, "21615");
      plotcenttow(c7, dcenttow[1][i], scenttow[1][i], i,i+1, 1, "IHCal", sc, sub, "21615");
      plotcenttow(c7, dcenttow[2][i], scenttow[2][i], i,i+1, 1, "OHCal", sc, sub, "21615");
    }
  string detstring[3] = {"EMCal","IHCal","OHCal"};
  for(int i=0; i<3; ++i)
    {
      plotcentet(c7, mc20etd[i], mc20ets[i], 8, centbins, 1, detstring[i], sc, sub, "21615");
      plotcentet(c7, lc30etd[i], lc30ets[i], 0, 3, 1, detstring[i], sc, sub, "21615");
      plotcenttow(c7, mc20twd[i], mc20tws[i], 8, centbins, 1, detstring[i], sc, sub, "21615");
      plotcenttow(c7, lc30twd[i], lc30tws[i], 0, 3, 1, detstring[i], sc, sub, "21615");
    }
    
  c5->cd(4);
  gPad->SetLogy();
  gPad->SetTicks(1);
  dedist->Scale(1./dedist->Integral());
  sedist->Scale(1./sedist->Integral());
  dedist->SetLineColor(kRed);
  sedist->SetLineColor(kBlue);
  maxval = max(dedist->GetMaximum(),sedist->GetMaximum());
  minval = min(dedist->GetBinContent(dedist->FindLastBinAbove(0,1)),sedist->GetBinContent(sedist->FindLastBinAbove(0,1)));
  dedist->GetXaxis()->SetRangeUser(0,10);
  dedist->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  //sedist->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  dedist->GetXaxis()->SetTitle("E_{T,EMCal tower} [GeV]");
  dedist->GetYaxis()->SetTitle("Counts");
  dedist->GetYaxis()->SetLabelSize(0.025);
  dedist->GetXaxis()->SetLabelSize(0.025);
  dedist->Draw();
  sedist->Draw("SAME");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText("Red: data",0.85,0.76,1,kBlack,0.025);
  drawText("Blue: sim",0.85,0.73,1,kBlack,0.025);
  drawText(("Sim scaled "+sc).c_str(),0.85,0.79,1,kBlack,0.025);
  drawText("Area normed to 1",0.85,0.82,1,kBlack,0.025);
  drawText("|z|<10cm",0.85,0.7,1,kBlack,0.025);
  drawText((sub+" MeV subtracted from each tower").c_str(),0.85,0.85,1,kBlack,0.025);
  doPlot(c7, dedist, sedist, sc, sub, 1, "emcal_et_tower" + sc + "_" + sub + ".png");


  c5->cd(2);
  gPad->SetLogy();
  gPad->SetTicks(1);
  ietdata->Scale(1./ietdata->Integral());
  ietsim->Scale(1./ietsim->Integral());
  ietdata->SetLineColor(kRed);
  ietsim->SetLineColor(kBlue);
  maxval = max(ietdata->GetMaximum(),ietsim->GetMaximum());
  minval = min(ietdata->GetBinContent(ietdata->FindLastBinAbove(0,1)),ietsim->GetBinContent(ietsim->FindLastBinAbove(0,1)));
  ietdata->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  ietdata->GetYaxis()->SetLabelSize(0.025);
  ietdata->GetXaxis()->SetLabelSize(0.025);
  ietdata->GetYaxis()->SetTitle("Counts");
  ietdata->GetXaxis()->SetTitle("E_{T,total IHCal} [GeV]");
  ietdata->Draw();
  ietsim->Draw("SAME");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText("Red: data",0.85,0.76,1,kBlack,0.025);
  drawText("Blue: sim",0.85,0.73,1,kBlack,0.025);
  drawText(("Sim scaled "+sc).c_str(),0.85,0.79,1,kBlack,0.025);
  drawText("Area normed to 1",0.85,0.82,1,kBlack,0.025);
  drawText("|z|<10cm",0.85,0.7,1,kBlack,0.025);
  drawText((sub+" MeV subtracted from each tower").c_str(),0.85,0.85,1,kBlack,0.025);

  doPlot(c7, ietdata, ietsim, sc, sub, 1, "ihcal_et_event_" + sc + "_" + sub + ".png");
  
  c5->cd(5);
  gPad->SetLogy();
  gPad->SetTicks(1);
  idedist->Scale(1./idedist->Integral());
  isedist->Scale(1./isedist->Integral());
  idedist->SetLineColor(kRed);
  isedist->SetLineColor(kBlue);
  maxval = max(idedist->GetMaximum(),isedist->GetMaximum());
  minval = min(idedist->GetBinContent(idedist->FindLastBinAbove(0,1)),isedist->GetBinContent(isedist->FindLastBinAbove(0,1)));
  idedist->GetXaxis()->SetRangeUser(0,10);
  idedist->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  idedist->GetXaxis()->SetTitle("E_{T,IHCal tower} [GeV]");
  idedist->GetYaxis()->SetTitle("Counts");
  idedist->GetYaxis()->SetLabelSize(0.025);
  idedist->GetXaxis()->SetLabelSize(0.025);
  idedist->Draw();
  isedist->Draw("SAME");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText("Red: data",0.85,0.76,1,kBlack,0.025);
  drawText("Blue: sim",0.85,0.73,1,kBlack,0.025);
  drawText(("Sim scaled "+sc).c_str(),0.85,0.79,1,kBlack,0.025);
  drawText("Area normed to 1",0.85,0.82,1,kBlack,0.025);
  drawText("|z|<10cm",0.85,0.7,1,kBlack,0.025);
  drawText((sub+" MeV subtracted from each tower").c_str(),0.85,0.85,1,kBlack,0.025);

  doPlot(c7, idedist, isedist, sc, sub, 1, "ihcal_et_tower" + sc + "_" + sub + ".png");


  c5->cd(3);
  gPad->SetLogy();
  gPad->SetTicks(1);
  oetdata->Scale(1./oetdata->Integral());
  oetsim->Scale(1./oetsim->Integral());
  oetdata->SetLineColor(kRed);
  oetsim->SetLineColor(kBlue);
  maxval = max(oetdata->GetMaximum(),oetsim->GetMaximum());
  minval = min(oetdata->GetBinContent(oetdata->FindLastBinAbove(0,1)),oetsim->GetBinContent(oetsim->FindLastBinAbove(0,1)));
  oetdata->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  oetdata->GetYaxis()->SetLabelSize(0.025);
  oetdata->GetXaxis()->SetLabelSize(0.025);
  oetdata->GetYaxis()->SetTitle("Counts");
  oetdata->GetXaxis()->SetTitle("E_{T,total OHCal} [GeV]");
  oetdata->Draw();
  oetsim->Draw("SAME");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText("Red: data",0.85,0.76,1,kBlack,0.025);
  drawText("Blue: sim",0.85,0.73,1,kBlack,0.025);
  drawText(("Sim scaled "+sc).c_str(),0.85,0.79,1,kBlack,0.025);
  drawText("Area normed to 1",0.85,0.82,1,kBlack,0.025);
  drawText("|z|<10cm",0.85,0.7,1,kBlack,0.025);
  drawText((sub+" MeV subtracted from each tower").c_str(),0.85,0.85,1,kBlack,0.025);
  doPlot(c7, oetdata, oetsim, sc, sub, 1, "ohcal_et_event_" + sc + "_" + sub + ".png");

  
  c5->cd(6);
  gPad->SetLogy();
  gPad->SetTicks(1);
  odedist->Scale(1./odedist->Integral());
  osedist->Scale(1./osedist->Integral());
  odedist->SetLineColor(kRed);
  osedist->SetLineColor(kBlue);
  maxval = max(odedist->GetMaximum(),osedist->GetMaximum());
  minval = min(odedist->GetBinContent(odedist->FindLastBinAbove(0,1)),osedist->GetBinContent(osedist->FindLastBinAbove(0,1)));
  odedist->GetXaxis()->SetRangeUser(0,10);
  odedist->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  odedist->GetXaxis()->SetTitle("E_{T,OHCal tower} [GeV]");
  odedist->GetYaxis()->SetTitle("Counts");
  odedist->GetYaxis()->SetLabelSize(0.025);
  odedist->GetXaxis()->SetLabelSize(0.025);
  odedist->Draw();
  osedist->Draw("SAME");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText("Red: data",0.85,0.76,1,kBlack,0.025);
  drawText("Blue: sim",0.85,0.73,1,kBlack,0.025);
  drawText(("Sim scaled "+sc).c_str(),0.85,0.79,1,kBlack,0.025);
  drawText("Area normed to 1",0.85,0.82,1,kBlack,0.025);
  drawText("|z|<10cm",0.85,0.7,1,kBlack,0.025);
  drawText((sub+" MeV subtracted from each tower").c_str(),0.85,0.85,1,kBlack,0.025);
  c5->SaveAs("simanddata.png");

  doPlot(c7, odedist, osedist, sc, sub, 1, "ohcal_et_tower" + sc + "_" + sub + ".png");
  
  TCanvas* c6 = new TCanvas("","",1500,1000);
  c6->cd();
  gPad->SetLogy();
  gPad->SetTicks(1);
  dsumtow->Scale(1./dsumtow->Integral());
  ssumtow->Scale(1./ssumtow->Integral());
  dsumtow->SetLineColor(kRed);
  ssumtow->SetLineColor(kBlue);
  maxval = max(dsumtow->GetMaximum(),ssumtow->GetMaximum());
  minval = min(dsumtow->GetBinContent(dsumtow->FindLastBinAbove(0,1)),ssumtow->GetBinContent(ssumtow->FindLastBinAbove(0,1)));
  //dsumtow->GetXaxis()->SetRangeUser(0,10);
  dsumtow->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  dsumtow->GetXaxis()->SetTitle("E_{T,tower by tower sum} [GeV]");
  dsumtow->GetYaxis()->SetTitle("Counts");
  dsumtow->GetYaxis()->SetLabelSize(0.025);
  dsumtow->GetXaxis()->SetLabelSize(0.025);
  dsumtow->Draw();
  ssumtow->Draw("SAME");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText("Red: data",0.85,0.76,1,kBlack,0.025);
  drawText("Blue: sim",0.85,0.73,1,kBlack,0.025);
  drawText(("Sim scaled "+sc).c_str(),0.85,0.79,1,kBlack,0.025);
  drawText("Area normed to 1",0.85,0.82,1,kBlack,0.025);
  drawText("|z|<10cm",0.85,0.7,1,kBlack,0.025);
  drawText((sub+" MeV subtracted from each tower").c_str(),0.85,0.85,1,kBlack,0.025);
  drawText("Overlapping calorimeter towers summed",0.85,0.67,1,kBlack,0.025);
  c6->SaveAs("towersum.png");
  doPlot(c7, dsumtow, ssumtow, sc, sub, 1, "sum_et_tower" + sc + "_" + sub + ".png");

  for(int i=1; i<max(dsumtow->FindLastBinAbove(0,1),ssumtow->FindLastBinAbove(0,1))+1; ++i)
    {
      if(dsumtow->GetBinContent(i) <= 0) continue;
      if(ssumtow->GetBinContent(i) <=0) continue;
      rsumtow->SetBinContent(i,dsumtow->GetBinContent(i)/ssumtow->GetBinContent(i));
    }
  rsumtow->SetMarkerStyle(39); rsumtow->SetMarkerColor(kRed);
  rsumtow->GetXaxis()->SetTitle("E_{T,tower by tower sum}");
  rsumtow->GetYaxis()->SetTitle("Data/sim");
  rsumtow->GetYaxis()->SetLabelSize(0.025);
  rsumtow->GetXaxis()->SetLabelSize(0.025);
  //rsumtow->GetYaxis()->SetRangeUser(0,2);
  gPad->SetLogy();
  rsumtow->Draw("P");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  c6->SaveAs("tsumrat.png");
  
  doPlot(c7, rsumtow, NULL, sc, sub, 1, "ratio_et_tower" + sc + "_" + sub + ".png");
  

  c6->cd();
  gPad->SetLogy();
  gPad->SetTicks(1);
  dsumhist->Scale(1./dsumhist->Integral());
  ssumhist->Scale(1./ssumhist->Integral());
  dsumhist->SetLineColor(kRed);
  ssumhist->SetLineColor(kBlue);
  maxval = max(dsumhist->GetMaximum(),ssumhist->GetMaximum());
  minval = min(dsumhist->GetBinContent(dsumhist->FindLastBinAbove(0,1)),ssumhist->GetBinContent(ssumhist->FindLastBinAbove(0,1)));
  //dsumhist->GetXaxis()->SetRangeUser(0,10);
  dsumhist->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  dsumhist->GetXaxis()->SetTitle("E_{T,calorimeter sum} [GeV]");
  dsumhist->GetYaxis()->SetTitle("Counts");
  dsumhist->GetYaxis()->SetLabelSize(0.025);
  dsumhist->GetXaxis()->SetLabelSize(0.025);
  dsumhist->Draw();
  ssumhist->Draw("SAME");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText("Red: data",0.85,0.76,1,kBlack,0.025);
  drawText("Blue: sim",0.85,0.73,1,kBlack,0.025);
  drawText(("Sim scaled "+sc).c_str(),0.85,0.79,1,kBlack,0.025);
  drawText("Area normed to 1",0.85,0.82,1,kBlack,0.025);
  drawText("|z|<10cm",0.85,0.7,1,kBlack,0.025);
  drawText((sub+" MeV subtracted from each tower").c_str(),0.85,0.85,1,kBlack,0.025);
  c6->SaveAs("calsum.png");
  doPlot(c7, dsumhist, ssumhist, sc, sub, 1, "sum_et_event" + sc + "_" + sub + ".png");
  

  for(int i=1; i<max(dsumhist->FindLastBinAbove(0,1),ssumhist->FindLastBinAbove(0,1))+1; ++i)
    {
      if(dsumhist->GetBinContent(i) <= 0) continue;
      if(ssumhist->GetBinContent(i) <=0) continue;
      rothist->SetBinContent(i,dsumhist->GetBinContent(i)/ssumhist->GetBinContent(i));
    }
  rothist->SetMarkerStyle(39); rothist->SetMarkerColor(kRed);
  rothist->GetXaxis()->SetTitle("E_{T,calorimeter sum}");
  rothist->GetYaxis()->SetTitle("Data/sim");
  rothist->GetYaxis()->SetLabelSize(0.025);
  rothist->GetXaxis()->SetLabelSize(0.025);
  //rothist->GetYaxis()->SetRangeUser(0,2);
  gPad->SetLogy(0);
  rothist->Draw("P");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  c6->SaveAs("sumrat.png");
  doPlot(c7, rothist, NULL, sc, sub, 1, "ratio_et_event" + sc + "_" + sub + ".png");

  c5->cd(1);
  for(int i=1; i<max(etdata->FindLastBinAbove(0,1),etsim->FindLastBinAbove(0,1))+1; ++i)
    {
      if(etdata->GetBinContent(i) <= 0) continue;
      if(etsim->GetBinContent(i) <=0) continue;
      rhist->SetBinContent(i,etdata->GetBinContent(i)/etsim->GetBinContent(i));
    }
  rhist->SetMarkerStyle(39); rhist->SetMarkerColor(kRed);
  rhist->GetXaxis()->SetTitle("E_{T,EMCal total}");
  rhist->GetYaxis()->SetTitle("Data/sim");
  rhist->GetYaxis()->SetLabelSize(0.025);
  rhist->GetXaxis()->SetLabelSize(0.025);
  //rhist->GetYaxis()->SetRangeUser(0,2);
  gPad->SetLogy();
  rhist->Draw("P");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);
  doPlot(c7, rhist, NULL, sc, sub, 1, "ratio_emcal_et_event" + sc + "_" + sub + ".png");

  
  c5->cd(4);
  for(int i=1; i<max(dedist->FindLastBinAbove(0,1),sedist->FindLastBinAbove(0,1))+1; ++i)
    {
      if(dedist->GetBinContent(i) <= 0) continue;
      if(sedist->GetBinContent(i) <=0) continue;
      rthist->SetBinContent(i,dedist->GetBinContent(i)/sedist->GetBinContent(i));
    }
  rthist->SetMarkerStyle(39); rthist->SetMarkerColor(kRed);
  rthist->GetXaxis()->SetTitle("E_{T,EMCal tower}");
  rthist->GetYaxis()->SetTitle("Data/sim");
  rthist->GetYaxis()->SetLabelSize(0.025);
  rthist->GetXaxis()->SetLabelSize(0.025);
  //rthist->GetYaxis()->SetRangeUser(0,2);
  gPad->SetLogy();
  rthist->Draw("P");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);

  doPlot(c7, rthist, NULL, sc, sub, 1, "ratio_emcal_et_tower" + sc + "_" + sub + ".png");
  
  c5->cd(2);
  for(int i=1; i<max(ietdata->FindLastBinAbove(0,1),ietsim->FindLastBinAbove(0,1))+1; ++i)
    {
      if(ietdata->GetBinContent(i) <= 0) continue;
      if(ietsim->GetBinContent(i) <=0) continue;
      rihist->SetBinContent(i,ietdata->GetBinContent(i)/ietsim->GetBinContent(i));
    }
  rihist->SetMarkerStyle(39); rihist->SetMarkerColor(kRed);
  rihist->GetXaxis()->SetTitle("E_{T,IHCal total}");
  rihist->GetYaxis()->SetTitle("Data/sim");
  rihist->GetYaxis()->SetLabelSize(0.025);
  rihist->GetXaxis()->SetLabelSize(0.025);
  //rihist->GetYaxis()->SetRangeUser(0,2);
  gPad->SetLogy();
  rihist->Draw("P");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);

  doPlot(c7, rihist, NULL, sc, sub, 1, "ratio_ihcal_et_event" + sc + "_" + sub + ".png");
  
  c5->cd(5);
  for(int i=1; i<max(idedist->FindLastBinAbove(0,1),isedist->FindLastBinAbove(0,1))+1; ++i)
    {
      if(idedist->GetBinContent(i) <= 0) continue;
      if(isedist->GetBinContent(i) <=0) continue;
      rithist->SetBinContent(i,idedist->GetBinContent(i)/isedist->GetBinContent(i));
    }
  rithist->SetMarkerStyle(39); rithist->SetMarkerColor(kRed);
  rithist->GetXaxis()->SetTitle("E_{T,IHCal tower}");
  rithist->GetYaxis()->SetTitle("Data/sim");
  rithist->GetYaxis()->SetLabelSize(0.025);
  rithist->GetXaxis()->SetLabelSize(0.025);
  //rithist->GetYaxis()->SetRangeUser(0,2);
  gPad->SetLogy();
  rithist->Draw("P");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);

  doPlot(c7, rithist, NULL, sc, sub, 1, "ratio_ihcal_et_tower" + sc + "_" + sub + ".png");
  
  c5->cd(3);
  for(int i=1; i<max(oetdata->FindLastBinAbove(0,1),oetsim->FindLastBinAbove(0,1))+1; ++i)
    {
      if(oetdata->GetBinContent(i) <= 0) continue;
      if(oetsim->GetBinContent(i) <=0) continue;
      rohist->SetBinContent(i,oetdata->GetBinContent(i)/oetsim->GetBinContent(i));
    }
  rohist->SetMarkerStyle(39); rohist->SetMarkerColor(kRed);
  rohist->GetXaxis()->SetTitle("E_{T,OHCal total}");
  rohist->GetYaxis()->SetTitle("Data/sim");
  rohist->GetYaxis()->SetLabelSize(0.025);
  rohist->GetXaxis()->SetLabelSize(0.025);
  //rohist->GetYaxis()->SetRangeUser(0,2);
  gPad->SetLogy();
  rohist->Draw("P");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);

  doPlot(c7, rohist, NULL, sc, sub, 1, "ratio_ohcal_et_event" + sc + "_" + sub + ".png");
  
  c5->cd(6);
  for(int i=1; i<max(odedist->FindLastBinAbove(0,1),osedist->FindLastBinAbove(0,1))+1; ++i)
    {
      if(odedist->GetBinContent(i) <= 0) continue;
      if(osedist->GetBinContent(i) <=0) continue;
      rothist->SetBinContent(i,odedist->GetBinContent(i)/osedist->GetBinContent(i));
    }
  rothist->SetMarkerStyle(39); rothist->SetMarkerColor(kRed);
  rothist->GetXaxis()->SetTitle("E_{T,OHCal tower}");
  rothist->GetYaxis()->SetTitle("Data/sim");
  rothist->GetYaxis()->SetLabelSize(0.025);
  rothist->GetXaxis()->SetLabelSize(0.025);
  //rothist->GetYaxis()->SetRangeUser(0,2);
  gPad->SetLogy();
  rothist->Draw("P");
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.93, 1, kBlack, 0.04);

  doPlot(c7, rothist, NULL, sc, sub, 1, "ratio_ohcal_et_tower" + sc + "_" + sub + ".png");
  
  c5->SaveAs("ratios.png");
  /*
  for(int i=0; i<10; ++i)
    {
      for(int j=0; j<96; ++j)
	{
	  cout << i << " " << j << " " << detacente[i][j] << endl;
	  if(j>7)
	    {
	      detacente[i][j] /= (nem[i][j]);
	    }
	  if(j<24)
	    {
	      detacenti[i][j] /= (nih[i][j]);
	      detacento[i][j] /= (noh[i][j]);
	    }
	}
    }
  */
  cout << "test5" << endl;
  /*
  for(int i=0; i<10; ++i)
    {
      for(int j=0; j<96; ++j)
	{
	  //cout << i << " " << j << " " << detacente[i][j] << endl;
	  centclass[i][0]->SetBinContent(j+1,detacente[i][j]);
	  if(j<24)
	    {
	      centclass[i][1]->SetBinContent(j+1,detacenti[i][j]);
	      centclass[i][2]->SetBinContent(j+1,detacento[i][j]);
	    }
	}
    }
  */
  hist->GetXaxis()->SetTitle("N_{Part}");
  hist->GetYaxis()->SetTitle("Counts");
  cout << "test2" << endl;
  //TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->SetLogy();
  c1->cd();
  c1->SetTicks(1);
  //hist->GetXaxis()->SetRangeUser(0.1,0000);
  hist->Draw();
  drawText("#bf{#it{sPHENIX}} internal", 0.85, 0.675, 1, kBlack, 0.04);
  //drawText("Simulated reco MBD charge sum", 0.85,0.825,1, kBlack,0.04);
  //drawText("Red lines represent 10\%-ile bins in",0.85,0.775,1,kBlack,0.04);
  //drawText("centrality, with most central at right",0.85,0.725,1,kBlack,0.04);
  for(int i=0; i<centbins; ++i) lines[i]->Draw();
  //c1->SaveAs("cent_edep.png");
  //TCanvas* c2 = new TCanvas("c2","c2",1000,1000);
  c2->SetLogy();
  c2->SetTicks(1);
  for(int i=0; i<centbins; ++i)
    {
      c2->cd();
      centclass[i][0]->Draw();
      centclass[i][1]->Draw("SAME");
      centclass[i][2]->Draw("SAME");
      c2->SaveAs(("centdETdeta_data_"+to_string(i)+".png").c_str());
    }
  TCanvas* c3 = new TCanvas("c3","c3",1000,500);
  c3->Divide(2,1);
  c3->cd(1);
  ectee->Draw("COLZ");
  c3->cd(2);
  ectea->Draw("COLZ");
  c3->SaveAs("emcalemvsall.png");
  TCanvas* c4 = new TCanvas("c4","c4");
  c4->cd();
  mthist->Draw();
  c4->SaveAs("mthist.png");

    for(int i=0; i<centbins; ++i)
    {
      cout << dcenttow[0][i]->GetEntries() << endl;
      cout << dcenttow[1][i]->GetEntries() << endl;
      cout << dcenttow[2][i]->GetEntries() << endl;

      cout << scentet[0][i]->GetEntries() << endl;
      cout << scentet[1][i]->GetEntries() << endl;
      cout << scentet[2][i]->GetEntries() << endl;

      cout << dcenttow[0][i]->GetBinContent(0) << endl;
      cout << dcenttow[1][i]->GetBinContent(0) << endl;
      cout << dcenttow[2][i]->GetBinContent(0) << endl;

      cout << scentet[0][i]->GetBinContent(0) << endl;
      cout << scentet[1][i]->GetBinContent(0) << endl;
      cout << scentet[2][i]->GetBinContent(0) << endl;
    }

  return 0;
}
