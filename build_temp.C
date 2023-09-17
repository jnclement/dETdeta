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

float fill_mbd_dat(int sectors, float* mbe, int* mbt, int* mbs, int* mbc, TH1* hist)
{
  float mbsum, ucmbd;
  int mbdnn, mbdns;
  float mbdtn, mbdts;
  for(int i=0; i<sectors; ++i)
    {
      ucmbd += mbenrgy[i];
      if(mbt[i] == 1)
	{
	  if(mbs[i] == 1)
	    {
	      mbsum += mbe[i]*gaincorr[mbc[i]+64];
	    }
	  else if(mbs[i] == 0)
	    {
	      mbsum += mbe[i]*gaincorr[mbc[i]];
	    }
	}
      else
	{
	  if(mbs[i] == 1 && mbe[i-8] > 10)
	    {
	      mbdnn++;
	      mbdtn += mbe[i]*9./5000-tq_t0_offsets[mbc[i]+64];
	    }
	  else if(mbs[i] == 0 && mbe[i-8] > 10)
	    {
	      mbdns++;
	      mbdts += mbe[i]*9./5000-tq_t0_offsets[mbc[i]];
	    }
	}
    }
  mbdts/=mbdns;
  mbdtn/=mbdnn;
  if(mbdsum <=0) return -1;
  if(isnan(mbdtn-mbdts)) return -1;
  if(abs(mbdtn-mbdts) > 10) return -1;
  if(hist) hist->Fill(mbsum);
  return mbsum;
}

void set_em_combined_towers_0(float* towers)
{
  for(int i=0; i<64; ++i)
    {
      for(int j = 0; j<24; ++j)
	{
	  towers[i][j] = 0;
	}
    }
}

float get_et(float energy, float subtract, int eta)
{
  return (energy-subtract)*sin(2*atan(exp(-(eta-48)*0.024)));
}



void set_cent_cuts(TH1* hist, float* cent, int centbins)
{
  int n=0;
  int nsum=0;
  while(n<centbins)
    {
      nsum = 0;
      for(int i=1; i<hist->GetNbinsX()+1; ++i)
	{
	  if(i==hist->GetNbinsX())
	    {
	      cout << "Your centrality bins cannot be set because your hist range is not wide enough" << endl;
	      exit();
	    }
	  nsum += hist->GetBinContent(i);
	  if(n*hist->GetEntries()/centbins < nsum)
	    {
	      cent[n] = hist->GetBinLowEdge(i+1);
	      ++n;
	      break;
	    }
	}
    }
  cent[0] = 0;
  cent[centbins] = 999999999;
}

int check_acceptance(int eta, int phi)
{
  if (eta < 8) return 1;
  if (eta >= 9 && eta <= 47 && phi >= 32 && phi <= 39) return 1; 
  if (eta >= 9 && eta <= 15 && phi >= 40 && phi <= 47) return 1; 
  if (eta >= 32 && eta <= 39 && phi >= 0 && phi <= 7) return 1; 
  if (eta >= 16 && eta <= 23 && phi >= 16 && phi <= 23) return 1; 
  if (eta >= 88 && eta <= 95 && phi >= 24 && phi <= 32) return 1; 
  if (eta >= 88 && eta <= 95 && phi >= 40 && phi <= 47) return 1; 
  if (eta >= 88 && eta <= 95 && phi >= 56 && phi <= 72) return 1; 
  if (eta >= 64 && eta <= 72 && phi >= 88 && phi <= 95) return 1; 
  if (eta >= 88 && eta <= 95 && phi >= 88 && phi <= 95) return 1; 
  if (eta >= 88 && eta <= 95 && phi >= 104 && phi <= 111) return 1; 
  if (eta >= 9 && eta <= 31 && phi >= 136 && phi <= 143) return 1; 
  if (eta >= 16 && eta <= 47 && phi >= 144 && phi <= 151) return 1; 
  if (eta >= 9 && eta <= 15 && phi >= 152 && phi <= 159) return 1; 
  if (eta >= 80 && eta <= 95 && phi >= 152 && phi <= 159) return 1; 
  if (eta >= 88 && eta <= 95 && phi >= 168 && phi <= 171) return 1; 
  if (eta >= 9 && eta <= 15 && phi >= 208 && phi <= 223) return 1; 
  if (eta >= 9 && eta <= 15 && phi >= 231 && phi <= 239) return 1; 
  if (eta >= 24 && eta <= 31 && phi >= 208 && phi <= 215) return 1; 
  if (eta >= 88 && eta <= 95 && phi >= 208 && phi <= 227) return 1; 
  if (eta >= 80 && eta <= 87 && phi >= 224 && phi <= 231) return 1; 
  if (eta >= 88 && eta <= 95 && phi >= 232 && phi <= 239) return 1; 
  if (eta >= 48 && eta <= 95 && phi >= 240 && phi <= 247) return 1;
  return 0;
}

//float fill_cal_hist(TH1* tower, TH1* event, TH1* most_tow, TH1* least_tow, TH1* cent_tow, TH1* cent_evt, TH1* most_evt, TH1* least_evt, float* eme, TH1* mbdhist, int* eta, int sectors, int* cents, int centbins, 

float get_E_T_em(float E, int eta, float sub)
{
  return (E-subtr)*sin(2*atan(exp(-(etabin-48)*0.024)));
}

float get_E_T_hc(float E, int eta, float sub)
{
  return (E-subtr)*sin(2*atan(exp(-(etabin-12)*0.096)));
}


int build_hists()
{
  mbd_init();
  //cout << gaincorr[0] << endl;
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  const int centbins = 20;
  float mbenrgy[25000], ohcalen[25000], ihcalen[25000], emcalen[25000];
  int emcalet[25000], emcalph[25000];
  int ihcalet[25000], ihcalph[25000];
  int ohcalet[25000], ohcalph[25000];
  int   mbdtype[25000], mbdside[25000], mbdchan[25000];
  float smbenrgy[25000], sohcalen[25000], sihcalen[25000], semcalen[25000];
  int semcalet[25000], semcalph[25000];
  int sihcalet[25000], sihcalph[25000];
  int sohcalet[25000], sohcalph[25000];
  float cents[centbins+1] = {0};
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
  TFile* outf = TFile::Open("savedhists.root","RECREATE");
  TH1D* hist = new TH1D("hist","",500,0,500);
  TH1D* dmbh = new TH1D("dmbh","",3000,0,3000);
  double dcent[centbins+1] = {0};
  dcent[10] = 999999;
  TH1D* centclass[10][3];
  TH1D* dcenttow[3][10];
  TH1D* scenttow[3][10];
  TH1D* dcentet[3][10];
  TH1D* scentet[3][10];
  for(int i=0; i<10; ++i)
    {
      dcenttow[0][i] = new TH1D(("dcenttowem" + to_string(i)).c_str(),"",50,0,10);
      dcenttow[1][i] = new TH1D(("dcenttowih" + to_string(i)).c_str(),"",50,0,2);
      dcenttow[2][i] = new TH1D(("dcenttowoh" + to_string(i)).c_str(),"",50,0,10);

      scenttow[0][i] = new TH1D(("scenttowem" + to_string(i)).c_str(),"",50,0,10);
      scenttow[1][i] = new TH1D(("scenttowih" + to_string(i)).c_str(),"",50,0,2);
      scenttow[2][i] = new TH1D(("scenttowoh" + to_string(i)).c_str(),"",50,0,10);


      dcentet[0][i] = new TH1D(("dcentetem" + to_string(i)).c_str(),"",100,0,1200);
      dcentet[1][i] = new TH1D(("dcentetih" + to_string(i)).c_str(),"",100,0,100);
      dcentet[2][i] = new TH1D(("dcentetoh" + to_string(i)).c_str(),"",100,0,300);

      scentet[0][i] = new TH1D(("scentetem" + to_string(i)).c_str(),"",100,0,1200);
      scentet[1][i] = new TH1D(("scentetih" + to_string(i)).c_str(),"",100,0,100);
      scentet[2][i] = new TH1D(("scentetoh" + to_string(i)).c_str(),"",100,0,300);
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
  for(int i=0; i<10; ++i)
    {
      centclass[i][0] = new TH1D(("centhiste" + to_string(i)).c_str(),"",96,-1.1,1.1);
      centclass[i][1] = new TH1D(("centhisti" + to_string(i)).c_str(),"",24,-1.1,1.1);
      centclass[i][2] = new TH1D(("centhisto" + to_string(i)).c_str(),"",24,-1.1,1.1);
    }
  float detacente[10][96] = {0};
  float detacento[10][24] = {0};
  float detacenti[10][24] = {0};
  int nem[10][96] = {0};
  int nih[10][24] = {0};
  int noh[10][24] = {0};
  int tpn = 0;
  float subtr = 0.018;
  float scale = 1.3;
  float mine = 0.005;
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
  TLine* lines[10];
  int fracdat = 10;
  int fracsim = 1;
  int counter = 0;
  int mbdsum = 0;
  int mbdsumunc = 0;
