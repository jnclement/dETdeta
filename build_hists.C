#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TMath.h"
#include <utility>  // for std::pair
#include <cstdio>
#include <iostream>
#include "TGraph.h"
#include "TH2I.h"
#include <iomanip>
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include "stdlib.h"
#include <typeinfo>
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
#include "TProfile.h"
#include <TSystem.h>

const float eta_hc[] = {-1.05417,-0.9625,-0.870833,-0.779167,-0.6875,-0.595833,-0.504167,-0.4125,-0.320833,-0.229167,-0.1375,-0.0458333,0.0458333,0.1375,0.229167,0.320833,0.4125,0.504167,0.595833,0.6875,0.779167,0.870833,0.9625,1.05417};


  //{-1.104, -1.008, -0.912, -0.816, -0.72, -0.624, -0.528, -0.432, -0.336, -0.24, -0.144, -0.048, 0.048,
//0.144, 0.24, 0.336, 0.432, 0.528, 0.624, 0.72, 0.816, 0.912, 1.008, 1.104};
const float em_eta[] = {-1.12318,-1.10191,-1.08023,-1.05875,-1.03686,-1.01518,-0.993083,-0.971197,-0.948893,-0.926804,-0.904295,-0.882004,-0.859292,
-0.836801,-0.813886,-0.791195,-0.768079,-0.745189,-0.721872,-0.698783,-0.675264,-0.651976,-0.628256,-0.604768,-0.580845,-0.557155,-0.533027,
-0.509133,-0.484798,-0.460697,-0.43615,-0.411838,-0.387076,-0.362547,-0.337564,-0.312813,-0.287601,-0.26262,-0.237171,-0.211952,-0.186257,
-0.160788,-0.135722,-0.111796,-0.0874366,-0.0633955,-0.0389474,-0.0148467,0.0148467,0.0389474,0.0633955,0.0874366,0.111796,0.135722,0.160788,
0.186257,0.211952,0.237171,0.26262,0.287601,0.312813,0.337564,0.362547,0.387076,0.411838,0.43615,0.460697,0.484798,0.509133,0.533027,0.557155,
0.580845,0.604768,0.628256,0.651976,0.675264,0.698783,0.721872,0.745189,0.768079,0.791195,0.813886,0.836801,0.859292,0.882004,0.904295,0.926804,
0.948893,0.971197,0.993083,1.01518,1.03686,1.05875,1.08023,1.10191,1.12318};

int fullregonly(int phi)
{
  if(phi > 205 || phi < 165) return 1;
  return 0;
}

float fill_mbd_dat(int sectors, float* mbe, int* mbt, int* mbs, int* mbc, TH1* hist, float zcut, float zval, TH1* zhist, int cut, int datsim)
{
  float mbsum;//, ucmbd;
  //int mbdnn, mbdns;
  //float mbdtn, mbdts;
  float zvtx;
  //mbdnn=0;
  //mbdns=0;
  //mbdtn=0;
  //mbdts=0;
  mbsum=0;
  //ucmbd=0;
  for(int i=0; i<sectors; ++i)
    {
      //ucmbd += mbe[i];
      if(datsim)
	{
	  if(mbt[i] == 1)
	    {
	      if(mbs[i] == 1)
		{
		  mbsum += mbe[i];//*gaincorr[mbc[i]+64];
		}
	      else if(mbs[i] == 0)
		{
		  mbsum += mbe[i];//*gaincorr[mbc[i]];
		}
	    }
	}
      else mbsum += mbe[i];
      /*
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
      */
    }
  //if(mbdns == 0 || mbdnn == 0) return -1;
  //mbdts/=mbdns;
  //mbdtn/=mbdnn;
  zvtx = zval;
  if(abs(zval) == 0) return -1;
  //zvtx = (mbdtn-mbdts)*15;
  if(mbsum <=0)
    {
      return -1;
    }
  
  if(isnan(zvtx))
    {
      return -1;
    }
  if(zhist) zhist->Fill(zvtx);
  if(abs(zvtx) > zcut && cut)
    {
      return -1;
    }
  if(hist) hist->Fill(mbsum);
  return mbsum;
}

void set_em_combined_towers_0(float towers[96][24])
{
  for(int i=0; i<64; ++i)
    {
      for(int j = 0; j<24; ++j)
	{
	  towers[i][j] = 0;
	}
    }
}

int set_cent_cuts(TH1* hist, float* cent, int centbins)
{
  int n=0;
  int nsum=0;
  while(n<centbins)
    {
      //cout << n << endl;
      nsum = 0;
      for(int i=0; i<hist->GetNbinsX()+1; ++i)
	{
	  //if(j>=16) cout << i << endl;
	  if(i==hist->GetNbinsX())
	    {
	      cout << "Your centrality bins cannot be set for hist " << hist->GetName() << " centrality bin " << n << " because your hist range is not wide enough. All greater bins will also fail." << endl;
	      cent[n] = 999999;
	      return -1;
	      //exit(1);
	    }
	  nsum += hist->GetBinContent(i);
	  //if(n>=16) cout << nsum<< " " << (n+1)*hist->GetEntries()/centbins << endl;
	  if((n+1)*hist->GetEntries()/centbins <= nsum)
	    {
	      cent[n] = hist->GetBinLowEdge(i+1);
	      ++n;
	      break;
	    }
	}
    }
  cent[centbins-1] = 999999;
  return 0;
}

int check_acceptance(int eta, int phi)
{
  if (eta < 9) return 1;
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

float get_E_T_em(float E, float eta, float sub)
{
  return (E-sub)/cosh(eta);//*abs(sin(atan(exp(-eta))));//cosh((eta-47.5)*0.024);//*sin(2*atan(exp(-(eta-47.5)*0.024)));
}

float get_E_T_hc(float E, float eta, float sub)
{
  return (E-sub)/cosh(eta);//*abs(sin(atan(exp(-eta))));//cosh((eta-11.5)*0.096);//*sin(2*atan(exp(-(eta-11.5)*0.096)));
}


int build_hists(int simfrac = 1, int datfrac = 1, float zcut = 30, float simscale = 1.3, float subtracted = 0, float mine = 0, string tag="", int cor = 1)
{
  cout << "Starting..." << endl;
  mbd_init();
  float dETrange = 1.;
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  const int centbins = 18;
  const int centoffs = 2;
  const int hcalbins = 24;
  const int ecalbins = 96;
  //int mbd_bins[centbins+1] = {0};
  float mbenrgy[25000], calen[2][3][25000];
  int calet[2][3][25000], calph[2][3][25000];
  int   mbdtype[25000], mbdside[25000], mbdchan[25000];
  float towercomb[64][hcalbins];
  float etacor[2][3][25000];
  float simmbe[256];
  int simsecmb;
  int nevt[2] = {0};
  int sector[2][3];
  int sectormb;
  int truthpar_n;
  int dETbins = 200;
  float truthpar_eta[100000];
  float truthpar_e[100000];
  TH1D* truthpar_et[centbins];
  int npart = 0;
  float z_v[2][3];
  TFile* file = TFile::Open(("datatemp/merged_dEdeta"+tag+"_data_"+(cor?"cor":"unc")+"_71.root").c_str());
  TTree* tree[2];
  tree[1] = file->Get<TTree>("ttree");
  TFile* simf = TFile::Open(("datatemp/merged_dEdeta"+tag+"_mc_"+(cor?"cor":"unc")+"_555.root").c_str());
  tree[0] = simf->Get<TTree>("ttree");  
  float cents[2][centbins+centoffs] = {0};
  float truth_vtx[3];
  TH1D* centtow[2][3][centbins];
  TH1D* centet[2][3][centbins];
  TH1D* ettotcent[2][centbins];
  TH1D* dET[2][3];
  TH1D* dETcent[2][3][centbins];
  TH1D* truthparnhist = new TH1D("truthparnhist","",1000,0,10000);
  TH1D* truthparncent[centbins];
  TH1D* truthparehist = new TH1D("truthparehist","",100,0,50);
  TH1D* truthparecent[centbins];
  TH1D* truthpareetac[centbins];
  TH1D* meandiff[3];
  TH1D* sigmu[2][3];
  TH1D* meancent[2][3];
  TH2D* deadmap[2][3][centbins];
  TH1D* zcent[2][centbins];
  TH2I* deadhits[2][3][centbins];
  int phibins[3] = {256,64,64};
  int etabins[3] = {96,24,24};
  bool hit[3][hcalbins] = {false};
  
  TH1D* fullcor[3][centbins];
  for(int j=0; j<3; ++j)
    {
      meandiff[j] = new TH1D(("md"+to_string(j)).c_str(),"",centbins,0,90);
      for(int i=0; i<2; ++i)
	{
	  meancent[i][j] = new TH1D(("meancent"+to_string(i)+to_string(j)).c_str(),"",centbins,0,90);
	  sigmu[i][j] = new TH1D(("sigmu"+to_string(i)+to_string(j)).c_str(),"",centbins,0,90);
	  dET[i][j] = new TH1D(("dET"+to_string(i)+to_string(j)).c_str(),"",dETbins,-dETrange,dETrange);
	  for(int k=0; k<centbins; ++k)
	    {
	      if(i==0) fullcor[j][k] = new TH1D(("fullcor_"+to_string(j)+to_string(k)).c_str(),"",dETbins,-dETrange,dETrange);
	      dETcent[i][j][k] = new TH1D(("dETcent"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str(),"",dETbins,-dETrange,dETrange);
	      deadmap[i][j][k] = new TH2D(("deadmap"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str(),"",etabins[j],-0.5,etabins[j]-0.5,phibins[j],-0.5,phibins[j]-0.5);
	      deadhits[i][j][k] = new TH2I(("deadhits"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str(),"",etabins[j],-0.5,etabins[j]-0.5,phibins[j],-0.5,phibins[j]-0.5);
	      if(j==0) zcent[i][k] = new TH1D(("zcent"+to_string(i)+"_"+to_string(k)).c_str(),"",120,-30,30);
	    }
	}
    }
  
  float et_em_range[centbins] = {150,150,200,200,275,275,350,350,400,400,600,600,800,800,1200,1200,1750,1750};
  float et_oh_range[centbins] = {35,35,50,50,80,80,100,100,140,140,175,175,225,225,300,300,400,400};
  float et_ih_range[centbins] = {10,10,15,15,25,25,35,35,50,50,75,75,100,100,120,120,150,150};
  float et_sm_range = 2000;
  float tw_em_range = 15;
  float tw_oh_range = 20;
  float tw_ih_range = 2.5;
  float tw_sm_range = 20;
  int bins_tw = 100;
  int bins_et = 100;
  for(int i=0; i<centbins; ++i)
    {
      truthparecent[i] = new TH1D(("truthparecent_"+to_string(i)).c_str(),"",100,0,50);
      truthparncent[i] = new TH1D(("truthparncent_"+to_string(i)).c_str(),"",1000,0,10000);
      truthpareetac[i] = new TH1D(("truthpareetac_"+to_string(i)).c_str(),"",dETbins,-dETrange,dETrange);
      truthpar_et[i] = new TH1D(("truthpar_et_"+to_string(i)).c_str(),"",dETbins,-dETrange,dETrange);
      ettotcent[0][i] = new TH1D(("ettotcent0_" + to_string(i)).c_str(),"",400,0,2000);//et_em_range[centbins-1]);
      ettotcent[1][i] = new TH1D(("ettotcent1_" + to_string(i)).c_str(),"",400,0,2000);//et_em_range[centbins-1];
      centtow[1][0][i] = new TH1D(("centtow10_" + to_string(i)).c_str(),"",bins_tw,0,tw_em_range*(10.+i)/20);
      centtow[1][1][i] = new TH1D(("centtow11_" + to_string(i)).c_str(),"",bins_tw,0,tw_ih_range*(10.+i)/20);
      centtow[1][2][i] = new TH1D(("centtow12_" + to_string(i)).c_str(),"",bins_tw,0,tw_oh_range*(10.+i)/20);

      centtow[0][0][i] = new TH1D(("centtow00_" + to_string(i)).c_str(),"",bins_tw,0,tw_em_range*(10.+i)/20);
      centtow[0][1][i] = new TH1D(("centtow01_" + to_string(i)).c_str(),"",bins_tw,0,tw_ih_range*(10.+i)/20);
      centtow[0][2][i] = new TH1D(("centtow02_" + to_string(i)).c_str(),"",bins_tw,0,tw_oh_range*(10.+i)/20);

      centet[1][0][i] = new TH1D(("centet10_" + to_string(i)).c_str(),"",et_em_range[i],0,et_em_range[i]);
      centet[1][1][i] = new TH1D(("centet11_" + to_string(i)).c_str(),"",et_ih_range[i],0,et_ih_range[i]);
      centet[1][2][i] = new TH1D(("centet12_" + to_string(i)).c_str(),"",et_oh_range[i],0,et_oh_range[i]);

      centet[0][0][i] = new TH1D(("centet00_" + to_string(i)).c_str(),"",et_em_range[i],0,et_em_range[i]);
      centet[0][1][i] = new TH1D(("centet01_" + to_string(i)).c_str(),"",et_ih_range[i],0,et_ih_range[i]);
      centet[0][2][i] = new TH1D(("centet02_" + to_string(i)).c_str(),"",et_oh_range[i],0,et_oh_range[i]);
    }

  TH1D* ET[2][3];
  TH1D* TW[2][3];
  TH1D* sumev[2];
  TH1D* sumtw[2];
  TH1D* mbh[2];
  
  mbh[0] = new TH1D("smbh","",1000,0,300000);
  mbh[1] = new TH1D("dmbh","",1000,0,300000);
  ET[1][0] = new TH1D("et10","",2000/5,0,2000);
  ET[0][0] = new TH1D("et00","",2000/5,0,2000);
  TW[0][0] = new TH1D("tw00","",bins_tw,0,tw_em_range);
  TW[1][0] = new TH1D("tw10","",bins_tw,0,tw_em_range);
  ET[1][1] = new TH1D("et11","",200,0,200);
  ET[0][1] = new TH1D("et01","",200,0,200);
  TW[0][1] = new TH1D("tw01","",bins_tw,0,tw_ih_range);
  TW[1][1] = new TH1D("tw11","",bins_tw,0,tw_ih_range);
  ET[1][2] = new TH1D("et12","",600/3,0,600);
  ET[0][2] = new TH1D("et02","",600/3,0,600);
  TW[0][2] = new TH1D("tw02","",bins_tw,0,tw_oh_range);
  TW[1][2] = new TH1D("tw12","",bins_tw,0,tw_oh_range);
  sumev[0] = new TH1D("sumev0","",400,0,et_sm_range);
  sumev[1] = new TH1D("sumev1","",400,0,et_sm_range);
  sumtw[0] = new TH1D("sumtw0","",bins_tw,0,tw_sm_range);
  sumtw[1] = new TH1D("sumtw1","",bins_tw,0,tw_sm_range);
  
  cout << "Hists initialized" << endl;
  tree[1]->SetBranchAddress("mbenrgy",mbenrgy);
  tree[1]->SetBranchAddress("emcalen",calen[1][0]);
  tree[1]->SetBranchAddress("emcaletabin",calet[1][0]);
  tree[1]->SetBranchAddress("emcalphibin",calph[1][0]);
  tree[1]->SetBranchAddress("ihcalen",calen[1][1]);
  tree[1]->SetBranchAddress("ihcaletabin",calet[1][1]);
  tree[1]->SetBranchAddress("ihcalphibin",calph[1][1]);
  tree[1]->SetBranchAddress("ohcalen",calen[1][2]);
  tree[1]->SetBranchAddress("ohcaletabin",calet[1][2]);
  tree[1]->SetBranchAddress("ohcalphibin",calph[1][2]);
  tree[1]->SetBranchAddress("mbdtype",mbdtype);
  tree[1]->SetBranchAddress("mbdside",mbdside);
  tree[1]->SetBranchAddress("mbdchan",mbdchan);
  tree[1]->SetBranchAddress("sectorem",&sector[1][0]);
  tree[1]->SetBranchAddress("sectorih",&sector[1][1]);
  tree[1]->SetBranchAddress("sectoroh",&sector[1][2]);
  tree[1]->SetBranchAddress("sectormb",&sectormb);
  tree[1]->SetBranchAddress("track_vtx",z_v[1]);
  tree[1]->SetBranchAddress("emetacor",etacor[1][0]);
  tree[1]->SetBranchAddress("ihetacor",etacor[1][1]);
  tree[1]->SetBranchAddress("ohetacor",etacor[1][2]);
  tree[0]->SetBranchAddress("sectormb",&simsecmb);
  tree[0]->SetBranchAddress("mbenrgy",simmbe);
  tree[0]->SetBranchAddress("truth_vtx",truth_vtx);
  tree[0]->SetBranchAddress("emetacor",etacor[0][0]);
  tree[0]->SetBranchAddress("ihetacor",etacor[0][1]);
  tree[0]->SetBranchAddress("ohetacor",etacor[0][2]);
  tree[0]->SetBranchAddress("track_vtx",z_v[0]);
  tree[0]->SetBranchAddress("truthpar_n",&truthpar_n);
  tree[0]->SetBranchAddress("truthpar_eta",truthpar_eta);
  tree[0]->SetBranchAddress("truthpar_e",truthpar_e);
  tree[0]->SetBranchAddress("npart",&npart);
  tree[0]->SetBranchAddress("emcalen",calen[0][0]);
  tree[0]->SetBranchAddress("emcaletabin",calet[0][0]);
  tree[0]->SetBranchAddress("emcalphibin",calph[0][0]);
  tree[0]->SetBranchAddress("ihcalen",calen[0][1]);
  tree[0]->SetBranchAddress("ihcaletabin",calet[0][1]);
  tree[0]->SetBranchAddress("ihcalphibin",calph[0][1]);
  tree[0]->SetBranchAddress("ohcalen",calen[0][2]);
  tree[0]->SetBranchAddress("ohcaletabin",calet[0][2]);
  tree[0]->SetBranchAddress("ohcalphibin",calph[0][2]);
  tree[0]->SetBranchAddress("sectorem",&sector[0][0]);
  tree[0]->SetBranchAddress("sectorih",&sector[0][1]);
  tree[0]->SetBranchAddress("sectoroh",&sector[0][2]);
  cout << "Branches set" << endl;
  TH1D* zhist = new TH1D("zhist","",120,-30,30);
  TH1D* f10h[2][3];
  for(int i=0; i<2; ++i)
    {
      f10h[i][0] = new TH1D(("f10h"+to_string(i)+"0").c_str(),"",400,0,400);
      f10h[i][1] = new TH1D(("f10h"+to_string(i)+"1").c_str(),"",50,0,50);
      f10h[i][2] = new TH1D(("f10h"+to_string(i)+"2").c_str(),"",120,0,120);
    }
  float subtr = subtracted;//0.018;
  float scale[2];
  scale[0] = simscale;
  scale[1] = 1;
  int frac[2];
  int nevtcent[2][centbins] = {0};
  frac[0] = simfrac;
  frac[1] = datfrac;
  int run = 21615;
  const int par = 4;
  float parval[par];
  float mult[par];
  mult[0] = 1;
  mult[1] = 1000;
  mult[2] = 1000;
  mult[3] = 1;
  parval[0] = scale[0];
  parval[1] = subtr;
  parval[2] = mine;
  parval[3] = zcut;
  string params[par];
  stringstream streams[par];
  int precision[par];
  precision[0] = 2;
  precision[1] = 0;
  precision[2] = 0;
  precision[3] = 0;
  float test = 0;
  int counter[2][3] = {0};
  for(int i=0; i<par; ++i)
    {
      streams[i] << std::fixed << std::setprecision(precision[i]) << parval[i]*mult[i];
      params[i] = streams[i].str();
    }
  cout << "Now opening output file..." << endl;
  string outname = "datatemp/savedhists_fracsim_" + to_string(simfrac) + "_fracdat_" + to_string(datfrac) + "_subtr_" + params[1] + "_minE_" + params[2] + "_scale_" + params[0] + "_zcut_" + params[3] + "_run_"+to_string(run)+tag+"_"+(cor?"cor":"unc")+ ".root";
  TFile* outf = TFile::Open(outname.c_str(),"RECREATE");
  TTree* outt = new TTree("ttree","");
  float dummy;
  float eval;
  float allsum;
  float esum;
  float mbsum;
  int toprint[2] = {10000/frac[0],10000/frac[1]};
  
  cout << "Events for sim:  " << tree[0]->GetEntries()/frac[0] << endl;
  cout << "Events for data: " << tree[1]->GetEntries()/frac[1] << endl;
  cout << "Beginning processing." << endl;
  cout << "Filling MBD sim hist." << endl;
  for(int i=0; i<tree[0]->GetEntries()/frac[0]; ++i)
    {
      if(i%toprint[0] == 0) cout << "Doing event " << i << endl;
      tree[0]->GetEntry(i);
      /*
      if(abs(z_v[0][2]) == 0) continue;
      //if(abs(z_v[0][2]) > zcut) continue;
      if(npart == 0) continue;
      mbh[0]->Fill(npart);
      */
      dummy = fill_mbd_dat(simsecmb, simmbe, NULL, NULL, NULL, mbh[0], zcut, truth_vtx[2], zhist, 0, 0);
    }
  cout << "Sim MBD histogram entries: " << mbh[0]->GetEntries() << endl;
  cout << "Done filling sim MBD hist." << endl;
  cout << "Filling MBD data hist." << endl;
  for(int i=0; i<tree[1]->GetEntries()/frac[1]; ++i)
    {
      if(i%toprint[1] == 0) cout << "Doing event " << i << endl;
      tree[1]->GetEntry(i);
      dummy = fill_mbd_dat(sectormb, mbenrgy, mbdtype, mbdside, mbdchan, mbh[1], zcut, z_v[1][2], NULL, 0, 1);
    }
  cout << "Data MBD histogram entries: " << mbh[1]->GetEntries() << endl;
  cout << "Done filling data MBD hist." << endl;
  cout << "Done filling all MBD hists." << endl;
  cout << "Setting centrality bins" << endl;
  dummy = set_cent_cuts(mbh[0], cents[0], centbins+centoffs);
  cout << "Sim minbias hist entries: " << mbh[0]->GetEntries() << endl;
  dummy = set_cent_cuts(mbh[1], cents[1], centbins);
  //for(int i=0; i<centbins+1; ++i) cents[1][i] = mbd_bins[i];
  cout << "Data minbias hist entries: " << mbh[1]->GetEntries() << endl;
  cout << "Done setting centrality bins." << endl;
  cout << "cent bins sim/dat:" << endl;
  for(int i=0; i<centbins; ++i) cout << cents[0][i+2] << " " << cents[1][i] << endl;
  for(int h=0; h<2; ++h)
    {
      cout << "Doing tree " << h << "." << endl;
      for(int i=0; i<tree[h]->GetEntries()/frac[h]; ++i)
	{
	  if(i%toprint[h]==0) cout << "Starting event " << i << endl;
	  set_em_combined_towers_0(towercomb);
	  tree[h]->GetEntry(i);
	  if(h==0)
	    {
	      if(abs(z_v[0][2]) == 0) continue;
	      if(abs(z_v[0][2]) > zcut) continue;
	      mbsum = fill_mbd_dat(simsecmb, simmbe, NULL, NULL, NULL, NULL, zcut, z_v[0][2], NULL, 1, 0);
	    }
	  else mbsum = fill_mbd_dat(sectormb, mbenrgy, mbdtype, mbdside, mbdchan, NULL, zcut, z_v[1][2], NULL, 1, 1);
	  if(mbsum < 0) continue;
	  for(int j=0; j<centbins; ++j)
	    {
	      if(mbsum < cents[h][j+centoffs*(1-h)])
		{
		  nevt[h]++;
		  nevtcent[h][j]++;
		  //if(h==1 && j==17) cout << "j = 17 reached " << z_v[h] << endl;
		  zcent[h][j]->Fill(z_v[h][2]);
		  for(int k=0; k<3; ++k)
		    {
		      esum = 0;
		      for(int l=0; l<sector[h][k]; ++l)
			{
			  if(calen[h][k][l] < mine) continue;
			  if(k==0)
			    {
			      if(check_acceptance(calet[h][k][l], calph[h][k][l])) continue;
			      //if(fullregonly(calph[h][k][l])) continue;
			      eval = scale[h]*get_E_T_em(calen[h][k][l], etacor[h][k][l], subtr);
			    }
			  else eval = scale[h]*get_E_T_hc(calen[h][k][l], etacor[h][k][l], subtr);
			  if(calen[h][k][l] > 0.03)
			    {
			      deadmap[h][k][j]->Fill(calet[h][k][l],calph[h][k][l],calen[h][k][l]);
			      deadhits[h][k][j]->Fill(calet[h][k][l],calph[h][k][l]);
			    }
			  esum += eval;
			  TW[h][k]->Fill(eval);
			  centtow[h][k][j]->Fill(eval);
			  if(k==0) towercomb[calph[h][k][l]/4][calet[h][k][l]/4] += eval;
			  else towercomb[calph[h][k][l]][calet[h][k][l]] += eval;
			  dETcent[h][k][j]->Fill(etacor[h][k][l],eval);
			  dET[h][k]->Fill(etacor[h][k][l],eval);
			}
		      allsum += esum;
		      ET[h][k]->Fill(esum);
		      centet[h][k][j]->Fill(esum);
		      
		    }
		  sumev[h]->Fill(allsum);
		  ettotcent[h][j]->Fill(allsum);
		  for(int k=0; k<3; ++k)
		    {
		      for(int l=0; l<hcalbins; ++l)
			{
			  hit[k][l] = false;
			}
		    }
		  for(int k=0; k<64; ++k)
		    {
		      for(int l=0; l<hcalbins; ++l)
			{
			  if(towercomb[k][l] < mine) continue;
			  sumtw[h]->Fill(towercomb[k][l]);
			}
		    }
		  allsum = 0;
		  if(h==0)
		    {
		      int gtp = 0;
		      for(int k=0; k<truthpar_n; ++k)
			{
			  if(truthpar_eta[k] == 0 || abs(truthpar_eta[k]) > dETrange || truthpar_e[k] < mine) continue;
			  if(j==centbins-1) test += get_E_T_em(truthpar_e[k],truthpar_eta[k],0);
			  truthparehist->Fill(truthpar_e[k]);
			  truthparecent[j]->Fill(truthpar_e[k]);
			  truthpar_et[j]->Fill(truthpar_eta[k],get_E_T_em(truthpar_e[k],truthpar_eta[k],0));
			  cout << typeid(get_E_T_em(truthpar_e[k],truthpar_eta[k],0)).name() << endl;
			  truthpareetac[j]->Fill(truthpar_eta[k],truthpar_e[k]);
			  gtp++;
			}
		      truthparncent[j]->Fill(gtp);
		      truthparnhist->Fill(gtp);
		    }
		  break;
		}
	    }
	}
      cout << "Done." << endl;
    }
  centet[1][0][17]->Draw();
  gPad->SaveAs("test.png");
  for(int i=0; i<3; ++i)
    {
      for(int k=0; k<centbins; ++k)
	{
	  meandiff[i]->SetBinContent(centbins-k,(centet[1][i][k]->GetMean()-centet[0][i][k]->GetMean())/(centet[1][i][k]->GetMean()+centet[0][i][k]->GetMean()));
	  for(int j=0; j<2; ++j)
	    {
	      centet[j][i][k]->Scale(1./centet[j][i][k]->Integral());
	      TFitResultPtr fit = centet[j][i][k]->Fit("gaus","S");
	      sigmu[j][i]->SetBinContent(centbins-k,fit->Parameter(2));
	      sigmu[j][i]->SetBinError(centbins-k,fit->Error(2));
	      meancent[j][i]->SetBinContent(centbins-k,fit->Parameter(1));
	      meancent[j][i]->SetBinError(centbins-k,fit->Error(1));
	    }
	}
    }
  /*
  for(int i=0; i<centbins; ++i)
    {
      truthpar_et[i]->Divide(truthpar_counts[i]);
    }
  */
  cout << "Doing a few histogram operations..." << endl;

  for(int i=0; i<centbins; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  fullcor[j][i]->Divide(dETcent[1][j][i],dETcent[0][j][i]);
	  fullcor[j][i]->Multiply(truthpar_et[i]);
	  fullcor[j][i]->Scale(1./(nevtcent[1][i]));
	  cout << dETcent[1][j][i]->GetBinContent(10) << " " << dETcent[0][j][i]->GetBinContent(10) << " " << truthpar_et[i]->GetBinContent(10) << " " << fullcor[j][i]->GetBinContent(10)<< " " << fullcor[j][i]->Integral()/(2*dETrange*175) <<endl;
	  outf->WriteObject(fullcor[j][i],fullcor[j][i]->GetName());
	  if(i==centbins-1) cout << test << " " << dETcent[0][j][i]->Integral() << " " << dETcent[1][j][i]->Integral() << " " << fullcor[j][i]->Integral() << " " << truthpar_et[i]->Integral() << endl;
	  cout << typeid(test).name() << typeid(truthpar_et[i]).name() << endl;
	}
    }
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<3; ++i)
	{
	  for(int j=0; j<centbins; ++j)
	    {
	      deadmap[h][i][j]->Divide(deadhits[h][i][j]);
	    }
	  for(int k=0; k<centbins; ++k)
	    {
	      /*
	      for(int l=0; l<hcalbins; l++)
		{
		  cout << "Bin " << l << " in: dETcent " << dETcent[h][i][k]->GetBinContent(l+1) << " dETcentcount " << dETcentcount[h][i][k]->GetBinContent(l+1) << " ratio: " << dETcent[h][i][k]->GetBinContent(l+1)/dETcentcount[h][i][k]->GetBinContent(l+1) << endl;
		}
	      */
	      //dETcent[h][i][k]->Divide(dETcentcount[h][i][k]);
	      continue;
	    }
	}
    }
  cout << "Saving hists to " << outname << endl;
  outf->WriteObject(zhist, zhist->GetName());
  outf->WriteObject(truthparehist,truthparehist->GetName());
  outf->WriteObject(truthparnhist,truthparnhist->GetName());
  for(int i=0; i<centbins; ++i)
    {
      outf->WriteObject(truthparecent[i],truthparecent[i]->GetName());
      outf->WriteObject(truthparncent[i],truthparncent[i]->GetName());
      outf->WriteObject(truthpar_et[i],truthpar_et[i]->GetName());
      outf->WriteObject(truthpareetac[i],truthpareetac[i]->GetName());
    }
  for(int i=0; i<3; ++i)
    {
      outf->WriteObject(meandiff[i], meandiff[i]->GetName());
    }
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<3; ++i)
	{
	  outf->WriteObject(dET[h][i], dET[h][i]->GetName());
	  outf->WriteObject(sigmu[h][i], sigmu[h][i]->GetName());
	  outf->WriteObject(ET[h][i], ET[h][i]->GetName());
	  outf->WriteObject(TW[h][i], TW[h][i]->GetName());
	  outf->WriteObject(meancent[h][i],meancent[h][i]->GetName());
	}
    }
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<3; ++i)
	{
	  for(int j=0; j<centbins; j++)
	    {
	      if(i==0)
		{
		  outf->WriteObject(ettotcent[h][j], ettotcent[h][j]->GetName());
		  outf->WriteObject(zcent[h][j],zcent[h][j]->GetName());
		}
	      outf->WriteObject(deadmap[h][i][j], deadmap[h][i][j]->GetName());
	      outf->WriteObject(centtow[h][i][j], centtow[h][i][j]->GetName());
	      outf->WriteObject(centet[h][i][j], centet[h][i][j]->GetName());
	      outf->WriteObject(dETcent[h][i][j], dETcent[h][i][j]->GetName());
	    }
	}
    }
  for(int h=0; h<2; ++h)
    {
      outf->WriteObject(sumev[h], sumev[h]->GetName());
      outf->WriteObject(sumtw[h], sumtw[h]->GetName());
      outf->WriteObject(mbh[h], mbh[h]->GetName());
    }

  cout << "Done writing hists." << endl;
  cout << "Saving parameters to file." << endl;

  outt->Branch("sub",&subtr,"sub/F");
  outt->Branch("scale",scale,"scale[2]/F");
  outt->Branch("frac",frac,"frac[2]/I");
  outt->Branch("mine",&mine,"mine/F");
  int nccb = centbins;
  outt->Branch("nevtcent",nevtcent,"nevtcent[2][centbins]/I");
  outt->Branch("nevt",nevt,"nevt[2]/I");
  outt->Branch("cbin",&nccb,"cbin/I");
  outt->Branch("zcut",&zcut,"zcut/F");
  outt->Branch("run",&run,"run/I");
  outt->Fill();
  outf->WriteObject(outt,outt->GetName());
  cout << "Done writing parameters to file" << endl;
  cout << "All done!" << endl;

  TCanvas* c1 = new TCanvas("","");
  c1->cd();
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  f10h[i][j]->Draw();
	  c1->SaveAs(("f10h"+to_string(i)+to_string(j)+".pdf").c_str());
	}
    }
  
  return 0;
}
