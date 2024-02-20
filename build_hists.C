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
#include <vector>
#include <TSystem.h>

std::set<std::tuple<int, int>> emcal_hot_dead_map_23696 = {{28,83},{39,115},{80, 228},{80, 236},{80, 237},{81, 236},{64, 164},{59,183},{72,188},{76,104},{88,104},{89, 104},{90, 104},{91, 104},
                                    {88,105},{89, 105},{90, 105},{91, 105},{88,108},{89, 108},{90, 107},{91, 107},{91,106},{88,111},{89,110},{89,111},{48,7},{58,7},{75,27},{9,200},{13,232},{22,157},
                                    {24,78},{24,208},{24,209},{24,210},{24,211},{24,212},{24,213},{24,214},{24,215},{25,208},{25,209},{25,210},{25,211},{25,212},{25,213},{25,214},{25,215},{26,208},
                                    {26,209},{26,210},{26,211},{26,212},{26,213},{26,214},{26,215},{27,208},{27,209},{27,210},{27,211},{27,212},{27,213},{27,214},{27,215},{28,83},{28,208},{28,209},
                                    {28,210},{28,211},{28,212},{28,213},{28,214},{28,215},{29,208},{29,209},{29,210},{29,211},{29,212},{29,213},{29,214},{29,215},{30,208},{30,209},{30,210},{30,211},
                                    {30,212},{30,213},{30,214},{30,215},{31,104},{31,208},{31,209},{31,210},{31,211},{31,212},{31,213},{31,214},{31,215},{32,144},{32,145},{32,146},{32,147},{32,148},
                                    {32,149},{32,150},{32,151},{33,144},{33,145},{33,146},{33,147},{33,148},{33,149},{33,150},{33,151},{34,144},{34,145},{34,146},{34,147},{34,148},{34,149},{34,150},
                                    {34,151},{35,144},{35,145},{35,146},{35,147},{35,148},{35,149},{35,150},{35,151},{36,144},{36,145},{36,146},{36,147},{36,148},{36,149},{36,150},{36,151},{37,144},
                                    {37,145},{37,146},{37,147},{37,148},{37,149},{37,150},{37,151},{38,144},{38,145},{38,146},{38,147},{38,148},{38,149},{38,150},{38,151},{38,219},{39,144},{39,145},
                                    {39,146},{39,147},{39,148},{39,149},{39,150},{39,151},{47,80},{47,96},{47,138},{48,215},{48,231},{48,253},{48,255},{63,77},{76,186},{76,187},{77,186},{77,187},
                                    {78,186},{78,187},{79,186},{79,187},{80,149},{88,40},{88,41},{88,42},{88,43},{88,168},{88,169},{88,170},{88,171},{88,224},{88,225},{88,226},{88,227},{89,40},
                                    {89,41},{89,42},{89,43},{89,168},{89,169},{89,170},{89,171},{89,224},{89,225},{89,226},{89,227},{90,40},{90,41},{90,42},{90,43},{90,111},{90,168},{90,169},
                                    {90,170},{90,171},{90,224},{90,225},{90,226},{90,227},{91,40},{91,41},{91,42},{91,43},{91,168},{91,169},{91,170},{91,171},{91,224},{91,225},{91,226},{91,227},
                                    {92,0},{92,40},{92,41},{92,42},{92,43},{92,168},{92,169},{92,170},{92,171},{92,224},{92,225},{92,226},{92,227},{93,40},{93,41},{93,42},{93,43},{93,168},
                                    {93,169},{93,170},{93,171},{93,224},{93,225},{93,226},{93,227},{94,40},{94,41},{94,42},{94,43},{94,96},{94,106},{94,157},{94,168},{94,169},{94,170},{94,171},
                                    {94,197},{94,219},{94,224},{94,225},{94,226},{94,227}};

int fullregonly(int phi)
{
  if(phi > 205 || phi < 165) return 1;
  return 0;
}

float fill_mbd_dat(int sectors, float* mbe, int* mbt, int* mbs, int* mbc, TH1* hist, float zcut, float zval, TH1* zhist, int cut, int datsim, int fillz)
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
      /*
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
	else*/
      mbsum += mbe[i];
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
  if(fillz) zhist->Fill(zvtx);
  if(abs(zvtx-zhist->GetMean()) > zcut && cut)
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
      for(int i=1; i<hist->GetNbinsX()+1; ++i)
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
  if (eta >= 64 && eta <= 71 && phi >= 64 && phi <= 71) return 1;
  if (eta < 9) return 1;
  if (eta < 48 && phi < 64) return 1;
  if (eta >= 48 && phi <=215 && phi >= 208) return 1;
  if (eta >= 32 && eta <=39 && phi <= 151 && phi >= 144) return 1;
  if (eta > 47 && phi <= 247 && phi >= 240) return 1;
  std::set<std::tuple<int,int>>::iterator it = emcal_hot_dead_map_23696.begin();
  while(it != emcal_hot_dead_map_23696.end())
    {
      if(std::get<0>(*it) == eta && std::get<1>(*it) == phi) return 1;
      ++it;
    }
  return 0;
  if ((eta == 80 || eta == 81) && phi == 237) return 1;
  if (phi >= 64 && phi <= 72 && eta <= 72 && eta >= 64) return 1;
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


int check_acc_map(TH2I* accmap, float mean, float std, int eta, int phi)
{
  float accbinc = accmap->GetBinContent(eta,phi);
  if(accbinc < (mean-2*std) || accbinc > (mean-2*std)) return 1;
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

bool check_eta_hit(float eta, vector<float> hits)
{
  float eps = 0.001;
  for(int i=0; i<hits.size(); ++i)
    {
      if(abs(hits.at(i) - eta) < eps)
	{
	  return true;
	}
    }
  return false;
}

int build_hists(int simfrac = 1, int datfrac = 1, float zcut = 30, float simscale = 1.3, float subtracted = 0, float mine = 0, string tag="", string tag2="", int cor = 1, float zlow = -30, float zup = 30)
{
  cout << "Starting..." << endl;
  //mbd_init();
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
  const int dETbins = 20;
  float truthpar_eta[100000];
  float truthpar_e[100000];
  TH1D* truthpar_et[centbins];
  int npart = 0;
  float z_v[2][3];
  TChain* tree[2];
  TChain* outt;
  TTree* outpt = new TTree("outpt","the output tree");

  tree[0] = new TChain("ttree");
  tree[1] = new TChain("ttree");
  outt = new TChain("outt");
  
  for(int i=0; i<1000/datfrac; ++i)
    {
      try
	{
      tree[1]->Add(("run/output/evt/events"+tag+"_data_cor_"+to_string(i)+".root").c_str());
	}
      catch(...)
	{
	  continue;
	}
    }

  for(int i=0; i<1000/simfrac; ++i)
    {
      try
	{
      tree[0]->Add(("run/output/evt/events"+tag2+"_mc_cor_"+to_string(i)+".root").c_str());
	}
      catch(...)
	{
	  continue;
	}
    }

  for(int i=0; i<1000; ++i)
    {
      try
	{
      outt->Add(("run/output/evt/events"+tag+"_data_cor_"+to_string(i)+".root").c_str());
	}
      catch(...)
	{
	  continue;
	}
    }
  
  float hotmap[3][96][256] = {0};

  outt->SetBranchAddress("hotmap",hotmap);
  outt->GetEntry(0);

  /*
  TFile* file = TFile::Open(("datatemp/merged_dEdeta"+tag+"_data_"+(cor?"cor":"unc")+"_600.root").c_str());
  TTree* tree[2];
  tree[1] = file->Get<TTree>("ttree");
  TFile* simf = TFile::Open(("datatemp/merged_dEdeta"+tag2+"_mc_"+(cor?"cor":"unc")+"_555.root").c_str());
  tree[0] = simf->Get<TTree>("ttree");
  TTree* hdtree = file->Get<TTree>("outt");
  */
  float cents[2][centbins+centoffs] = {0};
  float truth_vtx[3];
  TH1D* centtow[2][3][centbins];
  TH1D* centet[2][3][centbins];
  TH1D* ettotcent[2][centbins];
  TH1D* dET[2][3];
  TH1D* dETcent[2][3][centbins];
  TH1I* nfillcent[2][3][centbins];
  TH1D* dETcentrat[3][centbins];
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
  TH1D* dETcentsimunc[3][centbins];
  int phibins[3] = {256,64,64};
  int etabins[3] = {96,24,24};
  bool hit[3][hcalbins] = {false};
  vector<float> hits;
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
	      nfillcent[i][j][k] = new TH1I(("nfillcent_"+to_string(i)+"_"+to_string(j)+"_"+to_string(k)).c_str(),"",dETbins,-dETrange,dETrange);
	      if(i==0) dETcentsimunc[j][k] = new TH1D(("dETcentsimunc_"+to_string(j)+"_"+to_string(k)).c_str(),"",dETbins,-dETrange,dETrange);
	      if(i==0) fullcor[j][k] = new TH1D(("fullcor_"+to_string(j)+to_string(k)).c_str(),"",dETbins,-dETrange,dETrange);
	      dETcent[i][j][k] = new TH1D(("dETcent"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str(),"",dETbins,-dETrange,dETrange);
	      if(i==0) dETcentrat[j][k] = new TH1D(("dETcentrat"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str(),"",dETbins,-dETrange,dETrange);
	      deadmap[i][j][k] = new TH2D(("deadmap"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str(),"",etabins[j],-0.5,etabins[j]-0.5,phibins[j],-0.5,phibins[j]-0.5);
	      deadhits[i][j][k] = new TH2I(("deadhits"+to_string(i)+to_string(j)+"_"+to_string(k)).c_str(),"",etabins[j],-0.5,etabins[j]-0.5,phibins[j],-0.5,phibins[j]-0.5);
	      if(j==0) zcent[i][k] = new TH1D(("zcent"+to_string(i)+"_"+to_string(k)).c_str(),"",120,-30,30);
	    }
	}
    }
  
  float et_em_range[centbins] = {150,150,200,200,275,275,350,350,400,400,600,600,800,800,1200,1200,1750,1750};//{100,150,200,275,350,400,600,800,1200};//
  float et_oh_range[centbins] = {35,35,50,50,80,80,100,100,140,140,175,175,225,225,300,300,400,400};//{35,50,80,100,140,175,225,300,400};//
  float et_ih_range[centbins] = {10,10,15,15,25,25,35,35,50,50,75,75,100,100,120,120,150,150};//{10,15,25,35,50,75,100,120,150};//
  float et_sm_range = 1200;
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

      centet[1][0][i] = new TH1D(("centet10_" + to_string(i)).c_str(),"",et_em_range[i]+20,-20,et_em_range[i]);
      centet[1][1][i] = new TH1D(("centet11_" + to_string(i)).c_str(),"",et_ih_range[i]+10,-10,et_ih_range[i]);
      centet[1][2][i] = new TH1D(("centet12_" + to_string(i)).c_str(),"",et_oh_range[i]+20,-20,et_oh_range[i]);

      centet[0][0][i] = new TH1D(("centet00_" + to_string(i)).c_str(),"",et_em_range[i]+20,-20,et_em_range[i]);
      centet[0][1][i] = new TH1D(("centet01_" + to_string(i)).c_str(),"",et_ih_range[i]+10,-10,et_ih_range[i]);
      centet[0][2][i] = new TH1D(("centet02_" + to_string(i)).c_str(),"",et_oh_range[i]+20,-20,et_oh_range[i]);
    }

  TH1D* ET[2][3];
  TH1D* truthpar_total_ET = new TH1D("truthpar_total_ET","",2000,0,2000);
  TH1D* TW[2][3];
  TH1D* sumev[2];
  TH1D* sumtw[2];
  TH1D* mbh[2];
  TH1D* hcalraw[2][3][256];
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  for(int k=0; k<(j==0?256:64); ++k)
	    {
	      hcalraw[i][j][k] = new TH1D(("hcalraw_"+to_string(i)+"_"+to_string(j)+"_"+to_string(k)).c_str(),"",dETbins,-dETrange,dETrange);
	    }
	}
    }
  
  mbh[0] = new TH1D("smbh","",1000,0,300000);
  mbh[1] = new TH1D("dmbh","",1000,0,3000);
  ET[1][0] = new TH1D("et10","",1300,-100,1200);
  ET[0][0] = new TH1D("et00","",1300,-100,1200);
  TW[0][0] = new TH1D("tw00","",bins_tw,0,tw_em_range);
  TW[1][0] = new TH1D("tw10","",bins_tw,0,tw_em_range);
  ET[1][1] = new TH1D("et11","",220,-20,200);
  ET[0][1] = new TH1D("et01","",220,-20,200);
  TW[0][1] = new TH1D("tw01","",bins_tw,0,tw_ih_range);
  TW[1][1] = new TH1D("tw11","",bins_tw,0,tw_ih_range);
  ET[1][2] = new TH1D("et12","",660,-60,600);
  ET[0][2] = new TH1D("et02","",660,-60,600);
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
  TH1D* zhist[2];
  zhist[0] = new TH1D("zhist_0","",120,-30,30);
  zhist[1] = new TH1D("zhist_1","",120,-30,30);
  TH1D* f10h[2][3];
  for(int i=0; i<2; i++)
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
  int run = 23696;
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
  TH2I* accmaps[3];
  accmaps[0] = new TH2I("accmap0","",96,-0.5,95.5,256,-0.5,255.5);
  accmaps[1] = new TH2I("accmap1","",24,-0.5,23.5,96,-0.5,95.5);
  accmaps[2] = new TH2I("accmap2","",24,-0.5,23.5,96,-0.5,95.5);
  
  for(int i=0; i<par; ++i)
    {
      streams[i] << std::fixed << std::setprecision(precision[i]) << parval[i]*mult[i];
      params[i] = streams[i].str();
    }
  cout << "Now opening output file..." << endl;
  string outname = "datatemp/savedhists_fracsim_" + to_string(simfrac) + "_fracdat_" + to_string(datfrac) + "_subtr_" + params[1] + "_minE_" + params[2] + "_scale_" + params[0] + "_zcut_" + params[3]+tag+tag2+"_"+(cor?"cor":"unc")+"_zloup_"+to_string((int)zlow)+"_"+to_string((int)zup)+".root";
  TFile* outf = TFile::Open(outname.c_str(),"RECREATE");
  float dummy;
  float eval;
  float allsum;
  float esum;
  float mbsum;
  int toprint[2] = {10000/frac[0],10000/frac[1]};
  
  cout << "Events for sim:  " << tree[0]->GetEntries() << endl;
  cout << "Events for data: " << tree[1]->GetEntries() << endl;
  cout << "Beginning processing." << endl;
  cout << "Filling MBD sim hist." << endl;
  for(int i=0; i<tree[0]->GetEntries(); ++i)
    {
      if(i%toprint[0] == 0) cout << "Doing event " << i << endl;
      tree[0]->GetEntry(i);
      //if(z_v[0][2] < zlow || z_v[0][2] > zup) continue;
      /*
      if(abs(z_v[0][2]) == 0) continue;
      //if(abs(z_v[0][2]) > zcut) continue;
      if(npart == 0) continue;
      mbh[0]->Fill(npart);
      */
      dummy = fill_mbd_dat(simsecmb, simmbe, NULL, NULL, NULL, mbh[0], zcut, z_v[0][2], zhist[0], 0, 0, 1);
    }
  cout << "Sim MBD histogram entries: " << mbh[0]->GetEntries() << endl;
  cout << "Done filling sim MBD hist." << endl;
  cout << "Filling MBD data hist." << endl;
  for(int i=0; i<tree[1]->GetEntries(); ++i)
    {
      if(i%toprint[1] == 0) cout << "Doing event " << i << endl;
      tree[1]->GetEntry(i);
      //      if(z_v[1][2] < zlow || z_v[1][2] > zup) continue;
      dummy = fill_mbd_dat(sectormb, mbenrgy, mbdtype, mbdside, mbdchan, mbh[1], zcut, z_v[1][2], zhist[1], 0, 1, 1);
    }

  TFitResultPtr zfit[2];

  zfit[0] = zhist[0]->Fit("gaus","S");
  zfit[1] = zhist[1]->Fit("gaus","S");

  TF1* zfitf[2];
  zfitf[0] = zhist[0]->GetFunction("gaus");
  zfitf[1] = zhist[1]->GetFunction("gaus");
  
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
  for(int i=0; i<centbins; ++i) cout << cents[0][i+centoffs] << " " << cents[1][i] << endl;
  /*
  for(int i=0; i<tree[1]->GetEntries()/100; ++i)
    {
      if(i%toprint[1]==0) cout << "Filling accmap for event " << i << endl;
      tree[1]->GetEntry(i);
      mbsum = fill_mbd_dat(sectormb, mbenrgy, mbdtype, mbdside, mbdchan, NULL, zcut, z_v[1][2], zhist[0], 1, 1, 0);
      if(mbsum < 0) continue;
      for(int k=0; k<3; ++k)
	{
	  for(int l=0; l<sector[1][k]; ++l)
	    {
	      accmaps[k]->Fill(calet[1][k][l],calph[1][k][l]);
	    }
	}
    }
  float accrms[3] = {0};
  float accavg[3] = {0};
  
  for(int i=0; i<3; ++i)
    {
      accavg[i] = accmaps[i]->Integral()/(accmaps[i]->GetNbinsX()*accmaps[i]->GetNbinsY());
      for(int j=0; j<accmaps[i]->GetNbinsX(); ++j)
	{
	  for(int k=0; k<accmaps[i]->GetNbinsY(); k++)
	    {
	      accrms[i]+=pow(accmaps[i]->GetBinContent(j+1,k+1)-accavg[i],2)/(accmaps[i]->GetNbinsX()*accmaps[i]->GetNbinsY());
	    }
	}
      accrms[i] = sqrt(accrms[i]);
    }
  */
  /*
  TFile* hdf = TFile::Open("datatemp/hdm.root");
  TH2D* hdm[3];
  hdm[0] = (TH2D*)hdf->Get("hd0");
  hdm[1] = (TH2D*)hdf->Get("hd1");
  hdm[2] = (TH2D*)hdf->Get("hd2");
  */
  
  
  float totalweight;
  for(int h=0; h<2; ++h)
    {
      cout << "Doing tree " << h << "." << endl;
      for(int i=0; i<tree[h]->GetEntries(); ++i)
	{
	  if(i%toprint[h]==0) cout << "Starting event " << i << endl;
	  set_em_combined_towers_0(towercomb);
	  tree[h]->GetEntry(i);
	  if(h==0)
	    {
	      mbsum = fill_mbd_dat(simsecmb, simmbe, NULL, NULL, NULL, NULL, zcut, z_v[0][2], zhist[0], 1, 0, 0);
	    }
	  else mbsum = fill_mbd_dat(sectormb, mbenrgy, mbdtype, mbdside, mbdchan, NULL, zcut, z_v[1][2], zhist[0], 1, 1, 0);
	  if(mbsum < 0) continue;
	  if(z_v[h][2] < zlow || z_v[h][2] > zup) continue;
	  if(h==0 && mbsum < cents[0][centoffs-1]) continue;
	  for(int j=0; j<centbins; ++j)
	    {
	      if(mbsum < cents[h][j+centoffs*(1-h)])
		{
		  nevt[h]++;
		  nevtcent[h][j]++;
		  //if(h==1 && j==17) cout << "j = 17 reached " << z_v[h] << endl;
		  zcent[h][j]->Fill(z_v[h][2]);
		  float weight = 1;
		  if(h==0) weight = zfitf[1]->Eval(z_v[0][2])/zfitf[0]->Eval(z_v[0][2]);
		  if(h==0) totalweight += weight;
		  for(int k=0; k<3; ++k)
		    {
		      esum = 0;
		      for(int l=0; l<sector[h][k]; ++l)
			{
			  if(!check_eta_hit(etacor[h][k][l],hits))
			    {
			      nfillcent[h][k][j]->Fill(etacor[h][k][l],weight);
			      hits.push_back(etacor[h][k][l]);
			    }
			  if(calen[h][k][l] < mine) continue;
			  float eval_unc = scale[h]*get_E_T_em(calen[h][k][l], etacor[h][k][l], subtr);
			  if(h==0) dETcentsimunc[k][j]->Fill(etacor[h][k][l],eval_unc/(dETrange*2./dETbins));
			  //if(hotmap[k][calet[h][k][l]][calph[h][k][l]]) continue;
			  if(k==0)
			    {
			      if(check_acceptance(calet[h][k][l], calph[h][k][l])) continue;
			      //if(fullregonly(calph[h][k][l])) continue;
			      eval = scale[h]*get_E_T_em(calen[h][k][l], etacor[h][k][l], subtr);
			    }
			  else eval = scale[h]*get_E_T_hc(calen[h][k][l], etacor[h][k][l], subtr);

			  if(h==0)
			    {
			      eval *= weight;
			    }
			  
			  {
			    deadmap[h][k][j]->Fill(calet[h][k][l],calph[h][k][l],calen[h][k][l]);
			    deadhits[h][k][j]->Fill(calet[h][k][l],calph[h][k][l]);
			  }
			  if(j==centbins-1)
			    {
			      hcalraw[h][k][calph[h][k][l]]->Fill(etacor[h][k][l],eval/(dETrange*2./dETbins));
			    }
			  esum += eval;
			  TW[h][k]->Fill(eval);
			  centtow[h][k][j]->Fill(eval);
			  if(k==0) towercomb[calph[h][k][l]/4][calet[h][k][l]/4] += eval;
			  else towercomb[calph[h][k][l]][calet[h][k][l]] += eval;
			  dETcent[h][k][j]->Fill(etacor[h][k][l],eval/(dETrange*2./dETbins));
			  dET[h][k]->Fill(etacor[h][k][l],eval/(dETrange*2./dETbins));
			}
		      allsum += esum;
		      ET[h][k]->Fill(esum);
		      centet[h][k][j]->Fill(esum);
		      hits.clear();    
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
		      float total_et = 0;
		      for(int k=0; k<truthpar_n; ++k)
			{
			  if(truthpar_eta[k] == 0 || abs(truthpar_eta[k]) > dETrange || truthpar_e[k] < mine) continue;
			  if(j==centbins-1) test += get_E_T_em(truthpar_e[k],truthpar_eta[k],0);
			  truthparehist->Fill(truthpar_e[k]);
			  truthparecent[j]->Fill(truthpar_e[k]);
			  truthpar_et[j]->Fill(truthpar_eta[k],get_E_T_em(truthpar_e[k],truthpar_eta[k],0)/(dETrange*2./dETbins));
			  truthpareetac[j]->Fill(truthpar_eta[k],truthpar_e[k]);
			  total_et += get_E_T_em(truthpar_e[k], truthpar_eta[k],0);
			  gtp++;
			}
		      truthparncent[j]->Fill(gtp);
		      truthparnhist->Fill(gtp);
		      truthpar_total_ET->Fill(total_et);
		    }
		  break;
		}
	    }
	}
      cout << "Done." << endl;
    }

  for(int i=0; i<2; ++i)
    {
      sumtw[i]->Scale(1./totalweight);
      for(int j=0; j<3; ++j)
	{
	  TW[i][j]->Scale(1./totalweight);
	  ET[i][j]->Scale(1./totalweight);
	  for(int k=0; k<centbins; ++k)
	    {
	      centtow[i][j][k]->Scale(1./totalweight);
	      centet[i][j][k]->Scale(1./totalweight);
	    }
	}
    }
  
  //centet[1][0][17]->Draw();
  //gPad->SaveAs("test.png");
  for(int i=0; i<3; ++i)
    {
      for(int k=0; k<centbins; ++k)
	{
	  meandiff[i]->SetBinContent(centbins-k,(centet[1][i][k]->GetMean()-centet[0][i][k]->GetMean())/(centet[1][i][k]->GetMean()+centet[0][i][k]->GetMean()));
	  for(int j=0; j<2; ++j)
	    {
	      //centet[j][i][k]->Scale(1./centet[j][i][k]->Integral());
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
  outf->mkdir("global");
  outf->mkdir("dETcent");
  outf->mkdir("hcalraw");
  outf->mkdir("fullcor");
  outf->mkdir("dETcentrat");
  outf->mkdir("truthpar");
  outf->mkdir("meanstuff");
  outf->mkdir("centet");
  outf->mkdir("deadmap");
  outf->mkdir("zvtx");
  outf->mkdir("nfill");
  outf->cd("truthpar");
  gDirectory->WriteObject(truthpar_total_ET,truthpar_total_ET->GetName());
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  for(int k=0; k<(j==0?256:64); ++k)
	    {
	      hcalraw[i][j][k]->Divide(nfillcent[i][j][centbins-1]);
	      //if(i==0) hcalraw[i][j][k]->Scale(1./totalweight);
	      outf->cd("hcalraw");
	      gDirectory->WriteObject(hcalraw[i][j][k],hcalraw[i][j][k]->GetName());
	    }
	}
    }
  for(int i=0; i<centbins; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  outf->cd("nfill");
	  gDirectory->WriteObject(nfillcent[0][j][i],nfillcent[0][j][i]->GetName());
	  gDirectory->WriteObject(nfillcent[1][j][i],nfillcent[1][j][i]->GetName());
	  dETcent[0][j][i]->Divide(nfillcent[0][j][i]);
	  dETcent[1][j][i]->Divide(nfillcent[1][j][i]);
	  //dETcent[0][j][i]->Scale(1./totalweight);
	  //dETcent[1][j][i]->Scale(1./totalweight);
	  if(j==0)
	    {
	      dETcent[0][j][i]->Scale(4.);
	      dETcent[1][j][i]->Scale(4.);
	    }
	  //dETcent[0][j][i]->Scale(1./nevtcent[0][i]);
	  //dETcent[1][j][i]->Scale(1./nevtcent[1][i]);
	  if(j==0) truthpar_et[i]->Scale(1./nevtcent[0][i]);
	  fullcor[j][i]->Divide(dETcent[1][j][i],dETcent[0][j][i]);
	  dETcentrat[j][i]->Divide(dETcent[1][j][i],dETcent[0][j][i]);
	  fullcor[j][i]->Multiply(truthpar_et[i]);
	  //fullcor[j][i]->Scale(1./(nevtcent[1][i]));
	  cout <<"bin: " << i << " cal: " << j << " mean dET/deta/(0.5npart): " << fullcor[j][i]->Integral("width")/(2*dETrange*175) <<endl;
	  outf->cd("fullcor");
	  gDirectory->WriteObject(fullcor[j][i],fullcor[j][i]->GetName());
	  if(i==centbins-1) cout << test << " " << dETcent[0][j][i]->Integral("width") << " " << dETcent[1][j][i]->Integral("width") << " " << fullcor[j][i]->Integral("width") << " " << truthpar_et[i]->Integral("width") << endl;
	  dETcentsimunc[j][i]->Scale(1./nevtcent[0][i]);
	  outf->cd("dETcent");
	  gDirectory->WriteObject(dETcentsimunc[j][i],dETcentsimunc[j][i]->GetName());
	}
    }

  cout << nevtcent[0][centbins-1] << " " << nevtcent[1][centbins-1] << endl;
  for(int j=0; j<3; ++j)
    {
      for(int i=0; i<dETbins; ++i)
	{
	  cout << j << " " << i << " " << nfillcent[0][j][centbins-1]->GetBinContent(i+1) << " " << nfillcent[1][j][centbins-1]->GetBinContent(i+1) << endl;
	}
    }
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<3; ++i)
	{
	  for(int j=0; j<centbins; ++j)
	    {
	      outf->cd("dETcentrat");
	      if(h==0) gDirectory->WriteObject(dETcentrat[i][j],dETcentrat[i][j]->GetName());
	      //deadmap[h][i][j]->Divide(deadhits[h][i][j]);
	      deadmap[h][i][j]->Scale(1./nevt[h]);
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
  outf->cd("zvtx");
  gDirectory->WriteObject(zhist[0], zhist[0]->GetName());
  gDirectory->WriteObject(zhist[1], zhist[1]->GetName());
  outf->cd("truthpar");
  gDirectory->WriteObject(truthparehist,truthparehist->GetName());
  gDirectory->WriteObject(truthparnhist,truthparnhist->GetName());
  for(int i=0; i<centbins; ++i)
    {
      gDirectory->WriteObject(truthparecent[i],truthparecent[i]->GetName());
      gDirectory->WriteObject(truthparncent[i],truthparncent[i]->GetName());
      gDirectory->WriteObject(truthpar_et[i],truthpar_et[i]->GetName());
      gDirectory->WriteObject(truthpareetac[i],truthpareetac[i]->GetName());
    }
  outf->cd("meanstuff");
  for(int i=0; i<3; ++i)
    {
      gDirectory->WriteObject(meandiff[i], meandiff[i]->GetName());
    }
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<3; ++i)
	{
	  outf->cd("global");
	  gDirectory->WriteObject(dET[h][i], dET[h][i]->GetName());
	  outf->cd("meanstuff");
	  gDirectory->WriteObject(sigmu[h][i], sigmu[h][i]->GetName());
	  outf->cd("global");
	  gDirectory->WriteObject(ET[h][i], ET[h][i]->GetName());
	  gDirectory->WriteObject(TW[h][i], TW[h][i]->GetName());
	  outf->cd("meanstuff");
	  gDirectory->WriteObject(meancent[h][i],meancent[h][i]->GetName());
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
		  outf->cd("global");
		  gDirectory->WriteObject(ettotcent[h][j], ettotcent[h][j]->GetName());
		  outf->cd("zvtx");
		  gDirectory->WriteObject(zcent[h][j],zcent[h][j]->GetName());
		}
	      outf->cd("deadmap");
	      gDirectory->WriteObject(deadmap[h][i][j], deadmap[h][i][j]->GetName());
	      outf->cd("centet");
	      gDirectory->WriteObject(centtow[h][i][j], centtow[h][i][j]->GetName());
	      gDirectory->WriteObject(centet[h][i][j], centet[h][i][j]->GetName());
	      outf->cd("dETcent");
	      gDirectory->WriteObject(dETcent[h][i][j], dETcent[h][i][j]->GetName());
	    }
	}
    }
  outf->cd("global");
  for(int h=0; h<2; ++h)
    {
      gDirectory->WriteObject(sumev[h], sumev[h]->GetName());
      gDirectory->WriteObject(sumtw[h], sumtw[h]->GetName());
      gDirectory->WriteObject(mbh[h], mbh[h]->GetName());
    }


  cout << "Done writing hists." << endl;
  cout << "Saving parameters to file." << endl;

  
  outpt->Branch("sub",&subtr,"sub/F");
  outpt->Branch("scale",scale,"scale[2]/F");
  outpt->Branch("frac",frac,"frac[2]/I");
  outpt->Branch("mine",&mine,"mine/F");
  int nccb = centbins;
  outpt->Branch("nevtcent",nevtcent,"nevtcent[2][centbins]/I");
  outpt->Branch("nevt",nevt,"nevt[2]/I");
  outpt->Branch("cbin",&nccb,"cbin/I");
  outpt->Branch("zcut",&zcut,"zcut/F");
  outpt->Branch("run",&run,"run/I");
  outpt->Fill();
  outf->WriteObject(outpt,outpt->GetName());
  cout << "Done writing parameters to file" << endl;
  cout << "All done!" << endl;
  /*
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
  */
  return 0;
}
