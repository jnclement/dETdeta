#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TMath.h"
#include <utility>  // for std::pair
#include <cstdio>
#include <iostream>
#include "TGraph.h"
#include <iomanip>
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
#include "TProfile.h"
#include "/home/jocl/Documents/main/physics/projects/sphenix_macros/macros/macros/sPHENIXStyle/sPhenixStyle.h"
#include "/home/jocl/Documents/main/physics/projects/sphenix_macros/macros/macros/sPHENIXStyle/sPhenixStyle.C"

float fill_mbd_dat(int sectors, float* mbe, int* mbt, int* mbs, int* mbc, TH1* hist, float zcut, float zval, TH1* zhist)
{
  float mbsum, ucmbd;
  int mbdnn, mbdns;
  float mbdtn, mbdts;
  float zvtx;
  mbdnn=0;
  mbdns=0;
  mbdtn=0;
  mbdts=0;
  mbsum=0;
  ucmbd=0;
  for(int i=0; i<sectors; ++i)
    {
      ucmbd += mbe[i];
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
  if(mbdns == 0 || mbdnn == 0) return -1;
  mbdts/=mbdns;
  mbdtn/=mbdnn;
  zvtx = zval;
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
  if(abs(zvtx) > zcut)
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
  int n=1;
  int nsum=0;
  cout << hist->GetEntries() << endl;
  while(n<centbins)
    {
      //cout << n << endl;
      nsum = 0;
      for(int i=0; i<hist->GetNbinsX()+1; ++i)
	{
	  //cout << i << endl;
	  if(i==hist->GetNbinsX())
	    {
	      cout << "Your centrality bins cannot be set for hist " << hist->GetName() << " centrality bin " << n << " because your hist range is not wide enough. All greater bins will also fail." << endl;
	      return -1;
	      //exit(1);
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
  cent[centbins] = 999999;
  return 0;
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
  return (E-sub)/cosh(-(eta-47.5)*0.024);//*sin(2*atan(exp(-(eta-47.5)*0.024)));
}

float get_E_T_hc(float E, int eta, float sub)
{
  return (E-sub)/cosh(-(eta-11.5)*0.096);//*sin(2*atan(exp(-(eta-11.5)*0.096)));
}


int build_hists(float zcut = 30, float simscale = 1.3, float subtracted = 0, float mine = 0)
{
  mbd_init();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  const int centbins = 9;
  float mbenrgy[25000], calen[2][3][25000];
  int calet[2][3][25000], calph[2][3][25000];
  int   mbdtype[25000], mbdside[25000], mbdchan[25000];
  float towercomb[64][24];
  int sector[2][3];
  int sectormb;
  int npart = 0;
  float z_v[2];
  TFile* hottowers = TFile::Open("/home/jocl/datatemp/hot_towers_21518_1.root");
  TTree* hottree = hottowers->Get<TTree>("T_hot_tower");
  TFile* file = TFile::Open("/home/jocl/datatemp/merged_dEdeta_71.root");
  TTree* tree[2];
  tree[1] = file->Get<TTree>("ttree");
  TFile* simf = TFile::Open("/home/jocl/datatemp/merged_dEdeta_250.root");
  tree[0] = simf->Get<TTree>("ttree");  
  float cents[2][centbins+1] = {0};
  TH1D* centtow[2][3][centbins];
  TH1D* centet[2][3][centbins];
  TH1D* meandiff[3];
  TH1D* sigmu[2][3];
  for(int j=0; j<3; ++j)
    {
      meandiff[j] = new TH1D(("md"+to_string(j)).c_str(),"",9,0,90);
      for(int i=0; i<2; ++i)
	{
	  sigmu[i][j] = new TH1D(("sigmu"+to_string(i)+to_string(j)).c_str(),"",9,0,90);
	}
    }
  float et_em_range[centbins] = {100,150,275,350,400,600,800,1200,1750};
  float et_oh_range[centbins] = {35,50,80,100,140,175,225,300,400};
  float et_ih_range[centbins] = {10,15,25,35,50,75,100,120,150};
  float et_sm_range = 2000;
  float tw_em_range = 15;
  float tw_oh_range = 20;
  float tw_ih_range = 2.5;
  float tw_sm_range = 20;
  int bins_tw = 100;
  int bins_et = 100;
  for(int i=0; i<centbins; ++i)
    {
      centtow[1][0][i] = new TH1D(("centtow10_" + to_string(i)).c_str(),"",bins_tw,0,tw_em_range*(10.+i)/20);
      centtow[1][1][i] = new TH1D(("centtow11_" + to_string(i)).c_str(),"",bins_tw,0,tw_ih_range*(10.+i)/20);
      centtow[1][2][i] = new TH1D(("centtow12_" + to_string(i)).c_str(),"",bins_tw,0,tw_oh_range*(10.+i)/20);

      centtow[0][0][i] = new TH1D(("centtow00_" + to_string(i)).c_str(),"",bins_tw,0,tw_em_range*(10.+i)/20);
      centtow[0][1][i] = new TH1D(("centtow01_" + to_string(i)).c_str(),"",bins_tw,0,tw_ih_range*(10.+i)/20);
      centtow[0][2][i] = new TH1D(("centtow02_" + to_string(i)).c_str(),"",bins_tw,0,tw_oh_range*(10.+i)/20);

      centet[1][0][i] = new TH1D(("centet10_" + to_string(i)).c_str(),"",bins_et,0,et_em_range[i]);
      centet[1][1][i] = new TH1D(("centet11_" + to_string(i)).c_str(),"",bins_et,0,et_ih_range[i]);
      centet[1][2][i] = new TH1D(("centet12_" + to_string(i)).c_str(),"",bins_et,0,et_oh_range[i]);

      centet[0][0][i] = new TH1D(("centet00_" + to_string(i)).c_str(),"",bins_et,0,et_em_range[i]);
      centet[0][1][i] = new TH1D(("centet01_" + to_string(i)).c_str(),"",bins_et,0,et_ih_range[i]);
      centet[0][2][i] = new TH1D(("centet02_" + to_string(i)).c_str(),"",bins_et,0,et_oh_range[i]);
    }

  TH1D* ET[2][3];
  TH1D* TW[2][3];
  TH1D* sumev[2];
  TH1D* sumtw[2];
  TH1D* mbh[2];
  
  mbh[0] = new TH1D("smbh","",500,0,500);
  mbh[1] = new TH1D("dmbh","",3000,0,3000);
  ET[1][0] = new TH1D("et10","",bins_et,0,et_em_range[centbins-1]);
  ET[0][0]  = new TH1D("et00","",bins_et,0,et_em_range[centbins-1]);
  TW[0][0] = new TH1D("tw00","",bins_tw,0,tw_em_range);
  TW[1][0] = new TH1D("tw10","",bins_tw,0,tw_em_range);
  ET[1][1] = new TH1D("et11","",bins_et,0,et_ih_range[centbins-1]);
  ET[0][1] = new TH1D("et01","",bins_et,0,et_ih_range[centbins-1]);
  TW[0][1] = new TH1D("tw01","",bins_tw,0,tw_ih_range);
  TW[1][1] = new TH1D("tw11","",bins_tw,0,tw_ih_range);
  ET[1][2] = new TH1D("et12","",bins_et,0,et_oh_range[centbins-1]);
  ET[0][2] = new TH1D("et02","",bins_et,0,et_oh_range[centbins-1]);
  TW[0][2] = new TH1D("tw02","",bins_tw,0,tw_oh_range);
  TW[1][2] = new TH1D("tw12","",bins_tw,0,tw_oh_range);
  sumev[0] = new TH1D("sumev0","",bins_et,0,et_sm_range);
  sumev[1] = new TH1D("sumev1","",bins_et,0,et_sm_range);
  sumtw[0] = new TH1D("sumtw0","",bins_tw,0,tw_sm_range);
  sumtw[1] = new TH1D("sumtw1","",bins_tw,0,tw_sm_range);
  tree[1]->SetBranchAddress("mbenrgy",mbenrgy);
  tree[1]->SetBranchAddress("emcalen",calen[1][0]);
  tree[1]->SetBranchAddress("emcalet",calet[1][0]);
  tree[1]->SetBranchAddress("emcalph",calph[1][0]);
  tree[1]->SetBranchAddress("ihcalen",calen[1][1]);
  tree[1]->SetBranchAddress("ihcalet",calet[1][1]);
  tree[1]->SetBranchAddress("ihcalph",calph[1][1]);
  tree[1]->SetBranchAddress("ohcalen",calen[1][2]);
  tree[1]->SetBranchAddress("ohcalet",calet[1][2]);
  tree[1]->SetBranchAddress("ohcalph",calph[1][2]);
  tree[1]->SetBranchAddress("mbdtype",mbdtype);
  tree[1]->SetBranchAddress("mbdside",mbdside);
  tree[1]->SetBranchAddress("mbdchan",mbdchan);
  tree[1]->SetBranchAddress("sectorem",&sector[1][0]);
  tree[1]->SetBranchAddress("sectorih",&sector[1][1]);
  tree[1]->SetBranchAddress("sectoroh",&sector[1][2]);
  tree[1]->SetBranchAddress("sectormb",&sectormb);
  tree[1]->SetBranchAddress("zvtx",&z_v[1]);
  tree[0]->SetBranchAddress("zvtx",&z_v[0]);
  tree[0]->SetBranchAddress("npart",&npart);
  tree[0]->SetBranchAddress("emcalen",calen[0][0]);
  tree[0]->SetBranchAddress("emcalet",calet[0][0]);
  tree[0]->SetBranchAddress("emcalph",calph[0][0]);
  tree[0]->SetBranchAddress("ihcalen",calen[0][1]);
  tree[0]->SetBranchAddress("ihcalet",calet[0][1]);
  tree[0]->SetBranchAddress("ihcalph",calph[0][1]);
  tree[0]->SetBranchAddress("ohcalen",calen[0][2]);
  tree[0]->SetBranchAddress("ohcalet",calet[0][2]);
  tree[0]->SetBranchAddress("ohcalph",calph[0][2]);
  tree[0]->SetBranchAddress("sectorem",&sector[0][0]);
  tree[0]->SetBranchAddress("sectorih",&sector[0][1]);
  tree[0]->SetBranchAddress("sectoroh",&sector[0][2]);

  TH1D* zhist = new TH1D("zhist","",200,-100,100);
  float subtr = subtracted;//0.018;
  float scale[2];
  scale[0] = simscale;
  scale[1] = 1;
  int frac[2];
  frac[0] = 1;
  frac[1] = 1;
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
  for(int i=0; i<par; ++i)
    {
      streams[i] << std::fixed << std::setprecision(precision[i]) << parval[i]*mult[i];
      params[i] = streams[i].str();
    }
  string outname = "/home/jocl/datatemp/savedhists_subtr_" + params[1] + "_minE_" + params[2] + "_scale_" + params[0] + "_zcut_" + params[3] + "_run_"+to_string(run)+ ".root";
  TFile* outf = TFile::Open(outname.c_str(),"RECREATE");
  TTree* outt = new TTree("ttree","");
  float dummy;
  float eval;
  float allsum;
  float esum;
  float mbsum;
  int toprint = 10000/frac[0];
  cout << "Events for sim:  " << tree[0]->GetEntries()/frac[0] << endl;
  cout << "Events for data: " << tree[1]->GetEntries()/frac[1] << endl;
  cout << "Beginning processing." << endl;
  cout << "Filling MBD sim hist." << endl;
  for(int i=0; i<tree[0]->GetEntries()/frac[0]; ++i)
    {
      if(i%toprint == 0) cout << "Doing event " << i << endl;
      tree[0]->GetEntry(i);
      if(z_v[0] == 0) continue;
      if(abs(z_v[0]) > zcut) continue;
      if(npart == 0) continue;
      mbh[0]->Fill(npart);
    }
  cout << "Done filling sim MBD hist." << endl;
  cout << "Filling MBD data hist." << endl;
  for(int i=0; i<tree[1]->GetEntries()/frac[1]; ++i)
    {
      if(i%toprint == 0) cout << "Doing event " << i << endl;
      tree[1]->GetEntry(i);
      dummy = fill_mbd_dat(sectormb, mbenrgy, mbdtype, mbdside, mbdchan, mbh[1], zcut, z_v[1], zhist);
    }
  cout << "Done filling data MBD hist." << endl;
  cout << "Done filling all MBD hists." << endl;
  cout << "Setting centrality bins" << endl;
  dummy = set_cent_cuts(mbh[0], cents[0], centbins);
  dummy = set_cent_cuts(mbh[1], cents[1], centbins);
  cout << "Done setting centrality bins." << endl;
  for(int h=0; h<2; ++h)
    {
      cout << "Doing tree " << h << "." << endl;
      for(int i=0; i<tree[h]->GetEntries()/frac[h]; ++i)
	{
	  if(i%toprint==0) cout << "Starting event " << i << endl;
	  set_em_combined_towers_0(towercomb);
	  tree[h]->GetEntry(i);
	  if(h==0)
	    {
	      if(z_v[0] == 0) continue;
	      if(abs(z_v[0]) > zcut) continue;
	      mbsum = npart;
	    }
	  else mbsum = fill_mbd_dat(sectormb, mbenrgy, mbdtype, mbdside, mbdchan, NULL, zcut, z_v[1], NULL);
	  if(mbsum < 0) continue;
	  for(int j=0; j<centbins; ++j)
	  {
	    if(mbsum < cents[h][j+1])
	      {
		for(int k=0; k<3; ++k)
		    {
		      esum = 0;
		      for(int l=0; l<sector[h][k]; ++l)
			{
			  if(calen[h][k][l] < mine) continue;
			  if(k==0)
			    {
			      if(check_acceptance(calet[h][k][l], calph[h][k][l])) continue;
			      eval = scale[h]*get_E_T_em(calen[h][k][l], calet[h][k][l], subtr);
			    }
			  else eval = scale[h]*get_E_T_hc(calen[h][k][l],calet[h][k][l], subtr);
			  esum += eval;
			  TW[h][k]->Fill(eval);
			  centtow[h][k][j]->Fill(eval);
			  if(k==0) towercomb[calph[h][k][l]/4][calet[h][k][l]/4] += eval;
			  else towercomb[calph[h][k][l]][calet[h][k][l]] += eval;
			}
		      allsum += esum;
		      ET[h][k]->Fill(esum);
		      centet[h][k][j]->Fill(esum);
		    }
		  sumev[h]->Fill(allsum);
		  for(int k=0; k<64; ++k)
		    {
		      for(int l=0; l<24; ++l)
			{
			  if(towercomb[k][l] < mine) continue;
			  sumtw[h]->Fill(towercomb[k][l]);
			}
		    }
		  allsum = 0;
		  break;
		  
		  }
	  }
      }
      cout << "Done." << endl;
    }

  for(int i=0; i<3; ++i)
    {
      for(int k=0; k<centbins; ++k)
	{
	  meandiff[i]->SetBinContent(k+1,(centet[1][i][k]->GetMean()-centet[0][i][k]->GetMean())/(centet[1][i][k]->GetMean()+centet[0][i][k]->GetMean()));
	  for(int j=0; j<2; ++j)
	    {
	      sigmu[j][i]->SetBinContent(k+1,centet[j][i][k]->GetStdDev()/centet[j][i][k]->GetMean());
	    }
	}
    }
  
  cout << "Saving hists to " << outname << endl;

  outf->WriteObject(zhist, zhist->GetName());
  for(int i=0; i<3; ++i)
    {
      outf->WriteObject(meandiff[i], meandiff[i]->GetName());
    }
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<3; ++i)
	{
	  outf->WriteObject(sigmu[h][i], sigmu[h][i]->GetName());
	  outf->WriteObject(ET[h][i], ET[h][i]->GetName());
	  outf->WriteObject(TW[h][i], TW[h][i]->GetName());
	}
    }
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<3; ++i)
	{
	  for(int j=0; j<centbins; j++)
	    {
	      outf->WriteObject(centtow[h][i][j], centtow[h][i][j]->GetName());
	      outf->WriteObject(centet[h][i][j], centet[h][i][j]->GetName());
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
  outt->Branch("cbin",&nccb,"cbin/I");
  outt->Branch("zcut",&zcut,"zcut/F");
  outt->Branch("run",&run,"run/I");
  outt->Fill();
  outf->WriteObject(outt,outt->GetName());
  cout << "Done writing parameters to file" << endl;
  cout << "All done!" << endl;
  
  return 0;
}
