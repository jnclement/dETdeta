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

int simtow_check(int nfile)
{
  gStyle->SetOptStat(0);
  int sectors[3] = {24576,1536,1536};
  float sime[2][3][24576];
  float cale[2][3][24576];
  int calet[2][3][24576];
  int calph[2][3][24576];
  TChain* tree[2];
  tree[0] = new TChain("ttree");
  tree[1] = new TChain("ttree");
  for(int i=0; i<nfile; ++i)
    {
      if(i%100 == 0) cout << i << endl;
      try
	{
	  tree[1]->Add(("run/output/evt/events_20231207_shuhang_simtow_nouw_mc_cor_"+to_string(i)+".root").c_str());
	}
      catch(...)
	{
	  continue;
	}
    }
  for(int i=0; i<nfile; ++i)
    {
      if(i%100 == 0) cout << i << endl;
      try
	{
	  tree[0]->Add(("run/output/evt/events_20231207_shuhang_simtow_nomod_mc_cor_"+to_string(i)+".root").c_str());
	}
      catch(...)
	{
	  continue;
	}
    }

  TH2D* calcor[4][3];
  calcor[0][0] = new TH2D("calcor00","",1000,0,0.1,1000,0,5);
  calcor[0][1] = new TH2D("calcor01","",1000,0,0.2,1000,0,2);
  calcor[0][2] = new TH2D("calcor02","",1000,0,0.3,1000,0,10);

  calcor[1][0] = new TH2D("calcor10","",1000,0,0.1,1000,0,0.1);
  calcor[1][1] = new TH2D("calcor11","",1000,0,0.2,1000,0,0.2);
  calcor[1][2] = new TH2D("calcor12","",1000,0,0.3,1000,0,0.3);

  calcor[2][0] = new TH2D("calcor20","",1000,0,0.1,1000,0,5);
  calcor[2][1] = new TH2D("calcor21","",1000,0,0.2,1000,0,2);
  calcor[2][2] = new TH2D("calcor22","",1000,0,0.3,1000,0,10);

  calcor[3][0] = new TH2D("calcor30","",1000,0,5,1000,0,5);
  calcor[3][1] = new TH2D("calcor31","",1000,0,2,1000,0,2);
  calcor[3][2] = new TH2D("calcor32","",1000,0,10,1000,0,10);

  tree[0]->SetBranchAddress("emcalen",cale[0][0]);
  tree[0]->SetBranchAddress("ihcalen",cale[0][1]);
  tree[0]->SetBranchAddress("ohcalen",cale[0][2]);
  tree[0]->SetBranchAddress("emsime",sime[0][0]);
  tree[0]->SetBranchAddress("ihsime",sime[0][1]);
  tree[0]->SetBranchAddress("ohsime",sime[0][2]);
  tree[0]->SetBranchAddress("emcaletabin",calet[0][0]);
  tree[0]->SetBranchAddress("ihcaletabin",calet[0][1]);
  tree[0]->SetBranchAddress("ohcaletabin",calet[0][2]);
  tree[0]->SetBranchAddress("emcalphibin",calph[0][0]);
  tree[0]->SetBranchAddress("ihcalphibin",calph[0][1]);
  tree[0]->SetBranchAddress("ohcalphibin",calph[0][2]);

  tree[1]->SetBranchAddress("emcalen",cale[1][0]);
  tree[1]->SetBranchAddress("ihcalen",cale[1][1]);
  tree[1]->SetBranchAddress("ohcalen",cale[1][2]);
  tree[1]->SetBranchAddress("emsime",sime[1][0]);
  tree[1]->SetBranchAddress("ihsime",sime[1][1]);
  tree[1]->SetBranchAddress("ohsime",sime[1][2]);
  tree[1]->SetBranchAddress("emcaletabin",calet[1][0]);
  tree[1]->SetBranchAddress("ihcaletabin",calet[1][1]);
  tree[1]->SetBranchAddress("ohcaletabin",calet[1][2]);
  tree[1]->SetBranchAddress("emcalphibin",calph[1][0]);
  tree[1]->SetBranchAddress("ihcalphibin",calph[1][1]);
  tree[1]->SetBranchAddress("ohcalphibin",calph[1][2]);
  
  for(int i=0; i<tree[0]->GetEntries(); ++i)
    {
      tree[0]->GetEntry(i);
      tree[1]->GetEntry(i);
      if(i%1000 == 0) cout << "Doing event " << i << endl;
      for(int j=0; j<3; ++j)
	{
	  for(int k=0; k<sectors[j]; ++k)
	    {
	      calcor[0][j]->Fill(sime[0][j][k],cale[0][j][k]);
	      calcor[1][j]->Fill(sime[0][j][k],sime[1][j][k]);
	      calcor[2][j]->Fill(sime[1][j][k],cale[1][j][k]);
	      calcor[3][j]->Fill(cale[0][j][k],cale[1][j][k]);
	    }
	}
    }

  string calname[3] = {"EMCal ","IHCal ","OHCal "};
  string calx[4] = {"Sim Tower Energy from DSTs [GeV]","Sim Tower Energy from DSTs [GeV]","Sim Tower Energies from TowerBuilder [GeV]","Calib Tower Energy from DSTs [GeV]"};
  string caly[4] = {"Calib Tower Energy from DSTs [GeV]","Sim Tower Energy from TowerBuilder [GeV]","Calib Tower Energy from TowerBuilder [GeV]","Calib Tower Energy from TowerBuilder [GeV]"};
  
  TCanvas* c1 =  new TCanvas("c1","c1",1000,1000);
  c1->cd();
  for(int i=0; i<4; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  calcor[i][j]->GetXaxis()->SetTitle((calname[j]+calx[i]).c_str());
	  calcor[i][j]->GetYaxis()->SetTitle((calname[j]+caly[i]).c_str());
	  calcor[i][j]->Draw("COLZ");
	  c1->SaveAs(("simtowcheck/calcor"+to_string(i)+to_string(j)+".png").c_str());
	}
    }
  return 0;
}
