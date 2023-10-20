#include <TH1D.h>
#include "TROOT.h"
#include "TStyle.h"
#include <cstdio>
#include <TTree.h>
#include "TLatex.h"
#include <cstdlib>
#include <cmath>
#include <TAxis.h>
#include <TFile.h>
#include "dlUtility.h"

int quickroot()
{
  gStyle->SetOptStat(0);
  /*
  TFile* file = new TFile("/home/jocl/datatemp/savedhists_fracsim_1_fracdat_1_subtr_0_minE_0_scale_1.30_zcut_30_run_21615_20231011_nopileup_cor.root");
  TFile* file2 = new TFile("/home/jocl/datatemp/savedhists_fracsim_1_fracdat_1_subtr_0_minE_0_scale_1.00_zcut_30_run_21615_20231016_nopileup_21615_cor.root");
  //TTree* tree = (TTree*)file->Get("ttree");
  TH1D* phist = new TH1D("phist","",200,0,2);
  TH1D* khist = new TH1D("khist","",200,0,2);
  int npart;
  float pE[100000];
  int pid[100000];
  TH1D* dETcent00_17 = (TH1D*)file->Get("dETcent00_17");
  TH1D* dETcentcount00_17 = (TH1D*)file->Get("dETcentcount00_17");
  TH1D* dETcent10_17 = (TH1D*)file->Get("dETcent10_17");
  TH1D* dETcentcount10_17 = (TH1D*)file->Get("dETcentcount10_17");

  TH1D* dETcent00_172 = (TH1D*)file2->Get("dETcent00_17");
  TH1D* dETcent10_172 = (TH1D*)file2->Get("dETcent10_17");
  cout << dETcent10_17->Integral() << endl;
  //cout << dETcent10_172->Integral() << endl;
  cout << dETcent00_17->Integral() << endl;
  //cout << dETcent00_172->Integral() << endl;
  dETcent00_17->Multiply(dETcentcount00_17);
  dETcent10_17->Multiply(dETcentcount10_17);
  /*
  tree->SetBranchAddress("truthpar_nh",&npart);
  tree->SetBranchAddress("truthparh_e",pE);
  tree->SetBranchAddress("truthparh_id",pid);

  for(int i=0; i<tree->GetEntries(); ++i)
    {
      tree->GetEntry(i);
      for(int j=0; j<npart; ++j)
	{
	  if(pid[j] == 2212)
	    {
	      phist->Fill(pE[j]);
	    }
	  if(pid[j] == 130 || pid[j] == 310 || pid[j] == 311 || pid[j] == 321)
	    {
	      khist->Fill(pE[j]);
	    }
	}
    }
    

  phist->GetXaxis()->SetTitle("Proton Total Energy [GeV]");
  phist->GetYaxis()->SetTitle("Counts");
  khist->GetXaxis()->SetTitle("Combined Kaon Total Energy [GeV]");
  khist->GetYaxis()->SetTitle("Counts");
  TCanvas* c1 = new TCanvas("","");
  c1->cd();
  //phist->Draw();
  dETcent10_17->SetLineColor(kGreen);
  dETcent10_172->SetLineColor(kRed);
  dETcent00_172->SetLineColor(kBlue);
  dETcent10_172->Draw();
  dETcent00_172->Draw("SAME");
  dETcent10_17->Draw("SAME");
  dETcent00_17->Draw("SAME");
  sphenixtext();
  //c1->SaveAs("phist.png");
  //khist->Draw();
  //sphenixtext();
  c1->SaveAs("khist.png");
  */
  TFile* file = new TFile("/home/jocl/events_20231018_nopileup_mc_cor_36.root");
  float eta[100000];
  int npar;
  float E[100000];
  float hepphi[10000];
  float g4phi[10000];
  TTree* tree = (TTree*)file->Get("ttree");
  tree->SetBranchAddress("truthparh_eta",eta);
  tree->SetBranchAddress("truthpar_nh",&npar);
  tree->SetBranchAddress("truthparh_e",E);
  tree->SetBranchAddress("truthparh_phi",hepphi);
  tree->SetBranchAddress("truthpar_phi",g4phi);
  TH2D* hist = new TH2D("hist","",120,-6,6,100,0,10);
  TH1D* dEde = new TH1D("dEde","",20,-1,1);
  TH1D* tote = new TH1D("tote","",1000,0,80000);
  TH2D* phph = new TH2D("phph","",100,-M_PI,M_PI,100,-M_PI,M_PI);
  
  for(int i=0; i<tree->GetEntries(); ++i)
    {
      float total_e = 0;
      tree->GetEntry(i);
      for(int j=0; j<npar; ++j)
	{
	  hist->Fill(eta[j],E[j]/cosh(eta[j]));
	  dEde->Fill(eta[j],E[j]/cosh(eta[j]));
	  phph->Fill(hepphi[j],g4phi[j]);
	  total_e+=E[j];
	}
      tote->Fill(total_e);
    }
  hist->Scale(1./tree->GetEntries());
  dEde->Scale(1./tree->GetEntries());
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->cd();
  hist->Draw("COLZ");
  c1->SaveAs("2d.png");
  dEde->Draw();
  c1->SaveAs("1d.png");
  tote->Draw();
  c1->SaveAs("tote.png");
  phph->Draw("COLZ");
  c1->SaveAs("phph.png");
  return 0;
}
