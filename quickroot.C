#include <TH1D.h>
#include "TROOT.h"
#include "TStyle.h"
#include <cstdio>
#include <TTree.h>
#include "TLatex.h"
#include <cstdlib>
#include <TAxis.h>
#include <TFile.h>
#include "dlUtility.h"

int quickroot()
{
  gStyle->SetOptStat(0);
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
  */
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
  return 0;
}
