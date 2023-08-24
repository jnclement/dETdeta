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
int build_events()
{
  mbd_init();
  cout << gaincorr[0] << endl;
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
  float cents[11];
  int sectorem, sectoroh, sectorih, sectormb;
  int ssectorem, ssectoroh, ssectorih, ssectormb;
  TFile* hottowers = TFile::Open("/home/jocl/datatemp/hot_towers_21518_1.root");
  TTree* hottree = hottowers->Get<TTree>("T_hot_tower");
  TFile* file = TFile::Open("/home/jocl/datatemp/merged_dEdeta_19.root");
  TTree* tree = file->Get<TTree>("ttree");
  TFile* simf = TFile::Open("/home/jocl/datatemp/merged_dEdeta_250.root");
  TTree* simt = simf->Get<TTree>("ttree");
  TH1F* hist = new TH1F("hist","",394,0,394);
  TH1F* centclass[10][3];
  TH1F* etdata = new TH1F("etdata","",100,0,500);
  TH1F* etsim  = new TH1F("etsim","",100,0,500);
  TH1F* sedist = new TH1F("sedist","",1000,0,50);
  TH1F* dedist = new TH1F("dedist","",1000,0,50);
  TH1F* mthist = new TH1F("mthist","",100,-50,50);
  for(int i=0; i<10; ++i)
    {
      centclass[i][0] = new TH1F(("centhiste" + to_string(i)).c_str(),"",96,-1.1,1.1);
      centclass[i][1] = new TH1F(("centhisti" + to_string(i)).c_str(),"",24,-1.1,1.1);
      centclass[i][2] = new TH1F(("centhisto" + to_string(i)).c_str(),"",24,-1.1,1.1);
    }
  float detacente[10][96] = {0};
  float detacento[10][24] = {0};
  float detacenti[10][24] = {0};
  int nem[10][96] = {0};
  int nih[10][24] = {0};
  int noh[10][24] = {0};
  int tpn = 0;
  float tpe[25000] = {0};
  float tpet[25000] = {0};
  float tpph[25000] = {0};
  int tpem[25000] = {0};
  int npart = 0;
  int ncoll = 0;
  float bimp = 0;
  TH2F* ectee = new TH2F("ectee","",100,0,2000,100,0,3000);
  TH2F* ectea = new TH2F("ectea","",100,0,2000,100,0,3000);
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
  int fractouse = 10;
  int counter = 0;
  int mbdsum = 0;
  int mbdsumunc = 0;
  //MBD timing channels 54 and 56 are bad
  //cout << tree->GetEntries()/fractouse << endl;
  for(int i=0; i<simt->GetEntries()/fractouse; ++i)
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

  
  
  int n = 0;
  int nsum = 0;
  cout << "test" << endl;
  while(n<10)
    {
      nsum = 0;
      for(int i=1; i<hist->GetNbinsX()+1; ++i)
	{
	  //cout <<n <<" "<<i << " " << nsum << endl;
	  nsum += hist->GetBinContent(i);
	  if(n*hist->GetEntries()/10 < nsum)
	    {
	      lines[n] = new TLine(hist->GetBinCenter(i),0,hist->GetBinCenter(i), hist->GetBinContent(i));
	      lines[n]->SetLineColor(kRed);
	      cents[n] = hist->GetBinCenter(i);
	      cout << cents[n] << endl;
	      ++n;
	      break;
	    }
	}
    }
  cents[10] = 99999999999;
  cout << "test4" << endl;
  int kcodes[12] = {0,1,2,3,4,-1,-5,-8,-10,+1,+2,+3};
  for(int i=0; i<10; ++i)
    {
      setcolorcent(centclass[i][0], kcodes, kRed, 10, i);
      setcolorcent(centclass[i][1], kcodes, kGreen, 10, i);
      setcolorcent(centclass[i][2], kcodes, kBlue, 10, i);
    }

  for(int i=0; i<simt->GetEntries()/fractouse; ++i)
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











    
  for(int i=0; i<simt->GetEntries()/fractouse; ++i)
    {
      mbdsum = 0;
      float sime = 0;
      simt->GetEntry(i);
      mbdsum = npart;
      //cout << mbdsum << endl;
      for(int j=0; j<10; ++j)
	{
	  if(mbdsum < cents[j+1])
	    {
	      //cout <<sectorem << endl;
	      //cout << j << endl;
	      int etabin;
	      for(int k=0; k<ssectorem; ++k)
		{
		  etabin = semcalet[k];
		  sedist->Fill(semcalen[k]*sin(2*atan(exp(-(etabin-48)*0.024))));
		  if(semcalen[k] < 0) continue;
		  //floor(emcalet[k]/0.023)+48;
		  detacente[j][etabin] += semcalen[k]*sin(2*atan(exp(-(etabin-48)*0.024)));//-emcalet[k])));
		  sime += semcalen[k]*sin(2*atan(exp(-(etabin-48)*0.024)));
		  //cout << j << "etabin" << etabin << " detacente " << detacente[j][etabin] << endl;
		  //cout << emcalen[k] << endl;
		  if(semcalen[k] > 0)
		    {
		      nem[j][etabin]++;
		    }
		}
	      for(int k=0; k<ssectorih; ++k)
		{
		  if(sihcalen[k] < 0) continue;
		  etabin = sihcalet[k];//floor(ihcalet[k]/0.092)+12;
		  detacenti[j][etabin] += sihcalen[k]*sin(2*atan(exp(-(etabin-12)*0.096)));//-ihcalet[k])));
		  //sime += 1.3*sihcalen[k]*sin(2*atan(exp(-(etabin-12)*0.096)));
		  nih[j][etabin]++;
		}
	      for(int k=0; k<ssectoroh; ++k)
		{
		  if(sohcalen[k] < 0) continue;
		  etabin = sohcalet[k];//floor(ohcalet[k]/0.092)+12;
		  detacento[j][etabin] += sohcalen[k]*sin(2*atan(exp(-(etabin-12)*0.096)));//-ohcalet[k])));
		  if(sohcalen[k] > 0)
		    {
		      noh[j][etabin]++;
		    }
		}
	      etsim->Fill(sime);
	      //cout << detacenti[j][12] << endl;
	      break;
	    }
	}
    }
  cout << detacenti[5][12] << endl;
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
  for(int i=0; i<10; ++i)
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
  drawText("#bf{sPHENIX} internal", 0.85, 0.675, 1, kBlack, 0.04);
  //drawText("Simulated reco MBD charge sum", 0.85,0.825,1, kBlack,0.04);
  //drawText("Red lines represent 10\%-ile bins in",0.85,0.775,1,kBlack,0.04);
  //drawText("centrality, with most central at right",0.85,0.725,1,kBlack,0.04);
  for(int i=0; i<10; ++i) lines[i]->Draw();
  c1->SaveAs("cent_edep.png");
  TCanvas* c2 = new TCanvas("c2","c2",1000,1000);
  c2->SetLogy();
  c2->SetTicks(1);
  for(int j=0; j<3; ++j)
    {
      float min = 999999999;
      float max = 0;
      for(int i=0; i<10; ++i)
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
      for(int i=0; i<10; ++i)
	{
	  for(int k=0; k<8; ++k)
	    {
	      centclass[i][0]->SetBinContent(k+1,0);
	    }
	}
      for(int i=0; i<10; ++i)
	{
	  c2->cd();
	  if(i==0)
	    {
	      centclass[i][j]->GetYaxis()->SetRangeUser(min/2,max*2);
	      centclass[i][j]->Draw();
	    }
	  else centclass[i][j]->Draw("SAME");
	}
      c2->SaveAs(("centdETdeta"+to_string(j)+".pdf").c_str());
    }

  float date = 0;
  for(int i=0; i<tree->GetEntries()/fractouse; ++i)
    {
      float zvtx = 0;
      mbdsum = 0;
      float mbdtn = 0;
      float mbdts = 0;
      int mbdnn = 0;
      int mbdns = 0;
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
      if(isnan(mbdtn-mbdts)) continue;
      if(abs(mbdtn-mbdts)>10) continue;
      else
	{
	  mthist->Fill(mbdtn-mbdts);
	  //cout << mbdtn - mbdts << endl;
	}
      //cout << mbdsum << endl;
      //for(int j=0; j<10; ++j)
	{
	  //if(mbdsum < cents[j+1])
	    {
	      //cout <<sectorem << endl;
	      //cout << j << endl;
	      int etabin;
	      float subtr = 0.018;
	      for(int k=0; k<sectorem; ++k)
		{
		  etabin = emcalet[k];
		  dedist->Fill((emcalen[k]-subtr)*sin(2*atan(exp(-(etabin-48)*0.024))));
		  if(emcalen[k] < 0) continue;
		  //floor(emcalet[k]/0.023)+48;
		  //detacente[j][etabin] += emcalen[k]*sin(2*atan(exp(-(etabin-48)*0.024)));//-emcalet[k])));
		  date += (emcalen[k]-subtr)*sin(2*atan(exp(-(etabin-48)*0.024)));//-emcalet[k])));
		  //cout << j << "etabin" << etabin << " detacente " << detacente[j][etabin] << endl;
		  //cout << emcalen[k] << endl;
		  //if(emcalen[k] > 0)
		  //{
		  //nem[j][etabin]++;
		  //}
		}
	      for(int k=0; k<sectorih; ++k)
		{
		  etabin = ihcalet[k];//floor(ihcalet[k]/0.092)+12;
		  //detacenti[j][etabin] += ihcalen[k]*sin(2*atan(exp(-(etabin-12)*0.096)));//-ihcalet[k])));
		  //date += ihcalen[k]*sin(2*atan(exp(-(etabin-12)*0.096)));
		  //if(ihcalen[k] > 0)
		  //{
		  //nih[j][etabin]++;
		  //}
		}
	      for(int k=0; k<sectoroh; ++k)
		{
		  etabin = ohcalet[k];//floor(ohcalet[k]/0.092)+12;
		  //detacento[j][etabin] += ohcalen[k]*sin(2*atan(exp(-(etabin-12)*0.096)));//-ohcalet[k])));
		  //if(ohcalen[k] > 0)
		  //{
		  //noh[j][etabin]++;
		  //}
		}
	      etdata->Fill(date);
	      //break;
	    }
	}
	
    }
  c2->cd();
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
  etdata->GetXaxis()->SetTitle("E_{T,total IHCal} [GeV]");
  etdata->Draw();
  etsim->Draw("SAME");
  drawText("#bf{sPHENIX} internal", 0.85, 0.93, 1, kBlack, 0.04);
  drawText("Red: data",0.85,0.79,1,kBlack,0.025);
  drawText("Blue: sim",0.85,0.76,1,kBlack,0.025);
  drawText("Sim X scaled 1.3",0.85,0.82,1,kBlack,0.025);
  drawText("Area normed to 1",0.85,0.85,1,kBlack,0.025);
  c2->SaveAs("simanddata.pdf");
  dedist->Scale(1./dedist->Integral());
  sedist->Scale(1./sedist->Integral());
  dedist->SetLineColor(kRed);
  sedist->SetLineColor(kBlue);
  maxval = max(dedist->GetMaximum(),sedist->GetMaximum());
  minval = min(dedist->GetBinContent(dedist->FindLastBinAbove(0,1)),sedist->GetBinContent(sedist->FindLastBinAbove(0,1)));
  dedist->GetXaxis()->SetRangeUser(0,10);
  dedist->GetYaxis()->SetRangeUser(minval/2.,maxval*2.);
  dedist->Draw();
  sedist->Draw("SAME");
  c2->SaveAs("etdists.pdf");
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
  drawText("#bf{sPHENIX} internal", 0.85, 0.675, 1, kBlack, 0.04);
  //drawText("Simulated reco MBD charge sum", 0.85,0.825,1, kBlack,0.04);
  //drawText("Red lines represent 10\%-ile bins in",0.85,0.775,1,kBlack,0.04);
  //drawText("centrality, with most central at right",0.85,0.725,1,kBlack,0.04);
  for(int i=0; i<10; ++i) lines[i]->Draw();
  c1->SaveAs("cent_edep.png");
  //TCanvas* c2 = new TCanvas("c2","c2",1000,1000);
  c2->SetLogy();
  c2->SetTicks(1);
  for(int i=0; i<10; ++i)
    {
      c2->cd();
      centclass[i][0]->Draw();
      centclass[i][1]->Draw("SAME");
      centclass[i][2]->Draw("SAME");
      c2->SaveAs(("centdETdeta_data_"+to_string(i)+".pdf").c_str());
    }
  TCanvas* c3 = new TCanvas("c3","c3",1000,500);
  c3->Divide(2,1);
  c3->cd(1);
  ectee->Draw("COLZ");
  c3->cd(2);
  ectea->Draw("COLZ");
  c3->SaveAs("emcalemvsall.pdf");
  TCanvas* c4 = new TCanvas("c4","c4");
  c4->cd();
  mthist->Draw();
  c4->SaveAs("mthist.pdf");
  return 0;
}
