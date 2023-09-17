#include "TROOT.h"
#include <cstdio>
#include <iostream>
#include "TTree.h"
#include "stdlib.h"
#include <fstream>
#include <cstdlib>
#include <TH1D.h>
#include <TFile.h>
#include "mbd_info.h"
#include "TMath.h"

float fill_mbd_dat(int sectors, float* mbe, int* mbt, int* mbs, int* mbc, TH1* hist)
{
  float mbsum, ucmbd;
  int mbdnn, mbdns;
  float mbdtn, mbdts;
  mbdnn=0;
  mbdns=0;
  mbdtn=0;
  mbdts=0;
  mbsum=0;
  ucmbd=0;
  for(int i=0; i<sectors; ++i)
    {
      //cout << mbe[i] << endl;
      ucmbd += mbe[i];
      if(mbt[i] == 1)
	{
	  if(mbs[i] == 1)
	    {
	      //cout << gaincorr[mbc[i]+64] << endl;
	      //cout << i << " " << mbe[i]*gaincorr[mbc[i]+64] << endl;
	      mbsum += mbe[i]*gaincorr[mbc[i]+64];
	    }
	  else if(mbs[i] == 0)
	    {
	      //cout << i << " " << mbe[i]*gaincorr[mbc[i]] << endl;
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
  if(mbsum <=0)
    {
      //cout << "sum <= 0" << endl;
      return -1;
    }
  
  if(isnan(mbdtn-mbdts))
    {
      //cout << "nan time" << endl;
      return -1;
    }
  
  if(abs(mbdtn-mbdts) > 10)
    {
      //cout << "bad zvtx" << endl;
      return -1;
    }
  
  //if(mbsum > 3000) cout << "uh oh!" << endl;
  if(hist) hist->Fill(mbsum);
  return mbsum;
}

int build_test()
{
  mbd_init();
  float mbenrgy[25000];
  int   mbdtype[25000], mbdside[25000], mbdchan[25000];
  int sectormb;
  TFile* file = TFile::Open("/home/jocl/datatemp/merged_dEdeta_71.root");
  TTree* tree;
  tree = file->Get<TTree>("ttree");
  //TFile* simf = TFile::Open("/home/jocl/datatemp/merged_dEdeta_250.root");
  //tree[0] = simf->Get<TTree>("ttree");
  TH1D* mbh;
  cout << tq_t0_offsets[0] << endl;
  cout << gaincorr[0] << endl;
  //mbh[0] = new TH1D("smbh","",500,0,500);
  mbh = new TH1D("dmbh","",3000,0,3000);
  //int npart = 0;
  tree->SetBranchAddress("mbenrgy",mbenrgy);
  tree->SetBranchAddress("mbdtype",mbdtype);
  tree->SetBranchAddress("mbdside",mbdside);
  tree->SetBranchAddress("mbdchan",mbdchan);
  tree->SetBranchAddress("sectormb",&sectormb);
  //tree[0]->SetBranchAddress("npart",&npart);
  int frac[2];
  //frac[0] = 10;
  frac[1] = 100;
  float dummy;
  cout << "Beginning processing." << endl;
  for(int i=0; i<tree->GetEntries()/frac[1]; ++i)
    {
      tree->GetEntry(i);
      /*
      for(int j=0; j<sectormb; ++j)
	{
	  cout << mbenrgy[j] << endl;
	}
      */
      
      dummy = fill_mbd_dat(sectormb, mbenrgy, mbdtype, mbdside, mbdchan, mbh);
      /*
      if(dummy > 3000)
	{
	  cout << i << endl;
	}
      */
    }
  cout << mbh->FindLastBinAbove(0) << endl;
  cout << mbh->GetBinContent(0) << endl;
  cout << mbh->GetBinContent(1) << endl;
  cout << mbh->GetBinContent(mbh->FindLastBinAbove(0)) << endl;
  cout << mbh->GetBinContent(3001) << endl;
  exit(1);
  return 0;
}
