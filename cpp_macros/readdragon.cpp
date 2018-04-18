#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TProfile.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TLatex.h>

using namespace std;

#include "event.h"
#include "ffactor.h"
#include "pedestal.h"
#include "pixstat.h"
#include "timecal.h"
#include <TTimer.h>
#include <Getline.h>

Bool_t HandleInput()
{
    TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
    while (1){
        timer.TurnOn();TString input = Getline("Type 'q' to exit, <return> to go on: ");timer.TurnOff();

        if (input=="q\n") return kFALSE;
        if (input=="\n")  return kTRUE;
    };
    return kFALSE;
}

void FillHist(Event &evt, TH1 *h, int pix=0, int gain=0)
{
  for (int i=0; i<evt.roisize; i++)
    h->SetBinContent(i+1, evt.samples[pix][gain][i]);
}



//==================================================================
// rate - NSB rate [slices^-1]
// pulsewidth - RMS of the pulse shape [slices]
// onephemean - average integral 1 phe signal
// onepherms - rms of integral 1 phe signal
void AddNSB(Event &ev, float rate,  int pixid, int gain, float pulsewidth, float onephemean, float onepherms)
{
  Int_t nnsb=gRandom->Poisson(rate*ev.roisize);
  // cout<<"simulating "<<nnsb<<" NSB photons"<<endl;
  Double_t *poss = new Double_t [nnsb];
  Double_t *amps = new Double_t [nnsb];
  Double_t zeroing=rate*ev.roisize*onephemean/ev.roisize;
  for (int i=0; i<nnsb; i++)
    {
      poss[i]=gRandom->Uniform(0, ev.roisize);
      amps[i]=gRandom->Gaus(onephemean, onepherms);
      if (amps[i]<0)
	amps[i]=0;
    }

  for (int slice=0; slice<ev.roisize; slice++)
    {
      for (int i=0; i<nnsb; i++)
	ev.samples[pixid][gain][slice]+=amps[i]*TMath::Gaus(slice, poss[i], pulsewidth, kTRUE);
      ev.samples[pixid][gain][slice]-=zeroing;
    }
  delete[]poss;
  delete[]amps;

}

void checktime(ifstream &plik, int roi)
{
  Event evt(roi);
  double oldtime133=-1, oldtime10=-1;
  double time133=-1, time10=-1;
  TGraph *gr = new TGraph();
  TH1F *htime = new TH1F ("htime", "", 100, -3, 3);
  for (int ev=0; ev<400000; ev++)
    {
      evt.read(plik);
      oldtime133=time133;
      oldtime10=time10;

      time133=evt.GetDRSTime();
      time10=evt.GetDRSTime2();
      if (oldtime133>0)
	{
	  double timediff133=time133-oldtime133;
	  double timediff10=time10-oldtime10;
	  if (fabs(timediff10-timediff133)>20.e-3)
	    cout<<ev<<" "<<timediff133<<" "<<timediff10<<endl;

	  // if (timediff10>0)
	  gr->SetPoint(gr->GetN(), timediff133, timediff10);
	  htime->Fill(log10(timediff133));
	}
    }
  TCanvas *ctime = new TCanvas("ctime", "", 640, 480);
  ctime->SetLogy();
  // gr->Draw("AP");
  htime->Draw();
}

PedStat checkped(ifstream &plik, Event &evt, PedestalSimple &ped, int window, string tag="", TH1D *hh = 0)
// void checkped(ifstream &plik, Event &evt, PedestalUltimate &ped, int window)
{
  PedStat pedstat;
  TH1F *h[7][2];
  for (int i=0; i<7; i++)
    for (int j=0; j<2; j++)
      {
	h[i][j] = new TH1F(Form("h_pix%i_%s_gain", i, (j==0)?"high":"lo"), "", 200, -50, 150);
	h[i][j]->GetXaxis()->SetTitle(Form("sum of %i slices / sqrt(%i)", window, window));
      }

  for (int ev=0; ev<14999; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt.read(plik))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      removePed(evt, ped);
      evt.CorrTime();
      evt.InterpolatePseudoPulses();
      for (int i=0; i<7; i++)
	for (int j=0; j<2; j++)
	  for (int k=2; k<evt.roisize-2-window; k++)  // skipped the first and the last slice
	    h[i][j]->Fill(evt.SumSlices(i, j, k, window)/sqrt(1.*window));
    }
  TCanvas *can = new TCanvas ("can", "", 1600, 800);
  can->Divide(7, 2, 0.001, 0.001);
  for (int i=0; i<7; i++)
    for (int j=0; j<2; j++)
      {
	can->cd(i+7*j+1);
	can->GetPad(i+7*j+1)->SetLogy();
	h[i][j]->Draw();
	pedstat.fmean[i][j]=h[i][j]->GetMean();
	pedstat.frms[i][j]=h[i][j]->GetRMS();
      }
  // can->Print(Form("%s_%islices_corr.pdf", tag.data(), window));
  for (int i=0; i<7; i++)
    for (int j=0; j<2; j++)
      {
	if (hh)
	  hh->Add(h[i][j]);
	delete h[i][j];
      }
  return pedstat;
}

// function for making pedestal distribution plots
PedStat makepedplot(string namein, int window=1, int roi=40, TH1D *hh = 0)
{
  cout<<"opening "<<namein<<endl;
  ifstream plik(namein.data(), ios::binary);
  if (!plik.is_open())
    {
      cout<<"file: "<<namein<<" is not open, exiting"<<endl;
      PedStat dm;
      return dm;
    }
  Event evt(roi);

  PedestalSimple pedsim;
  // PedestalUltimate pedult(evt.roisize);
  for (int ev=0; ev<10000; ev++)
    {

      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt.read(plik))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      evt.CorrTime();
      pedsim.fillPedEvent(evt);
      // pedult.fillPedEvent(evt);
    }
  pedsim.finalizePed();
  // pedult.finalizePed();

  TString tag(namein);
  tag.ReplaceAll(".dat", "");

  return checkped(plik, evt, pedsim, window, (string)tag, hh);
}

// first get pedestal from first file, than show event one by one from the second
void showevents(string nameped, string namein, int roi, int mypix, int mygain)
{
  cout<<"opening pedestal"<<nameped<<endl;
  ifstream plik1(nameped.data(), ios::binary);
  if (!plik1.is_open())
    {
      cout<<"file: "<<nameped<<" is not open, exiting"<<endl;
      return;
    }
  Event evt1(roi);

  PedestalSimple ped;
  // PedestalUltimate ped(roi);
  for (int ev=0; ev<50000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt1.read(plik1))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      evt1.CorrTime();
      ped.fillPedEvent(evt1);
    }
  plik1.close();
  ped.finalizePed();

  TH1F *hh = new TH1F ("hh", ";RMS in ROI; Number of events", 400, -20, 20);

  ifstream plik2(namein.data(), ios::binary);
  if (!plik2.is_open())
    {
      cout<<"file: "<<namein<<" is not open, exiting"<<endl;
      return;
    }
  Event evt2(roi);

  TH1F *hwf = new TH1F ("hwf", ";sample;signal[counts]", roi, -0.5, roi-0.5);
  TH1F *hwf2 = new TH1F ("hwf2", ";sample;signal[counts]", roi, -0.5, roi-0.5);
  hwf2->SetLineColor(kRed);
  TCanvas *cwf = new TCanvas ("cwf", "", 640, 480);
  cwf->cd();
  hwf->Draw();
  hwf->SetMaximum(150);
  hwf->SetMinimum(-50);
  // hwf2->Draw("same");

  for (int ev=0; ev<20000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt2.read(plik2))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      int fc = evt2.firstcap[mypix][mygain];
      // for (int i=2; i<roi-2; i++)
      // 	cout<<evt2.GetDRSTime() - evt2.lasttime[mypix][mygain][(i+fc) % size4drs]<<", ";
      // cout<<endl;


      FillHist(evt2, hwf2, mypix, mygain);
      removePed(evt2, ped);
      evt2.CorrTime();
      evt2.InterpolatePseudoPulses();
      FillHist(evt2, hwf, mypix, mygain);

      // AddNSB(Event &ev, float rate,  int pixid, int gain, float pulsewidth, float onephemean, float onepherms)
      // AddNSB(evt2, 0.2, mypix, mygain, 1.2, 90, 35.);

      TH1F *hnoise = new TH1F("hnoise", "", 100, -100, 100);

      Float_t max = -1.e7;
      Int_t imax = -10000;
      for (int i=2; i<roi-2; i++)
	{
	  hnoise->Fill(evt2.samples[mypix][mygain][i]);
	  if (evt2.samples[mypix][mygain][i]>max)
	    {
	      max=evt2.samples[mypix][mygain][i];
	      imax = (i+fc+size4drs) % size4drs;
	    }
	  // cout<<evt2.samples[mypix][mygain][i]<<", ";
	}
      cout<<endl;
      cout<<"event "<<ev<<", max = "<<max<<" at "<<imax<<" %32 = "<<imax%32<<endl;
      cout<<"mean = "<<hnoise->GetMean()<<", RMS = "<<hnoise->GetRMS()<<endl;
      hh->Fill(hnoise->GetMean());
      delete hnoise;

      cwf->Modified();
      cwf->Update();
      if (!HandleInput())
      	{
      	  plik2.close();
      	  break;
      	}
    }
  plik2.close();

  cwf->cd();
  hh->Draw();
}


void testnoise(string nameped, string namein, int roi, int mypix, int mygain)
{
  cout<<"opening pedestal"<<nameped<<endl;
  ifstream plik1(nameped.data(), ios::binary);
  if (!plik1.is_open())
    {
      cout<<"file: "<<nameped<<" is not open, exiting"<<endl;
      return;
    }
  Event evt1(roi);

  PedestalSimple ped;
  for (int ev=0; ev<50000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt1.read(plik1))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}
      evt1.CorrTime();
      ped.fillPedEvent(evt1);
    }
  plik1.close();
  ped.finalizePed();

  ifstream plik2(namein.data(), ios::binary);
  if (!plik2.is_open())
    {
      cout<<"file: "<<namein<<" is not open, exiting"<<endl;
      return;
    }
  Event evt2(roi);

  // TH1F *h = new TH1F ("hwf", ";sample;signal[counts]", size4drs, -0.5, size4drs-0.5);
  TH1F *h = new TH1F ("h", ";signal[counts];", 200, -100, 100);
  TH1F *h1 = new TH1F ("h1", ";signal[counts];", 200, -100, 100);

  for (int ev=0; ev<40000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt2.read(plik2))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      removePed(evt2, ped);
      evt2.CorrTime();
      evt2.InterpolatePseudoPulses();

      if (ev<1000)
	continue;

      for (mypix = 0; mypix<7; mypix++)
	for (mygain =0; mygain<2; mygain++)
	  {

      int fc = evt2.firstcap[mypix][mygain];
      for (int i=10; i<roi-2; i++)
	// if ((i+fc) % size4drs % 32 == 31)
	  // h->Fill(evt2.samples[mypix][mygain][i]);
	// else
	   h1->Fill(evt2.samples[mypix][mygain][i]);
	// if (evt2.samples[mypix][mygain][i]>40)
	//   h->Fill(i);
	  // h->Fill((i+fc) % size4drs %64);
	  }
    }
  plik2.close();

  h->SetLineColor(kRed);
  TCanvas *cwf = new TCanvas ("cwf", "", 640, 480);
  cwf->cd();
  h1->Draw();
  // h->Draw("sames");

}

void prephist(TH1F *hist, int n=50)
{
  while ((hist->GetMaximum()<n) && (hist->GetNbinsX()%2 == 0))
    hist->Rebin(2);
  double rms = hist->GetRMS();
  double mean = hist->GetMean();
  hist->GetXaxis()->SetRangeUser(mean - 4 * rms, mean + 4 * rms);
}

bool fitHist(TH1F *hist, FFactorStat &fstat)
{
  double rms = hist->GetRMS();
  double mean = hist->GetMean();

  TF1 *fgaus = new TF1("fgaus", "gaus", mean-4*rms, mean+4*rms);
  hist->Fit(fgaus, "RQLL");
  TF1 *f = hist->GetFunction("fgaus");
  if (f)
    {
      fgaus->SetParameters(50, mean, rms); // const, mean, rms
      double prob = f->GetProb();
      if (prob<1.e-3) cout<<"prob = "<<prob<<endl;
      if (prob>1.e-5)
	{
	  fstat.meansig = f->GetParameter(1);
	  fstat.rmssig = f->GetParameter(2);
	  fstat.doFFactor();
	  fstat.print();
	  return true;
	}
    }
  fstat.isfit = 0;
  fstat.isgood = 0;
  return false;
}

void ffactor(string nameped, string namein, int roi, FFactorStat *fstat=0, int mygain=0, bool batch=false, int windowsize=6)
{
  const int nn=nch-1;
  if (fstat==0) fstat = new FFactorStat[nn];
  PedestalSimple &pedsim=*prepPedSimple(nameped, roi, 50000);

  PedStat *pedstat = getrms(nameped, roi, &pedsim, windowsize);
  pedstat->print();
  for (int i=0; i<nn; i++)
    fstat[i].pedrms=pedstat->frms[i][mygain];

  ifstream plik2(namein.data(), ios::binary);
  if (!plik2.is_open())
    {
      cout<<"file: "<<namein<<" is not open, exiting"<<endl;
      return;
    }
  Event evt2(roi);
  TH1F *hsig[nn];
  for (int i=0; i<nn; i++)
    hsig[i]=new TH1F (Form("hsig%i", i), ";signal[counts]", 4096, 300, 100000);
  TH1F *hsig1[nn];
  for (int i=0; i<nn; i++)
    hsig1[i]=new TH1F (Form("hsig1_%i", i), ";peak signal[counts]", 4096, 100, 10000);
  TH1F *hsigmean=new TH1F ("hsigmean", ";signal[counts]", 4096, 300, 100000);

  for (int ev=0; ev<1000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt2.read(plik2))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      removePed(evt2, pedsim);
      evt2.CorrTime();
      evt2.InterpolatePseudoPulses();
      float sumav=0;
      for (int i=0; i<nn; i++)
	{
	  float sum = evt2.Sliding(i, mygain, windowsize);
	  hsig[i]->Fill(sum);
	  sumav+=sum;
	  float peak = evt2.Sliding(i, mygain, 1);
	  hsig1[i]->Fill(peak);
	}
      sumav/=7;
      hsigmean->Fill(sumav);
      // cout<<"event "<<ev<<", sum = "<<sum<<endl;
    }
  prephist(hsigmean);
  for (int i=0; i<nn; i++)
    {
      prephist(hsig[i]);
      fitHist(hsig[i], fstat[i]);
      fstat[i].peak=hsig1[i]->GetMean()/fstat[i].nphe;
    }

  if (batch)
    {
      for (int i=0; i<nn; i++)
	{
	  delete hsig[i];
	  delete hsig1[i];
	}
      delete hsigmean;
    }
  else
    {
      TCanvas *csig = new TCanvas ("csig", "", 640, 480);
      csig->Divide(4,2,0.001, 0.001);
      for (int i=0; i<nn; i++)
	{
	  csig->cd(i+1);
	  hsig[i]->Draw();
	  TLatex *label = new TLatex(0.14, 0.8, Form("%.1f phe", fstat[i].nphe));
	  TLatex *label2 = new TLatex(0.14, 0.9, Form("1phe=~%.1f cnts", 1./fstat[i].conv));
	  label->SetTextSize(0.07);
	  label->SetNDC();
	  label->Draw();
	  label2->SetNDC();
	  label2->Draw();
	  label2->SetTextSize(0.07);
	}
      csig->cd(8);
      hsigmean->Draw();
    }
}

Double_t funcgain (Double_t *xx, Double_t *par)
{
  Double_t x=*xx;
  Double_t ped0=par[0];
  Double_t pedrms=par[1];
  Double_t nphe=par[2];
  Double_t pherms=par[3];
  Double_t conv = par[4];
  Double_t norm=par[5];
  Double_t suma=TMath::Gaus(x, ped0, pedrms, 1)*TMath::Poisson(0, nphe);
  for (int i=1; i<=5; i++)
    {
      Double_t mean=sqrt(ped0*ped0+i*conv*i*conv);
      Double_t rms = sqrt(pedrms*pedrms+i*pherms*pherms);
      suma+=TMath::Gaus(x, mean, rms,1)*TMath::Poisson(i, nphe);
    }
  return suma*norm;
}

// stacking up of single phe
// nameped - file to get the pedestal calibration of each slice
// namein - file with single phe pulses
// namecal - file with calibration pulses (for DRS4 time correction)
// nameped2 - if not empty a global correction of pedestal shift is taken from this one
void stack1phe(string nameped, string nameped2, string namein, int roi=40, int mygain=0)
{
  const int nn=nch-1;
  const int windowsize=6;
  const int timewindow=6;

  const float minsig=(90-35);
  const float maxsig=(90+35);
  PedestalSimple &pedsim=*prepPedSimple(nameped, roi, 10000);

  // int nfirst=2, nlast=35;
  int nfirst=15-3, nlast=26+3;

  TH1F *htime[nn];
  TH1F *htimecorr[nn];
  TH1F *hshapesw[nn]; // shape when single phe are aligned with sliding window
  TH1F *hshapenosw[nn]; // shape where single phe are NOT aligned with sliding window
  for (int i=0; i<nn; i++)
    {
      htime[i]=new TH1F (Form("htime%i", i), ";arrival time[slices]", 2*roi, -0.5, roi-0.5);
      htimecorr[i]=new TH1F (Form("htimecorr%i", i), ";arrival time[slices]", 2*roi, -0.5, roi-0.5);
      htimecorr[i]->SetLineColor(kRed);

      hshapesw[i]=new TH1F (Form("hshapesw%i", i), ";arrival time[slices]", 60, -10, 10);
      hshapenosw[i]=new TH1F (Form("hshapenosw%i", i), ";arrival time[slices]", 60, -10, 10);
      hshapesw[i]->SetLineColor(kGreen);
    }

  double meanped[nch]={0};
  if (nameped2!="")
    {
      cout<<" getting the pedestal correction from "<<nameped2<<endl;
      TH1F *hpeds[nch];
      for (int i=0;i<nch; i++)
	hpeds[i]=new TH1F (Form("hpeds_%i", i), "", 400, -100, 100);

      // loop over the second (newer) pedestal file
      ifstream plik1b(nameped2.data(), ios::binary);
      if (!plik1b.is_open())
	{
	  cout<<"file: "<<nameped2<<" is not open, exiting"<<endl;
	  return;
	}
      Event evt1b(roi);
      for (int ev=0; ev<10000; ev++)
	{
	  if (ev % 10000 == 0)
	    cout<<"Event "<<ev<<endl;
	  if (!evt1b.read(plik1b))
	    {
	      cout<<"event "<<ev<<"end of file"<<endl;
	      break;
	    }

	  removePed(evt1b, pedsim);
	  evt1b.CorrTime();
	  evt1b.InterpolatePseudoPulses();
	  for (int i=0; i<nn; i++)
	    for (int j=nfirst; j<=nlast; j++)
	      hpeds[i]->Fill(evt1b.samples[i][mygain][j]);
	}
      plik1b.close();
      for (int i=0; i<nn; i++)
	{
	  meanped[i]=hpeds[i]->GetMean();
	  cout<<"Ch "<<i<<", mean ped="<<meanped[i]<<", RMS="<<hpeds[i]->GetRMS()<<endl;
	}
    }

  TimeCalClu caltime(1024, 8, 7);
  // caltime.FillFile(namecal, pedsim, roi, 6);
  caltime.FillFile("../20170519_160606_MONI/MONI_Laser_NominalV_Filter7_Seq00_RD40_FEB1_IP346.dat", pedsim, roi, 6);
  caltime.FillFile("../DATA_Linearity_20170519_132355/Linearity_NominalHV_Filter5RD40_FEB1_IP346.dat", pedsim, roi, 6);
  caltime.FillFile("../DATA_Linearity_20170519_132355/Linearity_NominalHV_Filter6RD40_FEB1_IP346.dat", pedsim, roi, 6);
  caltime.FillFile("../DATA_Linearity_20170519_132355/Linearity_NominalHV_Filter7RD40_FEB1_IP346.dat", pedsim, roi, 6);
  caltime.FillFile("../DATA_Linearity_20170519_132355/Linearity_NominalHV_Filter8RD40_FEB1_IP346.dat", pedsim, roi, 6);

  caltime.Finalize();

  ifstream plik2(namein.data(), ios::binary);
  if (!plik2.is_open())
    {
      cout<<"file: "<<namein<<" is not open, exiting"<<endl;
      return;
    }
  Event evt2(roi);
  int numev[nn]={0};
  for (int ev=0; ev<10000; ev++)
    {
      if (ev % 10000 == 0)
  	cout<<"Event "<<ev<<endl;
      if (!evt2.read(plik2))
  	{
  	  cout<<"event "<<ev<<"end of file"<<endl;
  	  break;
  	}

      removePed(evt2, pedsim);
      evt2.CorrTime();
      evt2.InterpolatePseudoPulses();
      for (int i=0; i<nn; i++)
  	{
	  int n1=int(caltime.GetCorrTime(i, 0, evt2.firstcap[i][0]%1024)-(timewindow+windowsize)/2+0.5);
	  int n2=int(caltime.GetCorrTime(i, 0, evt2.firstcap[i][0]%1024)+(timewindow+windowsize)/2+0.5);
	  if (n2>=roi)
	    {
	      cout<<"n2 = "<<n2<<", >=ROI="<<roi<<endl;
	      n2=roi-1;
	    }
	  // cout<<n2-n1+1<<endl;
	  float time=0;
	  float sum = evt2.Sliding(i, 0, windowsize, &time, n1, n2) - meanped[i]*windowsize;
	  // time-=caltime.GetCorrTime(i, 0, evt2.firstcap[i][0]%1024);

	  float timecorr=caltime.GetCorrTime(i, 0, evt2.firstcap[i][0]%1024);
   	  // float sum = evt2.Sliding(i, mygain, windowsize, &time, nfirst, nlast) - meanped[i]*windowsize;
   	  if ((sum>minsig)&&(sum<maxsig))
	    {
	      htime[i]->Fill(time);
	      htimecorr[i]->Fill(time-timecorr+caltime.fan[i][0][0]/2);
	      numev[i]++;

	      for (int j=1; j<=hshapesw[i]->GetNbinsX(); j++)
		{
		  double posx=hshapesw[i]->GetBinCenter(j);
		  int slice=(int)(posx+timecorr+0.5);
		  hshapenosw[i]->Fill(posx, evt2.samples[i][mygain][slice]-meanped[i]);

		  int slicesw=(int)(posx+time+0.5);
		  hshapesw[i]->Fill(posx, evt2.samples[i][mygain][slicesw]-meanped[i]);
		}
	    }
   	}
    }
  plik2.close();

  for (int i=0; i<nn; i++)
    {
      hshapenosw[i]->Scale(1./numev[i]);
      hshapesw[i]->Scale(1./numev[i]);
    }

  // for (int i=0; i<nn; i++)
  //   {
  //     hsigped[i]->Scale(1.*numsig/numped); // just to equalize the number of pedestals;
  //     int ibinmean = hsigped[i]->FindBin(hsigped[i]->GetMean());
  //     float sumped=hsigped[i]->Integral(0, ibinmean);
  //     float sumsig=hsig[i]->Integral(0, ibinmean);
  //     float nphe = -log(sumped/sumsig);
  //     hsigped[i]->Scale(sumsig/sumped);
  //     cout<<"Nphe ~ "<<nphe<<endl;

  //   }
  // // prephist(hsigmean);
  // // for (int i=0; i<nn; i++)
  // //   prephist(hsig[i]);



  TCanvas *ctime = new TCanvas ("ctime", "", 640, 480);
  ctime->Divide(4,2,0.001, 0.001);
  for (int i=0; i<nn; i++)
    {
      ctime->cd(i+1);
      // fitHist(hsig[i], fstat[i]);
      htimecorr[i]->Draw();
      htime[i]->Draw("sames");
    }

  TCanvas *cshape = new TCanvas ("cshape", "", 640, 480);
  cshape->Divide(4,2,0.001, 0.001);
  for (int i=0; i<nn; i++)
    {
      cshape->cd(i+1);
      // fitHist(hsig[i], fstat[i]);
      hshapesw[i]->Draw();
      hshapenosw[i]->Draw("sames");
    }

}


// gain calibration using single phe for one file
// nameped - file to get the pedestal calibration of each slice
// namein - file with single phe pulses
// nameped2 - if not empty a global correction of pedestal shift is taken from this one
void gaincalib(string nameped, string namein, int roi, FFactorStat *fstat=0, int mygain=0, string nameped2="")
{
  const int nn=nch-1;
  const int windowsize=6;
  PedestalSimple &pedsim=*prepPedSimple(nameped, roi);

  // int nfirst=2, nlast=35;
  int nfirst=15-3, nlast=26+3;
  // PedStat *pedstat = getrms(nameped, roi, &pedsim);
  // pedstat->print();
  // for (int i=0; i<nn; i++)
  //   fstat[i].pedrms=pedstat->frms[i][mygain];

  TH1F *hsig[nn];
  TH1F *hsigped[nn];
  TH1F *htime[nn];
  for (int i=0; i<nn; i++)
    {
      htime[i]=new TH1F (Form("htime%i", i), ";arrival time[slices]", 2*roi, -0.5, roi-0.5);
      // hsig[i]=new TH1F (Form("hsig%i", i), ";signal[counts]", 100, -100, 500);
      // hsigped[i]=new TH1F (Form("hsigped_%i", i), ";signal[counts]", 100, -100, 500);
      hsig[i]=new TH1F (Form("hsig%i", i), ";signal[counts]", 100, -50, 250);
      hsigped[i]=new TH1F (Form("hsigped_%i", i), ";signal[counts]", 100, -50, 250);
      hsigped[i]->SetLineColor(kRed);
      hsigped[i]->Sumw2();
    }
  // first loop over the pedestal file
  ifstream plik1(nameped.data(), ios::binary);
  if (!plik1.is_open())
    {
      cout<<"file: "<<nameped<<" is not open, exiting"<<endl;
      return;
    }
  Event evt1(roi);
  int numped=0, numsig=0;
  for (int ev=0; ev<50000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt1.read(plik1))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      removePed(evt1, pedsim);
      evt1.CorrTime();
      evt1.InterpolatePseudoPulses();
      for (int i=0; i<nn; i++)
	hsigped[i]->Fill(evt1.Sliding(i, mygain, windowsize, 0, nfirst, nlast));
      numped++;
    }
  plik1.close();

  double meanped[nch]={0};
  if (nameped2!="")
    {
      cout<<" getting the pedestal correction from "<<nameped2<<endl;
      TH1F *hpeds[nch];
      for (int i=0;i<nch; i++)
	hpeds[i]=new TH1F (Form("hpeds_%i", i), "", 400, -100, 100);

      // loop over the second (newer) pedestal file
      ifstream plik1b(nameped2.data(), ios::binary);
      if (!plik1b.is_open())
	{
	  cout<<"file: "<<nameped2<<" is not open, exiting"<<endl;
	  return;
	}
      Event evt1b(roi);
      int numped=0, numsig=0;
      for (int ev=0; ev<10000; ev++)
	{
	  if (ev % 10000 == 0)
	    cout<<"Event "<<ev<<endl;
	  if (!evt1b.read(plik1b))
	    {
	      cout<<"event "<<ev<<"end of file"<<endl;
	      break;
	    }

	  removePed(evt1b, pedsim);
	  evt1b.CorrTime();
	  evt1b.InterpolatePseudoPulses();
	  for (int i=0; i<nn; i++)
	    for (int j=nfirst; j<=nlast; j++)
	      hpeds[i]->Fill(evt1b.samples[i][mygain][j]);
	}
      plik1b.close();
      for (int i=0; i<nn; i++)
	{
	  meanped[i]=hpeds[i]->GetMean();
	  cout<<"Ch "<<i<<", mean ped="<<meanped[i]<<", RMS="<<hpeds[i]->GetRMS()<<endl;
	}
    }


  ifstream plik2(namein.data(), ios::binary);
  if (!plik2.is_open())
    {
      cout<<"file: "<<namein<<" is not open, exiting"<<endl;
      return;
    }
  Event evt2(roi);

  for (int ev=0; ev<10000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt2.read(plik2))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      removePed(evt2, pedsim);
      evt2.CorrTime();
      evt2.InterpolatePseudoPulses();
      for (int i=0; i<nn; i++)
	{
	  float time=0;
	  float sum = evt2.Sliding(i, mygain, windowsize, &time, nfirst, nlast) - meanped[i]*windowsize;
	  hsig[i]->Fill(sum);
	  if (sum>100)
	    htime[i]->Fill(time);
	}
      numsig++;
    }
  plik2.close();

  for (int i=0; i<nn; i++)
    {
      hsigped[i]->Scale(1.*numsig/numped); // just to equalize the number of pedestals;
      int ibinmean = hsigped[i]->FindBin(hsigped[i]->GetMean());
      float sumped=hsigped[i]->Integral(0, ibinmean);
      float sumsig=hsig[i]->Integral(0, ibinmean);
      float nphe = -log(sumped/sumsig);
      hsigped[i]->Scale(sumsig/sumped);
      cout<<"Nphe ~ "<<nphe<<endl;

    }
  // prephist(hsigmean);
  // for (int i=0; i<nn; i++)
  //   prephist(hsig[i]);

  TCanvas *csig = new TCanvas ("csig", "", 1024, 768);
  csig->Divide(4,2,0.001, 0.001);
  for (int i=0; i<nn; i++)
    {
      csig->cd(i+1);
      // fitHist(hsig[i], fstat[i]);
      TF1 *f = new TF1("f", funcgain, -100, 500, 6);
      f->SetParNames("ped0", "pedrms", "nphe", "1pherms", "conv", "norm");
      f->SetParameters(30, 10, 0.1, 30, 150, 10000);
      // f->SetParameters(5, 10, 1, 20, 30, 1000);
      f->SetParLimits(1, 0, 100);
      f->SetParLimits(3, 0, 100);
      hsig[i]->Draw();
      hsigped[i]->Draw("same");
      hsig[i]->Fit(f, "R");
      TF1 *ffit = hsig[i]->GetFunction("f");
      TPaveText *pave = new TPaveText(0.5, 0.6, 0.99, 0.99, "NDC");
      pave->AddText(Form("<ped>=%.1fcnts", f->GetParameter(0)));
      pave->AddText(Form("RMS(ped)=%.1fcnts", f->GetParameter(1)));
      pave->AddText(Form("<gain>=%.1fcnts", f->GetParameter(4)));
      pave->AddText(Form("RMS(gain)=%.1fcnts", f->GetParameter(3)));
      pave->AddText(Form("<Nphe>=%.2f", f->GetParameter(2)));
      pave->AddText(Form("chi2/ndof=%.1f/%i",f->GetChisquare(), f->GetNDF()));
      pave->Draw();
      if (fstat)
	{
	  fstat[i].conv=1/f->GetParameter(4);
	  fstat[i].pedrms=f->GetParameter(1);
	}
    }


  TCanvas *ctime = new TCanvas ("ctime", "", 640, 480);
  ctime->Divide(4,2,0.001, 0.001);
  for (int i=0; i<nn; i++)
    {
      ctime->cd(i+1);
      // fitHist(hsig[i], fstat[i]);
      htime[i]->Draw();
    }

}

void gaincalall()
{
  const int nnames=19;
  string names[nnames]={
    "../data/Ped_Random1000_IC_ET_RD40_FEB0_IP211.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB1_IP217.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB2_IP227.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB3_IP230.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB4_IP233.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB5_IP244.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB6_IP265.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB7_IP268.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB8_IP269.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB9_IP274.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB10_IP280.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB11_IP284.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB12_IP285.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB13_IP292.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB14_IP308.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB15_IP248.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB16_IP317.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB17_IP333.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB18_IP337.dat"};

  string names2[nnames]={
    "../data/DATA_Filter16RD40_FEB0_IP211.dat",
    "../data/DATA_Filter16RD40_FEB1_IP217.dat",
    "../data/DATA_Filter16RD40_FEB2_IP227.dat",
    "../data/DATA_Filter16RD40_FEB3_IP230.dat",
    "../data/DATA_Filter16RD40_FEB4_IP233.dat",
    "../data/DATA_Filter16RD40_FEB5_IP244.dat",
    "../data/DATA_Filter16RD40_FEB6_IP265.dat",
    "../data/DATA_Filter16RD40_FEB7_IP268.dat",
    "../data/DATA_Filter16RD40_FEB8_IP269.dat",
    "../data/DATA_Filter16RD40_FEB9_IP274.dat",
    "../data/DATA_Filter16RD40_FEB10_IP280.dat",
    "../data/DATA_Filter16RD40_FEB11_IP284.dat",
    "../data/DATA_Filter16RD40_FEB12_IP285.dat",
    "../data/DATA_Filter16RD40_FEB13_IP292.dat",
    "../data/DATA_Filter16RD40_FEB14_IP308.dat",
    "../data/DATA_Filter16RD40_FEB15_IP248.dat",
    "../data/DATA_Filter16RD40_FEB16_IP317.dat",
    "../data/DATA_Filter16RD40_FEB17_IP333.dat",
    "../data/DATA_Filter16RD40_FEB18_IP337.dat"};

  FFactorStat fstat[nnames][nch-1];
  for (int i=0; i<nnames; i++)
    gaincalib(names[i], names2[i], 40, fstat[i], 0);

  TH1F *hconv = new TH1F ("hconv", ";<conversion factor, HG>;number of pixels ", 100, 0, 0.02);
  TH1F *hnoise = new TH1F ("hnoise", ";<noise RMS/[phe], HG>;number of pixels ", 100, 0, 0.6);

  for (int i=0; i<nnames; i++)
    for (int j=0; j<nch-1; j++)
      {
	hconv->Fill(fstat[i][j].conv);
	hnoise->Fill(fstat[i][j].pedrms * fstat[i][j].conv);
	cout<<fstat[i][j].conv<<" "<<fstat[i][j].pedrms<<endl;
      }

  TCanvas *c = new TCanvas ("c", "", 1200, 600);
  c->Divide(2,1,0.001,0.001);
  c->cd(1);
  hconv->Draw();
  c->cd(2);
  hnoise->Draw();
}

void ffactorall()
{
  const int nnames=19;
  string names[nnames]={
    "../data/Ped_Random1000_IC_ET_RD40_FEB0_IP211.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB1_IP217.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB2_IP227.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB3_IP230.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB4_IP233.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB5_IP244.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB6_IP265.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB7_IP268.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB8_IP269.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB9_IP274.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB10_IP280.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB11_IP284.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB12_IP285.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB13_IP292.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB14_IP308.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB15_IP248.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB16_IP317.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB17_IP333.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB18_IP337.dat"};

  string names2[nnames]={
    "../data/DATA_Filter6RD40_FEB0_IP211.dat",
    "../data/DATA_Filter6RD40_FEB1_IP217.dat",
    "../data/DATA_Filter6RD40_FEB2_IP227.dat",
    "../data/DATA_Filter6RD40_FEB3_IP230.dat",
    "../data/DATA_Filter6RD40_FEB4_IP233.dat",
    "../data/DATA_Filter6RD40_FEB5_IP244.dat",
    "../data/DATA_Filter6RD40_FEB6_IP265.dat",
    "../data/DATA_Filter6RD40_FEB7_IP268.dat",
    "../data/DATA_Filter6RD40_FEB8_IP269.dat",
    "../data/DATA_Filter6RD40_FEB9_IP274.dat",
    "../data/DATA_Filter6RD40_FEB10_IP280.dat",
    "../data/DATA_Filter6RD40_FEB11_IP284.dat",
    "../data/DATA_Filter6RD40_FEB12_IP285.dat",
    "../data/DATA_Filter6RD40_FEB13_IP292.dat",
    "../data/DATA_Filter6RD40_FEB14_IP308.dat",
    "../data/DATA_Filter6RD40_FEB15_IP248.dat",
    "../data/DATA_Filter6RD40_FEB16_IP317.dat",
    "../data/DATA_Filter6RD40_FEB17_IP333.dat",
    "../data/DATA_Filter6RD40_FEB18_IP337.dat"};

  FFactorStat fstat[nnames][nch-1];
  for (int i=0; i<nnames; i++)
    ffactor(names[i], names2[i], 40, fstat[i], 0, true);

  TH1F *hnphe = new TH1F ("hnphe", ";<Nphe>;number of pixels ", 100, 0, 150);
  TH1F *hconv = new TH1F ("hconv", ";<conversion factor, HG>;number of pixels ", 100, 0, 0.02);
  TH1F *hnoise = new TH1F ("hnoise", ";<noise RMS/[phe], HG>;number of pixels ", 100, 0, 0.6);
  TH1F *hpeak = new TH1F ("hpeak", ";1 phe equiv. peak;number of pixels ", 100, 0, 70);

  for (int i=0; i<nnames; i++)
    for (int j=0; j<nch-1; j++)
      {
	fstat[i][j].print();
	if (fstat[i][j].isgood)
	  {
	    hnphe->Fill(fstat[i][j].nphe);
	    hconv->Fill(fstat[i][j].conv);
	    hnoise->Fill(fstat[i][j].pedrms * fstat[i][j].conv);
	    hpeak->Fill(fstat[i][j].peak);
	  }
      }
  TCanvas *c = new TCanvas ("c", "", 1200, 600);
  c->Divide(3,1,0.001,0.001);
  c->cd(1);
  hnphe->Draw();
  c->cd(2);
  hconv->Draw();
  c->cd(3);
  hnoise->Draw();

  TCanvas *c2 = new TCanvas ("c2", "", 800, 600);
  hpeak->Draw();

  for (int i=0; i<nnames; i++)
    for (int j=0; j<nch-1; j++)
      cout<<fstat[i][j].conv<<" "<<fstat[i][j].pedrms<<endl;

}

void ffactorlin()
{
  float f2fact=1.2;
  float laserinstab=0.;
  int roi=40;
  const int nn=nch-1;
  const int nnames=19;
  string names[nnames]={
    "../data/Ped_Random1000_IC_ET_RD40_FEB0_IP211.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB1_IP217.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB2_IP227.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB3_IP230.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB4_IP233.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB5_IP244.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB6_IP265.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB7_IP268.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB8_IP269.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB9_IP274.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB10_IP280.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB11_IP284.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB12_IP285.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB13_IP292.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB14_IP308.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB15_IP248.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB16_IP317.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB17_IP333.dat",
    "../data/Ped_Random1000_IC_ET_RD40_FEB18_IP337.dat"};

  string names2[nnames]={
    "../data/DATA_Filter5RD40_FEB0_IP211.dat",
    "../data/DATA_Filter5RD40_FEB1_IP217.dat",
    "../data/DATA_Filter5RD40_FEB2_IP227.dat",
    "../data/DATA_Filter5RD40_FEB3_IP230.dat",
    "../data/DATA_Filter5RD40_FEB4_IP233.dat",
    "../data/DATA_Filter5RD40_FEB5_IP244.dat",
    "../data/DATA_Filter5RD40_FEB6_IP265.dat",
    "../data/DATA_Filter5RD40_FEB7_IP268.dat",
    "../data/DATA_Filter5RD40_FEB8_IP269.dat",
    "../data/DATA_Filter5RD40_FEB9_IP274.dat",
    "../data/DATA_Filter5RD40_FEB10_IP280.dat",
    "../data/DATA_Filter5RD40_FEB11_IP284.dat",
    "../data/DATA_Filter5RD40_FEB12_IP285.dat",
    "../data/DATA_Filter5RD40_FEB13_IP292.dat",
    "../data/DATA_Filter5RD40_FEB14_IP308.dat",
    "../data/DATA_Filter5RD40_FEB15_IP248.dat",
    "../data/DATA_Filter5RD40_FEB16_IP317.dat",
    "../data/DATA_Filter5RD40_FEB17_IP333.dat",
    "../data/DATA_Filter5RD40_FEB18_IP337.dat"};

  PedestalSimple *pedsims[nnames];
  for (int i=0; i<nnames; i++)
    pedsims[i]=prepPedSimple(names[i], roi, 5000);

  PedStat *pedstats[nnames];
  for (int i=0; i<nnames; i++)
    {
      pedstats[i] = getrms(names[i], roi, pedsims[i]);
      pedstats[i]->print();
    }

  const int nfilt=17;
  TH1F *hsigmeanhi[nfilt];
  TH1F *hsigmeanlo[nfilt];

  TH1F *hphehi[nfilt];
  TH1F *hphelo[nfilt];

  TGraph *grmslaserhi = new TGraph(); // RMS/Mean from the average pulse
  TGraph *grmslaserlo = new TGraph(); // RMS/Mean from the average pulse

  TGraph *gphehi = new TGraph(); // mean phe vs signal
  TGraph *gphelo = new TGraph(); // mean phe vs signal

  TGraph *gconvhi  = new TGraph(); // mean conv vs signal
  TGraph *gconvlo  = new TGraph(); // mean conv vs signal

  grmslaserhi->SetLineColor(kRed);
  grmslaserhi->SetMarkerColor(kRed);

  gphehi->SetLineColor(kRed);
  gphehi->SetMarkerColor(kRed);

  gconvhi->SetLineColor(kRed);
  gconvhi->SetMarkerColor(kRed);

  double meansighi[nfilt]={0}, meansiglo[nfilt]={0}, meanphehi[nfilt]={0}, meanphelo[nfilt]={0}, meanconvhi[nfilt]={0}, meanconvlo[nfilt]={0};

  for (int l = 0; l<nfilt; l++)
    {
      hsigmeanhi[l]=new TH1F (Form("hsigmeanhi%i",l), ";signal[counts]", 4096, 0, 100000);
      hsigmeanlo[l]=new TH1F (Form("hsigmeanlo%i",l), ";signal[counts]", 4096, 0, 100000);

      hphehi[l] = new TH1F (Form("hphehi%i",l), "", 4096, 0, 1000);
      hphelo[l] = new TH1F (Form("hphelo%i",l), "", 4096, 0, 1000);

      hsigmeanhi[l]->SetLineColor(kRed);
      hphehi[l]->SetLineColor(kRed);
      Event *evts[nnames];
      ifstream *files[nnames];
      for (int i=0; i<nnames; i++)
	{
	  evts[i] = new Event(roi);
	  TString str = names2[i];
	  str.ReplaceAll("Filter5", Form("Filter%i", l));
	  files[i] = new ifstream(str.Data(), ios::binary);
	}

      TH1F *hsig[nnames][nn][2];
      for (int i=0; i<nnames; i++)
	for (int j=0; j<nn; j++)
	  {
	    hsig[i][j][0]=new TH1F (Form("hsig_%i_%i_hi", i, j), ";signal[counts]", 4096, 0, 100000);
	    hsig[i][j][1]=new TH1F (Form("hsig_%i_%i_lo", i, j), ";signal[counts]", 4096, 0, 100000);
	  }

      for (int ev=0; ev<1000; ev++)
	{
	  if (ev % 10000 == 0)
	    cout<<"Event "<<ev<<endl;
	  if (!constructevent(evts, files, nnames))
	    {
	      cout<<"event "<<ev<<"end of file"<<endl;
	      break;
	    }

	  for (int i=0; i<nnames; i++)
	    {
	      removePed(*evts[i], *pedsims[i]);
	      evts[i]->CorrTime();
	      evts[i]->InterpolatePseudoPulses();
	    }

	  float sumavhi=0, sumavlo=0;
	  for (int i=0; i<nnames; i++)
	    for (int j=0; j<nn; j++)
	      for (int k=0; k<2; k++)
		{
		  float sum = evts[i]->Sliding(j, k);
		  hsig[i][j][k]->Fill(sum);
		  if (k==0)
		    sumavhi+=sum;
		  else
		    sumavlo+=sum;
		}
	  sumavhi/=(nn*nnames);
	  sumavlo/=(nn*nnames);
	  hsigmeanhi[l]->Fill(sumavhi);
	  hsigmeanlo[l]->Fill(sumavlo);
	}

      grmslaserhi->SetPoint(l, hsigmeanhi[l]->GetMean(), hsigmeanhi[l]->GetRMS()/hsigmeanhi[l]->GetMean());
      grmslaserlo->SetPoint(l, hsigmeanlo[l]->GetMean(), hsigmeanlo[l]->GetRMS()/hsigmeanlo[l]->GetMean());

      // now the Ffactor method
      int nall = nnames * nn;
      double meansig[nall][2], rmssig[nall][2], nphe[nall][2], conv[nall][2], pedrms[nall][2];

      TH1F *hconvhi = new TH1F ("hconvhi", "", 1000, 1.e-5, 1);
      TH1F *hconvlo = new TH1F ("hconvhi", "", 1000, 1.e-5, 1);
      for (int i=0; i<nnames; i++)
	for (int j=0; j<nn; j++)
	  for (int k=0; k<2; k++)
	    {
	      int ipix=i*nn+j;
	      meansig[ipix][k]= hsig[i][j][k]->GetMean();
	      rmssig[ipix][k]= hsig[i][j][k]->GetRMS();
	      pedrms[ipix][k]= pedstats[i]->frms[j][k];
	      double rmslaser = laserinstab*meansig[ipix][k];
	      nphe[ipix][k] = f2fact*meansig[ipix][k] * meansig[ipix][k] /(rmssig[ipix][k]*rmssig[ipix][k] - pedrms[ipix][k]*pedrms[ipix][k] - rmslaser*rmslaser);
	      conv[ipix][k]=nphe[ipix][k]/meansig[ipix][k];
	      cout<<"pix "<<ipix<<(k==0?"hi":"lo")<<", <signal>="<<meansig[ipix][k]<<", nphe="<<nphe[ipix][k]<<", conv="<<conv[ipix][k]<<endl;
	      if (k==0)
		{
		  hphehi[l]->Fill(nphe[ipix][k]);
		  hconvhi->Fill(conv[ipix][k]);
		  meansighi[l]+=meansig[ipix][k];
		}
	      else
		{
		  hphelo[l]->Fill(nphe[ipix][k]);
		  hconvlo->Fill(conv[ipix][k]);
		  meansiglo[l]+=meansig[ipix][k];
		}
	    }
      meansighi[l]/=(nn*nnames);
      meansiglo[l]/=(nn*nnames);
      meanphehi[l] = hphehi[l]->GetMean();
      meanphelo[l] = hphelo[l]->GetMean();
      meanconvhi[l] = hconvhi->GetMean();
      meanconvlo[l] = hconvlo->GetMean();

      gphehi->SetPoint(l, meansighi[l], meanphehi[l]);
      gphelo->SetPoint(l, meansighi[l], meanphelo[l]); // also as function of meansighi !

      gconvhi->SetPoint(l, meansighi[l], meanconvhi[l]);
      gconvlo->SetPoint(l, meansiglo[l], meanconvlo[l]);

      delete hconvhi;
      delete hconvlo;
      for (int i=0; i<nnames; i++)
	for (int j=0; j<nn; j++)
	  for (int k=0; k<2; k++)
	    delete hsig[i][j][k];
    }
  for (int l = 0; l<nfilt; l++)
    cout<<"Filter "<<l<<", HI, <signal>"<<meansighi[l]<<", <phe>="<<meanphehi[l]<<", <conv>="<<meanconvhi[l]<<
      " LO, <signal>"<<meansiglo[l]<<", <phe>="<<meanphelo[l]<<", <conv>="<<meanconvlo[l]<<endl;

  TCanvas *claserhi = new TCanvas ("claserhi", "", 1024, 768);
  claserhi->Divide(5, 4, 0.001, 0.001);
  for (int l=0; l<nfilt; l++)
    {
      claserhi->cd(l+1);
      prephist(hsigmeanhi[l]);
      hsigmeanhi[l]->Draw();
    }

  TCanvas *claserlo = new TCanvas ("claserlo", "", 1024, 768);
  claserlo->Divide(5, 4, 0.001, 0.001);
  for (int l=0; l<nfilt; l++)
    {
      claserlo->cd(l+1);
      prephist(hsigmeanlo[l]);
      hsigmeanlo[l]->Draw();
    }

  TCanvas *cphe = new TCanvas ("cphe", "", 1024, 768);
  cphe->Divide(5, 4, 0.001, 0.001);
  for (int l=0; l<nfilt; l++)
    {
      cphe->cd(l+1);
      prephist(hphehi[l]);
      hphehi[l]->Draw();
      prephist(hphelo[l]);
      hphelo[l]->Draw("sames");
    }


  TCanvas *cc = new TCanvas("cc", "", 1024, 768);
  cc->Divide(2, 2, 0.001, 0.001);
  for (int i=1; i<=4; i++)
    cc->GetPad(i)->SetLogx();
  cc->cd(1);
  grmslaserhi->Draw("A*L");
  grmslaserhi->GetXaxis()->SetTitle("average signal [counts]");
  grmslaserhi->GetYaxis()->SetTitle("average RMS/average signal");
  grmslaserlo->Draw("*L");
  cc->cd(2);
  cc->GetPad(2)->SetLogy();
  gphehi->Draw("A*L");
  gphehi->GetXaxis()->SetTitle("average signal in hi gain[counts]");
  gphehi->GetYaxis()->SetTitle("average Nphe");
  gphelo->Draw("*L");
  cc->cd(3);
  gconvhi->Draw("A*L");
  gconvhi->GetXaxis()->SetTitle("average signal [counts]");
  gconvhi->GetYaxis()->SetTitle("average conversion factor [phe/counts]");
  cc->cd(4);
  gconvlo->Draw("A*L");
  gconvlo->GetXaxis()->SetTitle("average signal [counts]");
  gconvlo->GetYaxis()->SetTitle("average conversion factor [phe/counts]");

}


//==================================================================
// checking charge resolution vs Nphe for one pixel
// nameped - pedestal file
// namecal - calibration file to get coeff
// namelin - linearity files
void CheckChargeResVsPheOne(int roi=40)
{
  const int nn=nch-1;
  const float NSBrate=-0.2;
  int windowsize=10;
  const Double_t timewindow=-1;
  // const Double_t timewindow=3;
  string nameped="../TimeResStudy/PedestalData/PedTable_20170901_160812_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB0_IP237.dat";
  string namecallo="../TimeResStudy/CalibrationData/DATA_Linearity_20170901_174631/Linearity_NominalHV_Filter5RD40_FEB0_IP237.dat";
  string namecalhi="../TimeResStudy/CalibrationData/DATA_Linearity_20170901_174631/Linearity_NominalHV_Filter7RD40_FEB0_IP237.dat";
  string namelin0="../TimeResStudy/LinearityData/DATA_Linearity_20170901_173134/Linearity_NominalHV_Filter0RD40_FEB0_IP237.dat";

  // string nameped="../TimeResStudy/PedestalData/PedTable_20170901_160812_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB2_IP76.dat";
  // string namecallo="../TimeResStudy/CalibrationData/DATA_Linearity_20170901_174631/Linearity_NominalHV_Filter5RD40_FEB2_IP76.dat";
  // string namecalhi="../TimeResStudy/CalibrationData/DATA_Linearity_20170901_174631/Linearity_NominalHV_Filter7RD40_FEB2_IP76.dat";
  // string namelin0="../TimeResStudy/LinearityData/DATA_Linearity_20170901_173134/Linearity_NominalHV_Filter0RD40_FEB2_IP76.dat";

  TF1 *resreq = new TF1("resreq", "sqrt(x*pow(1+[0],2)+[1]+pow([2]*x,2))/x", 1, 1000);
  resreq->SetParameters(0.2,sqrt(7), 0.1);
  resreq->SetLineWidth(3);
  resreq->SetLineColor(kRed);
  TF1 *resgoal = new TF1("resgoal", "sqrt(x*pow(1+[0],2)+[1]+pow([2]*x,2))/x", 1, 3000);
  resgoal->SetParameters(0.1152,2, 0.05);
  resgoal->SetLineWidth(3);
  resgoal->SetLineColor(kGreen);

  FFactorStat ffhi[nch];
  ffactor(nameped, namecalhi, roi, ffhi, 0, true, windowsize);
  FFactorStat fflo[nch];
  ffactor(nameped, namecallo, roi, fflo, 1, true, windowsize);

  for (int i=0; i<nn; i++)
    cout<<"Hi/Lo gain for pix"<<i<<" = "<<fflo[i].conv/ffhi[i].conv<<endl;

  PedestalSimple &pedsim=*prepPedSimple(nameped, roi);

  ifstream plik2(namecallo.data(), ios::binary);
  if (!plik2.is_open())
    {
      cout<<"file: "<<namecallo<<" is not open, exiting"<<endl;
      return;
    }
  Event evt2(roi);
  TimeCalClu caltime(1024, 8, 7);
  for (int ev=0; ev<19000; ev++)
    {
      if (ev % 10000 == 0)
   	cout<<"Event "<<ev<<endl;
      if (!evt2.read(plik2))
   	{
   	  cout<<"event "<<ev<<"end of file"<<endl;
   	  break;
   	}

      removePed(evt2, pedsim);
      evt2.CorrTime();
      evt2.InterpolatePseudoPulses();
      float times[nn][2];
      for (int i=0; i<nn; i++)
   	for (int j=0; j<2; j++)//gain
   	  {
   	    evt2.Sliding(i, j, windowsize, &times[i][j]);
   	    caltime.Fill(i, j, evt2.firstcap[i][j]%1024, times[i][j]);
   	  }
    }
  plik2.close();

  caltime.Finalize();
  caltime.Print();

  const Int_t filmin=0;
  const Int_t filmax=18;

  const Int_t nfilt=filmax-filmin+1; // number of filters
  TH1F *hcharge[nn][nfilt];
  TH1F *hchargelo[nn][nfilt];
  for (int i=0; i<nn; i++)
    for (int j=0; j<nfilt; j++)
      {
	hcharge[i][j]=new TH1F (Form("hchargehi%i_%i", i, j), "", 4096, -3, 1000);
	hchargelo[i][j]=new TH1F (Form("hchargelo%i_%i", i, j), "", 4096, -3, 2000);
      }

  // second loop over linearity files
  for (int filter=filmin; filter<=filmax; filter++)
    {
      TString namelin=namelin0;
      namelin.ReplaceAll("Filter0", Form("Filter%i", filter));
      ifstream plik3(namelin.Data(), ios::binary);
      if (!plik3.is_open())
	{
	  cout<<"file: "<<namelin<<" is not open, exiting"<<endl;
	  return;
	}
      Event evt3(roi);

      for (int ev=0; ev<2000; ev++)
	{
	  if (ev % 10000 == 0)
	    cout<<"Event "<<ev<<endl;
	  if (!evt3.read(plik3))
	    {
	      cout<<"event "<<ev<<"end of file"<<endl;
	      break;
	    }

	  removePed(evt3, pedsim);
	  evt3.CorrTime();
	  evt3.InterpolatePseudoPulses();

	  if (NSBrate>0)
	    {
	      for (int i=0; i<nn; i++)
		for (int j=0; j<2; j++)
		  {
		    Double_t onephe=0;
		    if (j==0) onephe = 1/ffhi[i].conv;
		    else onephe = 1/fflo[i].conv;
		    // AddNSB(Event &ev, float rate,  int pixid, int gain, float pulsewidth, float onephemean, float onepherms)
		    AddNSB(evt3, NSBrate, i, j, 1.2, onephe, onephe*0.35);
		  }
	    }

	  // high gain
	  for (int i=0; i<nn; i++)
	    {
	      float time;
	      int n1=int(caltime.GetCorrTime(i, 0, evt3.firstcap[i][0]%1024)-(timewindow+windowsize)/2+0.5);
	      int n2=int(caltime.GetCorrTime(i, 0, evt3.firstcap[i][0]%1024)+(timewindow+windowsize)/2+0.5);
	      if (n2>=roi)
		{
		  cout<<"n2 = "<<n2<<", >=ROI="<<roi<<endl;
		  n2=roi-1;
		}
	      cout<<n2-n1+1<<endl;
	      float sum = evt3.Sliding(i, 0, windowsize, &time, n1, n2);
	      float nphe = sum*ffhi[i].conv;
	      time-=caltime.GetCorrTime(i, 0, evt3.firstcap[i][0]%1024);

	      hcharge[i][filter-filmin]->Fill(nphe);
	    }

	  // lo gain
	  for (int i=0; i<nn; i++)
	    {
	      float time;
	      int n1=int(caltime.GetCorrTime(i, 1, evt3.firstcap[i][1]%1024)-(timewindow+windowsize)/2+0.5);
	      int n2=int(caltime.GetCorrTime(i, 1, evt3.firstcap[i][1]%1024)+(timewindow+windowsize)/2+0.5);
	      cout<<n2-n1+1<<endl;
	      if (n2>=roi)
		{
		  cout<<"LO gain n2 = "<<n2<<", >=ROI="<<roi<<endl;
		  n2=roi-1;
		}
	      float sum = evt3.Sliding(i, 1, windowsize, &time, n1, n2);
	      float nphe = sum*fflo[i].conv;
	      time-=caltime.GetCorrTime(i, 1, evt3.firstcap[i][0]%1024);

	      hchargelo[i][filter-filmin]->Fill(nphe);
	    }
	}
      plik3.close();
    }

  TGraph *gcharge[nn];
  TGraph *gcharge2[nn];
  for (int i=0; i<nn; i++)
    {
      gcharge[i] = new TGraph();
      gcharge2[i] = new TGraph();
    }

  for (int i=0; i<nn; i++)
    for (int ifil=0; ifil<nfilt; ifil++)
      {
	Double_t mean=hcharge[i][ifil]->GetMean();
	Double_t rms=hcharge[i][ifil]->GetRMS();
	if (mean>200) //switching to low gain
	  {
	    mean=hchargelo[i][ifil]->GetMean();
	    rms=hchargelo[i][ifil]->GetRMS();
	  }

	gcharge[i]->SetPoint(gcharge[i]->GetN(), mean, rms/mean);
	gcharge2[i]->SetPoint(gcharge2[i]->GetN(), mean, rms);

	TH1F *h=hcharge[i][ifil];
	Int_t nentr=h->GetEntries();
	if (nentr>10)
	  while(h->GetMaximum()<TMath::Sqrt(nentr))
	    h->Rebin(2);
	h->GetXaxis()->SetRangeUser(mean-3*rms, mean+3*rms);
      }


  TF1 *fpoisson=new TF1 ("fpoisson", "1/sqrt(x)", 1, 2000);
  fpoisson->SetLineWidth(3);
  fpoisson->SetLineColor(kGray);

  TLegend *leg = new TLegend (0.5, 0.7, 0.99, 0.99);
  TString header=Form("%s, %s",
		      (NSBrate>0)?Form("NSB %i MHz", (int)(NSBrate*1000)):"no NSB",
		      (timewindow>0)?Form("sl. window %i of %i", windowsize, windowsize+(int)timewindow+1):Form("fixed window of %i", windowsize));
  leg->SetHeader(header);
  leg->AddEntry(resreq, "Requirement", "L");
  leg->AddEntry(resgoal, "Goal", "L");
  leg->AddEntry(fpoisson, "Poisson limit", "L");
  leg->AddEntry(gcharge[0], "individual channels", "L*");

  TCanvas *c2 = new TCanvas ("c2", "", 800, 600);
  TH1F *hosie2 = c2->DrawFrame(1, 0.01, 2000, 2.1);
  c2->SetLogx();
  c2->SetLogy();
  for (int i=0; i<nn; i++)
    gcharge[i]->Draw("*L");
  fpoisson->Draw("same");
  hosie2->GetXaxis()->SetTitle("Nphe");
  hosie2->GetYaxis()->SetTitle("#Delta Q / Q");
  resreq->Draw("same");
  resgoal->Draw("same");
  leg->Draw();

  TF1 *fpoisson2=new TF1 ("fpoisson2", "sqrt(x)", 1, 2000);
  fpoisson2->SetLineWidth(3);
  fpoisson2->SetLineColor(kGray);

  TCanvas *c2b = new TCanvas ("c2b", "", 800, 600);
  TH1F *hosie2b = c2b->DrawFrame(1, 0.5, 2000, 60);
  c2b->SetLogx();
  c2b->SetLogy();
  for (int i=0; i<nn; i++)
    gcharge2[i]->Draw("*L");
  fpoisson2->Draw("same");
  hosie2b->GetXaxis()->SetTitle("Nphe");
  hosie2b->GetYaxis()->SetTitle("#Delta Nphe");

  // TF1 *ftimefit = new TF1("ftimefit", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+[2]*[2])", 1, 300);
  // ftimefit->SetParameters(1, 1, 1);
  // gtimerms->Fit(ftimefit);

  TCanvas *c3 = new TCanvas ("c3", "", 800, 600);
  c3->Divide(5, 4, 0.001, 0.001);
  for (int i=0; i<nfilt; i++)
    {
      c3->cd(i+1);
      hcharge[0][i]->Draw();
    }
  TCanvas *c3lo = new TCanvas ("c3lo", "", 800, 600);
  c3lo->Divide(5, 4, 0.001, 0.001);
  for (int i=0; i<nfilt; i++)
    {
      c3lo->cd(i+1);
      hchargelo[0][i]->Draw();
    }

}


void readdragon()
{
  gStyle->SetOptStat(111111);
  int nwindow=6;
  TH1F *hhi = new TH1F ("hhi", Form(";noise RMS in %i cap./sqrt(%i) ; number of pixels", nwindow, nwindow), 100, 0, 20);
  TH1F *hlo = new TH1F ("hlo", "", 100, 0, 20);
  const int nnames=19;
  // string names[nnames]={
  //   "Ped_Random1000_IC_ET_RD40_FEB0_IP211.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB1_IP217.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB2_IP227.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB3_IP230.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB4_IP233.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB5_IP244.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB6_IP265.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB7_IP268.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB8_IP269.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB9_IP274.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB10_IP280.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB11_IP284.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB12_IP285.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB13_IP292.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB14_IP308.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB15_IP248.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB16_IP317.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB17_IP333.dat",
  //   "Ped_Random1000_IC_ET_RD40_FEB18_IP337.dat"};

// string names[nnames]= {
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB0_IP147.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB1_IP323.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB2_IP123.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB3_IP134.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB4_IP255.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB5_IP34.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB6_IP60.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB7_IP131.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB8_IP303.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB9_IP245.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB10_IP63.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB11_IP322.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB12_IP143.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB13_IP72.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB14_IP283.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB15_IP297.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB16_IP250.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB17_IP95.dat",
// "../data2/PedTable_20170508_121832_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB18_IP253.dat"  };

//   string names[nnames]= {
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB0_IP147.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB1_IP323.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB2_IP123.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB3_IP134.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB4_IP255.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB5_IP34.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB6_IP60.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB7_IP131.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB8_IP303.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB9_IP245.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB10_IP63.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB11_IP322.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB12_IP143.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB13_IP72.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB14_IP283.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB15_IP297.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB16_IP250.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB17_IP95.dat",
// "../data2/PedTable_20170508_122227_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB18_IP253.dat"
//   };

 string names[nnames]= {
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB0_IP62.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB1_IP84.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB2_IP76.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB3_IP315.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB4_IP25.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB5_IP216.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB6_IP116.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB7_IP120.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB8_IP263.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB9_IP65.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB10_IP205.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB11_IP270.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB12_IP87.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB13_IP331.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB14_IP7.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB15_IP258.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB16_IP237.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB17_IP110.dat",
"../data3/PedTable_20170512_115831_IC_ET_RD40_NomHV/Ped_Random1000_IC_ETRD40_FEB18_IP257.dat"};

//  string names[nnames]= {
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB0_IP62.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB1_IP84.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB2_IP76.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB3_IP315.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB4_IP25.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB5_IP216.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB6_IP116.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB7_IP120.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB8_IP263.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB9_IP65.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB10_IP205.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB11_IP270.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB12_IP87.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB13_IP331.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB14_IP7.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB15_IP258.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB16_IP237.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB17_IP110.dat",
// "../data3/PedTable_20170512_120758_IC_ET_RD40/Ped_Random1000_IC_ETRD40_FEB18_IP257.dat"};


  PedStat pedstat[nnames];
  TH1D *hped = new TH1D ("hped", Form(";Signal in %i cap/sqrt(%i) [counts];Number of entries",nwindow, nwindow), 200, -50, 150);
  for (int i=0; i<nnames; i++)
  // int i=0;
    pedstat[i] = makepedplot(names[i], nwindow, 40, hped);

  for (int i=0; i<nnames; i++)
    {
      cout<<names[i]<<", RMS = "<<endl;
      for (int j=0; j<nch-1; j++)
	{
	  cout<<pedstat[i].frms[j][0]<<", ";
	  hhi->Fill(pedstat[i].frms[j][0]);
	  hlo->Fill(pedstat[i].frms[j][1]);
	}
      cout<<endl;
      for (int j=0; j<nch-1; j++)
	cout<<pedstat[i].frms[j][1]<<", ";
      cout<<endl;
    }
  hhi->SetLineColor(kRed);
  TCanvas *canped = new TCanvas ("canped", "", 640, 480);
  hhi->Draw();
  hlo->Draw("sames");

  TCanvas *canped2 = new TCanvas ("canped2", "", 640, 480);
  canped2->SetLogy();
  hped->Draw();

  TFile plikout("plikout.root", "RECREATE");
  plikout.cd();
    {
      canped->Write("canped");
      canped2->Write("canped2");
    }
  plikout.Close();

}
//==================================================================
// checking the calibration of high to low gain for one
void hilowone(string nameped, string namein, int roi, float *sigshi=0, float *sigslo=0, float *ratios=0, bool batch=false)
{
  const int nn=nch-1;
  const int windowsize=6;
  PedestalSimple &pedsim=*prepPedSimple(nameped, roi);

  int nfirst=15-3, nlast=26+3;

  TH1F *hsighi[nn];
  TH1F *hsiglo[nn];
  TH1F *hratio[nn];
  for (int i=0; i<nn; i++)
    {
      hsighi[i]=new TH1F (Form("hsighi%i", i), ";signal[counts]", 4096, 0, 100000);
      hsiglo[i]=new TH1F (Form("hsiglo%i", i), ";signal[counts]", 4096, 0, 100000);
      hratio[i]=new TH1F (Form("hratio%i", i), ";Hi/Lo ratio", 100, 0, 50.);
    }

  ifstream plik2(namein.data(), ios::binary);
  if (!plik2.is_open())
    {
      cout<<"file: "<<namein<<" is not open, exiting"<<endl;
      return;
    }
  Event evt2(roi);

  for (int ev=0; ev<1000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt2.read(plik2))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      removePed(evt2, pedsim);
      evt2.CorrTime();
      evt2.InterpolatePseudoPulses();
      for (int i=0; i<nn; i++)
	{
	  float time=0;
	  float sumhi = evt2.Sliding(i, 0, windowsize, &time, nfirst, nlast);
	  float sumlo = evt2.Sliding(i, 1, windowsize, &time, nfirst, nlast);
	  hsighi[i]->Fill(sumhi);
	  hsiglo[i]->Fill(sumlo);
	  hratio[i]->Fill(sumhi/sumlo);
	}
    }
  plik2.close();

  for (int i=0; i<nn; i++)
    {
      prephist(hsighi[i]);
      prephist(hsiglo[i]);
      cout<<i<<" high="<<hsighi[i]->GetMean()<<", lo="<<hsiglo[i]->GetMean()<<", ratio = "<<hratio[i]->GetMean()<<endl;
      if (ratios)
	{
	  sigshi[i]=hsighi[i]->GetMean();
	  sigslo[i]=hsiglo[i]->GetMean();
	  ratios[i]=hratio[i]->GetMean();
	}
    }

  if (batch)
    {
      for (int i=0; i<nn; i++)
	{
	  delete hsighi[i];
	  delete hsiglo[i];
	  delete hratio[i];
	}
    }
  else
    {
      TCanvas *csighi = new TCanvas ("csighi", "", 1024, 768);
      csighi->Divide(4,2,0.001, 0.001);
      for (int i=0; i<nn; i++)
	{
	  csighi->cd(i+1);
	  hsighi[i]->Draw();
	}

      TCanvas *csiglo = new TCanvas ("csiglo", "", 1024, 768);
      csiglo->Divide(4,2,0.001, 0.001);
      for (int i=0; i<nn; i++)
	{
	  csiglo->cd(i+1);
	  hsiglo[i]->Draw();
	}

      TCanvas *cratio = new TCanvas ("cratio", "", 1024, 768);
      cratio->Divide(4,2,0.001, 0.001);
      for (int i=0; i<nn; i++)
	{
	  cratio->cd(i+1);
	  hratio[i]->Draw();
	}
    }

}

//==================================================================
// checking linearity of the high/low gain calibration for one module
void hilowonelin()
{
  int filtmin=1, filtmax=8;
  const int nn=nch-1;
  TGraph *ghi[nn];
  TGraph *glo[nn];
  for (int i=0; i<nn; i++)
    {
      ghi[i]=new TGraph();
      glo[i]=new TGraph();
      ghi[i]->SetLineColor(i+1);
      glo[i]->SetLineColor(i+1);
      ghi[i]->SetLineWidth(2);
      glo[i]->SetLineWidth(2);
      ghi[i]->SetMarkerStyle(20);
      glo[i]->SetMarkerStyle(20);
    }

  for (int fi=filtmin; fi<=filtmax; fi++)
    {
      float sigshi[nn], sigslo[nn], ratios[nn];
      hilowone("../data/Ped_Random1000_IC_ET_RD40_FEB0_IP211.dat",
	       Form("../data/DATA_Filter%iRD40_FEB0_IP211.dat", fi),
	       40, sigshi, sigslo, ratios, true);
      for (int i=0; i<nn; i++)
	{
	  ghi[i]->SetPoint(ghi[i]->GetN(), sigshi[i], ratios[i]);
	  glo[i]->SetPoint(glo[i]->GetN(), sigslo[i], ratios[i]);
	}
    }
  TCanvas *chi = new TCanvas ("chi", "", 640, 480);
  chi->SetLogx();
  TH1 *hosie1 = chi->DrawFrame(50, 5, 100000, 30);
  hosie1->GetXaxis()->SetTitle("Hi gain signal [counts]");
  hosie1->GetYaxis()->SetTitle("Hi to lo gain ratio");
  for (int i=0; i<nn; i++)
    ghi[i]->Draw("PL");

  TCanvas *clo = new TCanvas ("clo", "", 640, 480);
  clo->SetLogx();
  TH1 *hosie2 = clo->DrawFrame(50, 5, 100000, 30);
  hosie2->GetXaxis()->SetTitle("Lo gain signal [counts]");
  hosie2->GetYaxis()->SetTitle("Hi to lo gain ratio");
  for (int i=0; i<nn; i++)
    glo[i]->Draw("PL");

}

//==================================================================
// averaging points of similar position
TGraph *average (TGraph *g, int numsl=16, int minnum=5)
{
  int nbin = size4drs/numsl + 1;
  float *x = new float[nbin];
  float *y = new float[nbin];
  int *n = new int[nbin];
  for (int i=0; i<nbin; i++)
    {
      x[i]=(i+0.5)*numsl;
      y[i]=0; n[i]=0;
    }
  for (int i=0; i<g->GetN(); i++)
    {
      double xx, yy;
      g->GetPoint(i, xx, yy);
      int pos = ((int) xx) / numsl;
      y[pos]+=yy;
      n[pos]++;
    }
  TGraph *gav = new TGraph();
  for (int i=0; i<nbin; i++)
    if (n[i]>=minnum)
      gav->SetPoint(gav->GetN(), x[i], y[i]/n[i]);
  gav->SetMarkerColor(kRed);
  gav->SetLineColor(kRed);
  gav->SetLineWidth(2);
  return gav;

}

//==================================================================
// checking arrival times for one module
// nameped - pedestal file
// namecal - calibration file to get coeff
// nametest - another calibration file to test time
void CheckArrTimeOne(string nameped, string namecal, string nametest, int roi=40)
{
  const int nn=nch-1;
  int windowsize=8;
  // int windowsize=12;
  TGraph *gtime[nn];
  TH1F *htime0[nn];
  TH1F *htime1[nn];
  TH1F *htime2[nn];
  for (int i=0; i<nn; i++)
    {
      gtime[i]=new TGraph();
      gtime[i]->SetMarkerStyle(7);
      htime0[i]=new TH1F (Form("htime0%i", i), "", 100, -6,6);
      htime2[i]=new TH1F (Form("htime2%i", i), "", 100, -6,6);
      htime2[i]->SetLineColor(kRed);
    }

  PedestalSimple &pedsim=*prepPedSimple(nameped, roi);
  TimeCalClu caltime2(1024, 8, 7);
  caltime2.FillFile(namecal, pedsim, roi, windowsize);
  caltime2.Finalize();
  caltime2.Print();

  // second loop over an independent file
  ifstream plik3(nametest.data(), ios::binary);
  if (!plik3.is_open())
    {
      cout<<"file: "<<nametest<<" is not open, exiting"<<endl;
      return;
    }
  Event evt3(roi);

  for (int ev=0; ev<19000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt3.read(plik3))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      removePed(evt3, pedsim);
      evt3.CorrTime();
      evt3.InterpolatePseudoPulses();
      for (int i=0; i<nn; i++)
	{
	  float time;
	  float sumhi = evt3.Sliding(i, 0, windowsize, &time, -1, -1);

	  htime0[i]->Fill(time-caltime2.GetMeanTime(i, 0));
	  htime2[i]->Fill(time-caltime2.GetCorrTime(i, 0, evt3.firstcap[i][0]%1024));
	}
    }
  plik3.close();

  TCanvas *c2 = new TCanvas ("c2", "", 800, 600);
  c2->Divide(4,2, 0.001, 0.001);
  for (int i=0; i<nn; i++)
    {
      c2->cd(i+1);
      htime2[i]->Draw();
      htime0[i]->Draw("sames");
    }
}

//==================================================================
// checking pulse shapes for one module
// nameped - pedestal file
// namecal - calibration file to get coeff
// nametest - another calibration file to test
void CheckPulseShape(string nameped, string namecal, string nametest, int roi=40)
{
  const int nn=nch-1;
  int windowsize=8;

  PedestalSimple &pedsim=*prepPedSimple(nameped, roi);

  // TimeCalClu caltime2(1024, 8, 7);
  // caltime2.FillFile(namecal, pedsim, roi, windowsize);
  // caltime2.Finalize();

  TH2F *hshape[nn];
  for (int i=0; i<nn; i++)
    hshape[i]=new TH2F (Form("htime0%i", i), "", 12*4, -6,6, 100, -300, 3000);


  // second loop over an independent file
  ifstream plik3(nametest.data(), ios::binary);
  if (!plik3.is_open())
    {
      cout<<"file: "<<nametest<<" is not open, exiting"<<endl;
      return;
    }
  Event evt3(roi);

  for (int ev=0; ev<19000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt3.read(plik3))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      removePed(evt3, pedsim);
      evt3.CorrTime();
      evt3.InterpolatePseudoPulses();
      int skip=2;
      for (int ich=0; ich<nn; ich++)
	{
	  // float timeoffset=caltime2.GetCorrTime(ich, 0, evt3.firstcap[ich][0]%1024);
	  // for (int i=skip; i<roi-skip; i++)
	  //   hshape[ich]->Fill(i-timeoffset, evt3.samples[ich][0][i]);

	  float time;
	  float sumhi = evt3.Sliding(ich, 0, windowsize, &time, -1, -1);
	  for (int i=skip; i<roi-skip; i++)
	    hshape[ich]->Fill(i-time, evt3.samples[ich][0][i]);

	}
    }
  plik3.close();

  TCanvas *c2 = new TCanvas ("c2", "", 800, 600);
  c2->Divide(4,2, 0.001, 0.001);
  for (int i=0; i<nn; i++)
    {
      c2->cd(i+1);
      hshape[i]->Draw("colz");
      hshape[i]->ProfileX()->Draw("same");
      // gaussian with RMS ~1.7 slice
    }
}


//==================================================================
// checking time resolution vs Nphe for one pixel
// nameped - pedestal file
// namecal - calibration file to get coeff
// namelin - linearity files
void CheckTimeResVsPheOne(int roi=40)
{
  const int nn=nch-1;
  int windowsize=6;
  string nameped="../TimeResStudy/PedestalData/PedTable_20170901_160812_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB0_IP237.dat";
  string namecal="../TimeResStudy/CalibrationData/DATA_Linearity_20170901_174631/Linearity_NominalHV_Filter5RD40_FEB0_IP237.dat";
  string namelin0="../TimeResStudy/LinearityData/DATA_Linearity_20170901_173134/Linearity_NominalHV_Filter0RD40_FEB0_IP237.dat";
  FFactorStat ff[nch];
  ffactor(nameped, namecal, roi, ff, 0, true, windowsize);

  TH2F *htimevsnphe = new TH2F("htimevsnphe", "", 100, 0, 3, 100, -5, 5);

  // TH1F *htime0[nn];
  // TH1F *htime1[nn];
  // TH1F *htime2[nn];
  // for (int i=0; i<nn; i++)
  //   {
  //     gtime[i]=new TGraph();
  //     gtime[i]->SetMarkerStyle(7);
  //     // htime[i]=new TH1F (Form("htime%i", i), "", 160, 0, 40);
  //     htime0[i]=new TH1F (Form("htime0%i", i), "", 100, -6,6);
  //     htime1[i]=new TH1F (Form("htime1%i", i), "", 100, -6,6);
  //     htime2[i]=new TH1F (Form("htime2%i", i), "", 100, -6,6);
  //     htime1[i]->SetLineColor(kGreen);
  //     htime2[i]->SetLineColor(kRed);
  //   }

  PedestalSimple &pedsim=*prepPedSimple(nameped, roi);

  ifstream plik2(namecal.data(), ios::binary);
  if (!plik2.is_open())
    {
      cout<<"file: "<<namecal<<" is not open, exiting"<<endl;
      return;
    }
  Event evt2(roi);
  TimeCalClu caltime(1024, 8, 7);
  for (int ev=0; ev<19000; ev++)
    {
      if (ev % 10000 == 0)
   	cout<<"Event "<<ev<<endl;
      if (!evt2.read(plik2))
   	{
   	  cout<<"event "<<ev<<"end of file"<<endl;
   	  break;
   	}

      removePed(evt2, pedsim);
      evt2.CorrTime();
      evt2.InterpolatePseudoPulses();
      float times[nn][2];
      for (int i=0; i<nn; i++)
   	for (int j=0; j<2; j++)//gain
   	  {
   	    float sumhi = evt2.Sliding(i, j, windowsize, &times[i][j]);
   	    caltime.Fill(i, j, evt2.firstcap[i][j]%1024, times[i][j]);
   	  }
    }
  plik2.close();

  caltime.Finalize();
  caltime.Print();

  int nhist=20;
  double ranges[nhist+1];
  for (int i=0; i<=nhist; i++)
    ranges[i]=2*pow(100, 1.*i/nhist);
  TH1F *hists[nhist];
  for (int i=0; i<nhist; i++)
    {
      cout<<ranges[i]<<" - "<<ranges[i+1]<<endl;
      hists[i]=new TH1F (Form("htime%i", i), Form("Nphe = %.1f -%.1f", ranges[i], ranges[i+1]), 100, -6, 6);
    }

  // second loop over an independent file
  for (int filter=3; filter<=18; filter++)
    {
      TString namelin=namelin0;
      namelin.ReplaceAll("Filter0", Form("Filter%i", filter));
      ifstream plik3(namelin.Data(), ios::binary);
      if (!plik3.is_open())
	{
	  cout<<"file: "<<namelin<<" is not open, exiting"<<endl;
	  return;
	}
      Event evt3(roi);

      for (int ev=0; ev<2000; ev++)
	{
	  if (ev % 10000 == 0)
	    cout<<"Event "<<ev<<endl;
	  if (!evt3.read(plik3))
	    {
	      cout<<"event "<<ev<<"end of file"<<endl;
	      break;
	    }

	  removePed(evt3, pedsim);
	  evt3.CorrTime();
	  evt3.InterpolatePseudoPulses();
	  // for (int i=0; i<nn; i++)
	  int i=0;
	    {
	      float time;
	      float sumhi = evt3.Sliding(i, 0, windowsize, &time);
	      float nphe = sumhi*ff[i].conv;
	      time-=caltime.GetCorrTime(i, 0, evt3.firstcap[i][0]%1024);
	      for (int bin=0; bin<nhist; bin++)
		{
		  if (nphe>ranges[bin] && nphe<ranges[bin+1])
		    hists[bin]->Fill(time);
		}

	      // htime0[i]->Fill(time-caltime1.GetMeanTime(i, 0));
	      // 	  htime1[i]->Fill(time-caltime1.GetCorrTime(i, 0, evt3.firstcap[i][0]));
	      // 	  htime2[i]->Fill(time-caltime2.GetCorrTime(i, 0, evt3.firstcap[i][0]%1024));
	    }
	}
      plik3.close();
    }

  TGraphErrors *gtimerms = new TGraphErrors();
  for (int i=0; i<nhist; i++)
    if (hists[i]->GetEntries()>10)
      {
	gtimerms->SetPoint(gtimerms->GetN(), sqrt(ranges[i]*ranges[i+1]), hists[i]->GetRMS());
	gtimerms->SetPointError(gtimerms->GetN()-1, 0, hists[i]->GetRMSError());
      }

  TCanvas *c2 = new TCanvas ("c2", "", 800, 600);
  c2->SetLogx();
  gtimerms->Draw("APE");
  gtimerms->GetXaxis()->SetTitle("Nphe");
  gtimerms->GetYaxis()->SetTitle("Time resolution");
  TF1 *ftimefit = new TF1("ftimefit", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+[2]*[2])", 1, 300);
  ftimefit->SetParameters(1, 1, 1);
  gtimerms->Fit(ftimefit);

  TCanvas *c3 = new TCanvas ("c3", "", 800, 600);
  c3->Divide(5, 4, 0.001, 0.001);
  for (int i=0; i<nhist; i++)
    {
      c3->cd(i+1);
      hists[i]->Draw();
    }
  // c2->Divide(4,2, 0.001, 0.001);
  // for (int i=0; i<nn; i++)
  //   {
  //     c2->cd(i+1);
  //     htime1[i]->Draw();
  //     htime0[i]->Draw("sames");
  //     htime2[i]->Draw("sames");
  //   }
}


//==================================================================
// checking time resolution vs Nphe for all pixels
// nameped - 0th pedestal file
// namecal - 0th calibration file to get coeff
// namelin - 0th linearity files
void CheckTimeResVsPheAll(int roi=40)
{
  const int nclu2=19; // number of clusters
  int ips[nclu2]={237, 331, 76, 216, 25, 110, 116, 258, 62, 270, 263, 120, 315, 205, 65, 87, 84, 7, 257}; // IP numbers in file names
  // const int nclu=3; // for tests
  const int nclu=nclu2;

  const int nn=nch-1;
  int windowsize=6;

  string nameped0="../TimeResStudy/PedestalData/PedTable_20170901_160812_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB0_IP237.dat";
  string namecal0="../TimeResStudy/CalibrationData/DATA_Linearity_20170901_174631/Linearity_NominalHV_Filter6RD40_FEB0_IP237.dat";
  string namelin00="../TimeResStudy/LinearityData/DATA_Linearity_20170901_173134/Linearity_NominalHV_Filter0RD40_FEB0_IP237.dat";


  PedestalSimple *pedsims[nclu];
  FFactorStat ffs[nclu][nch];
  TimeCalClu *caltimes[nclu];
  for (int iclu=0; iclu<nclu; iclu++)
    {
      cout<<"=== processing cluster "<<iclu<<" ==="<<endl;
      TString nameped(nameped0);
      nameped.ReplaceAll("_FEB0_", Form("_FEB%i_", iclu));
      nameped.ReplaceAll(Form("_IP%i.dat", ips[0]), Form("_IP%i.dat", ips[iclu]));
      TString namecal(namecal0);
      namecal.ReplaceAll("_FEB0_", Form("_FEB%i_", iclu));
      namecal.ReplaceAll(Form("_IP%i.dat", ips[0]), Form("_IP%i.dat", ips[iclu]));

      // calibration
      ffactor(nameped.Data(), namecal.Data(), roi, ffs[iclu], 0, true, windowsize);

      // pedestals (again)
      pedsims[iclu]=prepPedSimple(nameped.Data(), roi);

      // time calibration
      caltimes[iclu]= new TimeCalClu(1024, 8, 7);

      ifstream plik2(namecal.Data(), ios::binary);
      if (!plik2.is_open())
	{
	  cout<<"file: "<<namecal<<" is not open, exiting"<<endl;
	  return;
	}
      Event evt2(roi);
      for (int ev=0; ev<19000; ev++)
	{
	  if (ev % 10000 == 0)
	    cout<<"Event "<<ev<<endl;
	  if (!evt2.read(plik2))
	    {
	      cout<<"event "<<ev<<"end of file"<<endl;
	      break;
	    }

	  removePed(evt2, *pedsims[iclu]);
	  evt2.CorrTime();
	  evt2.InterpolatePseudoPulses();
	  float times[nn][2];
	  for (int i=0; i<nn; i++)
	    for (int j=0; j<2; j++)//gain
	      {
		float sumhi = evt2.Sliding(i, j, windowsize, &times[i][j]);
		caltimes[iclu]->Fill(i, j, evt2.firstcap[i][j]%1024, times[i][j]);
	      }
	}
      plik2.close();
      caltimes[iclu]->Finalize();
      caltimes[iclu]->Print();
    }


  int nhist=20;
  double ranges[nhist+1];
  for (int i=0; i<=nhist; i++)
    ranges[i]=2*pow(100, 1.*i/nhist);
  TH1F *hists[nclu][nch][nhist];
  for (int iclu=0; iclu<nclu; iclu++)
    for (int ich=0; ich<nch; ich++)
      for (int i=0; i<nhist; i++)
	hists[iclu][ich][i]=new TH1F (Form("htime_%i_%i%i", iclu, ich, i), "", 100, -6, 6);

  TH1F *hpulsemean[nhist]; // distribution of the avarage pulse position (to check that it does not depend on the intensity)
  for (int i=0; i<=nhist; i++)
    hpulsemean[i]=new TH1F (Form("hpulsemean_%i", i), "", 100, -6, 6);

  // second loop, over linearity files
  cout<<"==="<<endl<<" Starting loop over linearity files "<<endl;
  for (int filter=3; filter<=18; filter++)
    {
      cout<<"filter "<<filter<<endl;
      TString namelin0=namelin00;
      namelin0.ReplaceAll("Filter0", Form("Filter%i", filter));

      Event *evts[nclu];
      ifstream *files[nclu];
      for (int iclu=0; iclu<nclu; iclu++)
	{
	  evts[iclu] = new Event(roi);

	  TString namelin(namelin0);
	  namelin.ReplaceAll("_FEB0_", Form("_FEB%i_", iclu));
	  namelin.ReplaceAll(Form("_IP%i.dat", ips[0]), Form("_IP%i.dat", ips[iclu]));

	  files[iclu] = new ifstream(namelin.Data(), ios::binary);
	  if (!files[iclu]->is_open())
	    {
	      cout<<"file: "<<namelin<<" is not open, exiting"<<endl;
	      return;
	    }
	}

      for (int ev=0; ev<2000; ev++)
	{
	  if (ev % 1000 == 0)
	    cout<<"Event "<<ev<<endl;
	  if (!constructevent(evts, files, nclu))
	    {
	      cout<<"event "<<ev<<"end of file"<<endl;
	      break;
	    }

	  TH1F htime("htime", "", 100, -6, 6); // all the arrival times from this event
	  TH1F hsignal ("hsignal", "", 100, -1, 3000); // all the nphe from this event
	  float times[nclu][nch]={0};
	  float nphes[nclu][nch]={0};
	  // the first loop is to get the average arrival time of this event (to correct jitters)
	  for (int iclu=0; iclu<nclu; iclu++)
	    {
	      removePed(*evts[iclu], *pedsims[iclu]);
	      evts[iclu]->CorrTime();
	      evts[iclu]->InterpolatePseudoPulses();
	      for (int ich=0; ich<nn; ich++)
		{
		  float sumhi = evts[iclu]->Sliding(ich, 0, windowsize, &times[iclu][ich]);
		  nphes[iclu][ich] = sumhi*ffs[iclu][ich].conv;
		  times[iclu][ich]-=caltimes[iclu]->GetCorrTime(ich, 0, evts[iclu]->firstcap[ich][0]%1024);
		  htime.Fill(times[iclu][ich]);
		  hsignal.Fill(nphes[iclu][ich]);
		}
	    }
	  double meantime= htime.GetMean();
	  double meanphe = hsignal.GetMean();
	  for (int bin=0; bin<nhist; bin++)
	    if (meanphe>ranges[bin] && meanphe<ranges[bin+1])
	      hpulsemean[bin]->Fill(meantime);

	  // the second loop fills the actual histograms after all the corrections
	  for (int iclu=0; iclu<nclu; iclu++)
	    for (int ich=0; ich<nn; ich++)
	      for (int bin=0; bin<nhist; bin++)
		if (nphes[iclu][ich]>ranges[bin] && nphes[iclu][ich]<ranges[bin+1])
		  hists[iclu][ich][bin]->Fill(times[iclu][ich]-meantime);
	}

      // now closing everything
      for (int iclu=0; iclu<nclu; iclu++)
	{
	  delete evts[iclu];
	  files[iclu]->close();
	  delete files[iclu];
	}
    }
  cout<<"Average pulse estimation: "<<endl;
  for (int bin=0; bin<nhist; bin++)
    cout<<ranges[bin]<<" - "<<ranges[bin+1]<<"phe: "<<hpulsemean[bin]->GetRMS()<<endl;

  TH2F *ht0vst1 = new TH2F("ht0vst1", ";T1 (reconstruction)[ns]; T0 (poissonian)[ns]", 40, 0, 3.5, 40, 0, 2);
  TH1F *ht2 = new TH1F ("ht2", ";T2 (constant) [ns]; # of pixels", 40, 0, 1);

  TH1F *hres5phe = new TH1F("hres5phe", ";Time resolution at 5 phe[ns]; # of channels", 40, 0, 2);

  // end of reading of files, now time to do graphs and fits
  TGraphErrors *gtimerms[nclu][nch];
  for (int iclu=0; iclu<nclu; iclu++)
    for (int ich=0; ich<nn; ich++)
      {
	gtimerms[iclu][ich] = new TGraphErrors();
	for (int i=0; i<nhist; i++)
	  if (hists[iclu][ich][i]->GetEntries()>10)
	    {
	      int ipoint=gtimerms[iclu][ich]->GetN();
	      gtimerms[iclu][ich]->SetPoint(ipoint, sqrt(ranges[i]*ranges[i+1]), hists[iclu][ich][i]->GetRMS());
	      gtimerms[iclu][ich]->SetPointError(ipoint, 0, hists[iclu][ich][i]->GetRMSError());
	    }
	gtimerms[iclu][ich]->GetXaxis()->SetTitle("Nphe");
	gtimerms[iclu][ich]->GetYaxis()->SetTitle("Time res. [ns]");

	TF1 *ftimefit = new TF1("ftimefit", "sqrt(pow([0]/sqrt(x),2)+pow([1]/x,2)+[2]*[2])", 1, 300);
	ftimefit->SetParameters(1, 1, 1);
	gtimerms[iclu][ich]->Fit(ftimefit, "Q");
	double t0=ftimefit->GetParameter(0);
	double t1=ftimefit->GetParameter(1);
	double t2=ftimefit->GetParameter(2);
	cout<<iclu<<"FEB, "<<ich<<"ch:  "<<t0<<" "<<t1<<" "<<t2<<endl;
	ht0vst1->Fill(t1, t0);
	ht2->Fill(t2);
	hres5phe->Fill(ftimefit->Eval(5.));
      }
  TCanvas *ct10 = new TCanvas ("ct10", "", 640, 480);
  ht0vst1->Draw("COLZ");
  TCanvas *ct2 = new TCanvas ("ct2", "", 640, 480);
  ht2->Draw();
  TCanvas *cres5 = new TCanvas ("cres5", "", 640, 480);
  hres5phe->Draw();

  TFile plikout("time_calib_out.root", "RECREATE");
  plikout.cd();
  for (int iclu=0; iclu<nclu; iclu++)
    for (int ich=0; ich<nn; ich++)
      gtimerms[iclu][ich]->Write(Form("gtimerms_%i_%i", iclu, ich));
  ct10->Write();
  ct2->Write();
  hres5phe->Write();
  plikout.Close();


  // TCanvas *c2 = new TCanvas ("c2", "", 800, 600);
  // c2->SetLogx();
  // gtimerms->Draw("APE");
  // gtimerms->GetXaxis()->SetTitle("Nphe");
  // gtimerms->GetYaxis()->SetTitle("Time resolution");

  // TCanvas *c3 = new TCanvas ("c3", "", 800, 600);
  // c3->Divide(5, 4, 0.001, 0.001);
  // for (int i=0; i<nhist; i++)
  //   {
  //     c3->cd(i+1);
  //     hists[i]->Draw();
  //   }
  // // c2->Divide(4,2, 0.001, 0.001);
  // // for (int i=0; i<nn; i++)
  // //   {
  // //     c2->cd(i+1);
  // //     htime1[i]->Draw();
  // //     htime0[i]->Draw("sames");
  // //     htime2[i]->Draw("sames");
  // //   }
}

void CheckTimeCalTaka()
{
  int roi=40;

  string nameped="../TimeResStudy/PedestalData/PedTable_20170901_160812_IC_ET_RD40/Ped_Random1000_IC_ET_RD40_FEB0_IP237.dat";
  string namecal5="../TimeResStudy/CalibrationData/DATA_Linearity_20170901_174631/Linearity_NominalHV_Filter5RD40_FEB0_IP237.dat";
  string namecal6="../TimeResStudy/CalibrationData/DATA_Linearity_20170901_174631/Linearity_NominalHV_Filter6RD40_FEB0_IP237.dat";
  string namecal7="../TimeResStudy/CalibrationData/DATA_Linearity_20170901_174631/Linearity_NominalHV_Filter7RD40_FEB0_IP237.dat";
  string namecal8="../TimeResStudy/CalibrationData/DATA_Linearity_20170901_174631/Linearity_NominalHV_Filter8RD40_FEB0_IP237.dat";
  string namecal9="../TimeResStudy/CalibrationData/DATA_Linearity_20170901_174631/Linearity_NominalHV_Filter9RD40_FEB0_IP237.dat";

  PedestalSimple &pedsim=*prepPedSimple(nameped, roi);

  TimeCalClu caltime(1024, 8, 7);
  caltime.FillFile(namecal6, pedsim, roi, 6);
  caltime.Finalize();


  // TimeTakaClu cal(4096);
  TimeTakaClu cal(1024);
  cal.FillFile(namecal5, pedsim, roi);
  cal.FillFile(namecal6, pedsim, roi);
  cal.FillFile(namecal7, pedsim, roi);
  cal.FillFile(namecal8, pedsim, roi);
  cal.FillFile(namecal9, pedsim, roi);
  cal.Finalize();


  TGraphErrors *gr = new TGraphErrors();
  gr->SetMarkerStyle(7);
  for (int i=0; i<1024; i++)
    {
      cout<<i<<" "<<cal.width[0][i]<<"+/-"<<cal.dwidth[0][i]<<endl;
      gr->SetPoint(i, i, cal.width[0][i]);
      gr->SetPointError(i, 0, cal.dwidth[0][i]);
    }
  TCanvas *c = new TCanvas ("c", "", 640, 480);
  c->cd();
  gr->Draw("APE");

  TGraph *gmarscorr[nch];
  TGraph *gtakacorr[nch];

  for (int ich=0; ich<7; ich++)
    {
      // int ich=1;
      int pulsepos=30;
      // now we do comparisons of both methods
      cout<<"finding the integration time in Taka's method"<<endl;
      int bestsize=-1;
      double bestdiff=1.e7;
      for (int size=0; size<1024; size++) // size of the integration
      	{
      	  double totaldiff=0;
      	  for (int pos0=0; pos0<1024; pos0++)
      	    {
      	      float marscorr = caltime.GetCorrTime(ich, 0, pos0%1024)-caltime.GetMeanTime(ich,0) ;
      	      float takacorr=0;
      	      for (int i=-pulsepos; i<size; i++)
      		{
      		  int pos=(pos0-i + 1024)%1024;
      		  takacorr+=cal.width[ich][pos]-1;
      		}
      	      float diff=takacorr-marscorr;
      	      // totaldiff+=diff*diff;
      	      totaldiff+=fabs(diff);
      	    }
      	  // totaldiff=sqrt(totaldiff/1024);
      	  totaldiff/=1024;
      	  if (totaldiff<bestdiff)
      	    {
      	      bestdiff=totaldiff;
      	      bestsize=size;
      	    }
      	  // cout<<"size="<<size<<", averagediff="<<totaldiff<<endl;
      	}
      cout<<"channel"<<ich<<", best size="<<bestsize<<", averagediff="<<bestdiff<<endl;

      // int bestsize=958; //fixed by hand

      gmarscorr[ich] = new TGraph();
      gtakacorr[ich] = new TGraph();
      for (int pos0=0; pos0<1024; pos0++)
	{
	  gmarscorr[ich]->SetPoint(pos0, pos0,  caltime.GetCorrTime(ich, 0, pos0%1024)-caltime.GetMeanTime(ich,0));
	  float takacorr=0;
	  for (int i=-pulsepos; i<bestsize; i++)
	    {
	      int pos=(pos0-i + 1024)%1024;
	      takacorr+=cal.width[ich][pos]-1;
	    }
	  gtakacorr[ich]->SetPoint(pos0, pos0,  takacorr);
	}
      gmarscorr[ich]->SetLineColor(kRed);
    }

  TCanvas *cc =new TCanvas ("cc", "", 640, 480);
  cc->Divide(4, 2, 0.001, 0.001);
  for (int i=0; i<7; i++)
    {
      cc->cd(i+1);
      gmarscorr[i]->Draw("AL");
      gtakacorr[i]->Draw("L");
    }

}

// macro for testing the Tsutomu pattern, i.e. dependence of the baseline position in the first 9 slices on the stopcell
void TestTsutomu(string nameped, string namein, int roi, int mypix, int mygain)
{
  cout<<"opening pedestal"<<nameped<<endl;
  ifstream plik1(nameped.data(), ios::binary);
  if (!plik1.is_open())
    {
      cout<<"file: "<<nameped<<" is not open, exiting"<<endl;
      return;
    }
  Event evt1(roi);

  PedestalSimple ped;
  // PedestalUltimate ped(roi);
  for (int ev=0; ev<20000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt1.read(plik1))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      evt1.CorrTime();
      // evt1.CorrTsutomu();

      ped.fillPedEvent(evt1);
    }
  plik1.close();
  ped.finalizePed();

  // second loop for Tsutomu
  PedestalTsutomu pedtsu;
  ifstream plik1b(nameped.data(), ios::binary);
  for (int ev=0; ev<20000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt1.read(plik1b))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      removePed(evt1, ped);
      evt1.CorrTime();

      pedtsu.fillPedEvent(evt1);
    }
  plik1b.close();
  pedtsu.finalizePed();

  TH1F *htime = new TH1F ("htime", ";timediff[ms]; Number of events", 100, 0, 10);


  ifstream plik2(namein.data(), ios::binary);
  if (!plik2.is_open())
    {
      cout<<"file: "<<namein<<" is not open, exiting"<<endl;
      return;
    }
  Event evt2(roi);

  const Int_t nbits=9;
  TH2F *hbits[nbits];
  for (int i=0; i<nbits; i++)
    {
      Int_t nc=1024/(1<<i);
      hbits[i] = new TH2F(Form("hbits_%i", i), Form(";stop cell%%%i ; signal[counts]", nc), nc, -0.5, nc-0.5, 100, -50, 50);
    }

  TH2F *h2 = new TH2F ("h2", ";sample in ROI; signal[counts]", roi, -0.5, roi-0.5, 100, -50, 50);
  TH2F *h3 = new TH2F ("h3", ";stop cell; signal[counts]", size4drs/4, -0.5, size4drs/4-0.5, 100, -50, 50);
  Double_t lasttime=-1;
  for (int ev=0; ev<20000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt2.read(plik2))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}


      int fc = evt2.firstcap[mypix][mygain];
      // for (int i=2; i<roi-2; i++)
      // 	cout<<evt2.GetDRSTime() - evt2.lasttime[mypix][mygain][(i+fc) % size4drs]<<", ";
      // cout<<endl;


      removePed(evt2, ped);
      evt2.CorrTime();
      // pedtsu.CorrEvt(evt2);
      evt2.InterpolatePseudoPulses();

      for (int i=0; i<roi; i++)
	h2->Fill(i, evt2.samples[mypix][mygain][i]);
      h3->Fill(fc%1024, evt2.samples[mypix][mygain][2]);

      for (int i=0; i<nbits; i++)
	{
	  Int_t nc=1024/(1<<i);
	  hbits[i]->Fill(fc%nc, evt2.samples[mypix][mygain][i]);
	}

    }
  plik2.close();

  TCanvas *cc = new TCanvas ("cc", "", 640, 480);
  h2->Draw("colz");
  h2->ProfileX("_pfx", 1, -1, "s")->Draw("same");
  TCanvas *cc2 = new TCanvas ("cc2", "", 640, 480);
  h3->Draw("colz");
  TCanvas *cc3 = new TCanvas ("cc3", "", 640, 480);
  cc3->Divide(3,3,0.001, 0.001);
  for (int i=0; i<nbits; i++)
    {
      cc3->cd(i+1);
      hbits[i]->Draw("colz");
    }
}


// makes pedestal distribution of all the pixels in one cluster after all the corrections
void PedestalDist(string nameped, string namein, int roi)
{
  const int window=6;
  cout<<"opening pedestal"<<nameped<<endl;
  ifstream plik1(nameped.data(), ios::binary);
  if (!plik1.is_open())
    {
      cout<<"file: "<<nameped<<" is not open, exiting"<<endl;
      return;
    }
  Event evt1(roi);

  PedestalSimple ped;
  // PedestalUltimate ped(roi);
  for (int ev=0; ev<20000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt1.read(plik1))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}
      cout<<"Event: "<<ev<<endl;
      evt1.CorrTime();
      // evt1.CorrTsutomu();

      ped.fillPedEvent(evt1);
    }
  plik1.close();
  ped.finalizePed();

  // second loop for Tsutomu
  PedestalTsutomu pedtsu;
  ifstream plik1b(nameped.data(), ios::binary);
  for (int ev=0; ev<20000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt1.read(plik1b))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      removePed(evt1, ped);
      evt1.CorrTime();

      pedtsu.fillPedEvent(evt1);
    }
  plik1b.close();
  pedtsu.finalizePed();


  ifstream plik2(namein.data(), ios::binary);
  if (!plik2.is_open())
    {
      cout<<"file: "<<namein<<" is not open, exiting"<<endl;
      return;
    }
  Event evt2(roi);

  const int nn=7;
  TH1F *hped[nn][2];
  // TH1F *hpeds[nn][2];
  for (int i=0; i<nn; i++)
    for (int j=0; j<2; j++)
      {
	hped[i][j] = new TH1F(Form("hped_%i_%s", i, (j==0)?"hi":"lo"), ";signal[counts];number of entries", 4000, -100-0.5, 3900-0.5);
	// hpeds[i][j] = new TH1F(Form("hpeds_%i_%s", i, (j==0)?"hi":"lo"), ";signal[counts];number of entries", 4000, -100, 3900);
	// hpeds[i][j]->SetLineColor(kRed);
	if (window>1)
	  hped[i][j]->GetXaxis()->SetTitle(Form("signal/sqrt(%i slices) [counts]", window));
      }

  Double_t lasttime=-1;
  for (int ev=0; ev<20000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt2.read(plik2))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}



      removePed(evt2, ped);
      evt2.CorrTime();
      pedtsu.CorrEvt(evt2);
      evt2.InterpolatePseudoPulses();
      if (ev>0)
	{
	  for (int i=0; i<nn; i++)
	    for (int j=0; j<2; j++)
	      {
		for (int k=3; k<roi-3; k++)
		  // for (int k=0; k<1; k++)
		  // if ((evt2.firstcap[i][j]+k)%32 != 31)
		  // if ((evt2.firstcap[i][j]+k+1)%size4drs != evt2.oldfirstcap[i][j]) // spike B
		  //   hped[i][j]->Fill(evt2.samples[i][j][k]);
		  // else
		  //   hpeds[i][j]->Fill(evt2.samples[i][j][k]);
		  for (int k=3; k<evt2.roisize-3-window; k++)  // skipped the first and the last slices
		    hped[i][j]->Fill(evt2.SumSlices(i, j, k, window)/sqrt(1.*window));
        cout<<"hped"<<hped[i][j] << endl;
		}

	}
    }
  plik2.close();

  for (int j=0; j<2; j++)
    {
      Double_t meanrms=0;
      for (int i=0; i<nn; i++)
	{
	  cout<<"Pix "<<i<<(j==0?"hi":"lo")<<", RMS= "<<hped[i][j]->GetRMS()<<", Mean="<<hped[i][j]->GetMean()<<endl;
	  meanrms+=hped[i][j]->GetRMS();
	}
      meanrms/=nn;
      cout<<"Average RMS for "<<(j==0?"hi":"lo")<<" gain="<<meanrms<<endl;
    }


  TCanvas *cchi = new TCanvas ("cchi", "", 1024, 768);
  cchi->Divide(4,2,0.001, 0.001);
  for (int i=0; i<nn; i++)
    {
      cchi->cd(i+1);
      cchi->GetPad(i+1)->SetLogy();
      hped[i][0]->Draw();
      hped[i][0]->GetXaxis()->SetRangeUser(hped[i][0]->GetMean()-100, hped[i][0]->GetMean()+100);
      // hpeds[i][0]->Draw("sames");
    }

  TCanvas *cclo = new TCanvas ("cclo", "", 1024, 768);
  cclo->Divide(4,2,0.001, 0.001);
  for (int i=0; i<nn; i++)
    {
      cclo->cd(i+1);
      cclo->GetPad(i+1)->SetLogy();
      hped[i][1]->Draw();
      hped[i][1]->GetXaxis()->SetRangeUser(hped[i][1]->GetMean()-100, hped[i][1]->GetMean()+100);
    }
}

int main(int argv, char **argc)
{
  readdragon();
  return 0;
}
