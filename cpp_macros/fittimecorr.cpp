#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TF1.h>
#include <TProfile.h>

#include "event.h"

using namespace std;

TH1F *h1, *h2, *h3, *h4; 

void makeone (TH2F *h, string namein, Int_t roi = 40)
{
  cout<<"opening "<<namein<<endl;
  ifstream plik(namein.data(), ios::binary);
  if (!plik.is_open())
    {
      cout<<"file: "<<namein<<" is not open, exiting"<<endl;
      return;
    }
  Event evt(roi);

  const int npix=nch-1;
  double ped0[npix][2][size4drs]={{{0}}};
  int ped0n[npix][2][size4drs]={{{0}}};

  // first loop to get rough pedestal from long delays
  for (int ev=0; ev<9000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt.read(plik))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      double timenow=evt.GetDRSTime(); // [ms]
      for (int ipix=0; ipix<npix; ipix++)
	for (int igain=0; igain<2; igain++)
	  {
	    int fc=evt.firstcap[ipix][igain];

	    for (int i=10; i<roi-2; i++) // 2 & 2 capacitors dropped !
	      {
		int posabs = (i+fc)%size4drs;
		if ((evt.lasttime[ipix][igain][posabs]>0)&&(timenow-evt.lasttime[ipix][igain][posabs]>20)) // > 20 ms
		  {
		    ped0[ipix][igain][posabs]+=evt.samples[ipix][igain][i];
		    ped0n[ipix][igain][posabs]++;
		  }
	      }
	  }
      evt.CorrTime(); // just to reset the times 
    }
  for (int ipix=0; ipix<npix; ipix++)
    for (int igain=0; igain<2; igain++)
      for (int i=0; i<size4drs; i++)
	ped0[ipix][igain][i]/=ped0n[ipix][igain][i];
  
  // second loop to get the time correction function
  plik.seekg(0); // rewind file


  double timeold=-1;
  for (int ev=0; ev<90000; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt.read(plik))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}

      // apply pedestal
      for (int ipix=0; ipix<npix; ipix++)
	for (int igain=0; igain<2; igain++)
	  {
	    int fc=evt.firstcap[ipix][igain];
	    for (int i=0; i<roi; i++) 
	      evt.samples[ipix][igain][i]-=ped0[ipix][igain][(i+fc)%size4drs];
	  }

      evt.InterpolatePseudoPulses();

      double timenow=evt.GetDRSTime(); // [ms]
      if (timenow-timeold<0)
      cout<<"time diff="<<timenow-timeold<<endl;
      timeold=timenow;
      for (int ipix=0; ipix<npix; ipix++)
	if (ipix%2==1)
	for (int igain=0; igain<2; igain++)
	  {
	    int fc=evt.firstcap[ipix][igain];

	    for (int i=3; i<roi-3; i++) // 2 & 2 capacitor dropped !
	      {
		int posabs = (i+fc)%size4drs;
		if (evt.lasttime[ipix][igain][posabs]>0) 
		  {
		    double timediff=timenow - evt.lasttime[ipix][igain][posabs];
		    float val = evt.samples[ipix][igain][i];
		    //if (posabs %32 !=31 ) 
		    // if (posabs!=(evt.oldoldoldfirstcap[ipix][igain]+roi-1)%size4drs) 
		    // if (posabs!=(evt.oldoldfirstcap[ipix][igain]+roi-1)%size4drs) 
			// if (evt.oldfirstcap[ipix][igain]%1024>767)
		    // if (posabs!=(evt.oldfirstcap[ipix][igain]+roi-1)%size4drs) //Seiya's capacitor
		      // if (posabs%1024<=766)
		      // if (posabs%1024==1)
		      // if (evt.samplestag[ipix][igain][posabs]==1)
		      {
			if (timediff>0)
			  {
			    // cout<<posabs<<" "<<evt.oldfirstcap[ipix][igain]<<endl;
			    h->Fill(log10(timediff), val);
			    if (log10(timediff)<-2)
			      cout<<ipix<<" "<<igain<<" "<<val<<" "<<timediff<<" "<<i<<" "<<posabs<<" "<<evt.samplestag[ipix][igain][posabs]<<endl;
			    if ((timediff>0) && (timediff<0.1))
				if (val<20)
				  {
				    // cout<<timediff<<" "<<val<<" "<<ipix<<" "<<igain<<" "<<fc<<" "<<posabs<<" "<<evt.oldfirstcap[ipix][igain]<<endl;
				    // h1->Fill((evt.oldoldfirstcap[ipix][igain]-evt.firstcap[ipix][igain]+size4drs)%size4drs);
				    h1->Fill((posabs - evt.oldoldoldfirstcap[ipix][igain]+size4drs)%size4drs);
				    // h1->Fill((posabs - evt.oldfirstcap[ipix][igain]+size4drs)%size4drs);
				    // h2->Fill((fc + evt.oldfirstcap[ipix][igain]+size4drs+i)%size4drs);
				    h3->Fill((evt.oldfirstcap[ipix][igain])%1024);
				    h4->Fill((fc + evt.oldfirstcap[ipix][igain]+2*size4drs-i)%size4drs);
				    // h2->Fill((evt.oldoldfirstcap[ipix][igain]-evt.firstcap[ipix][igain]+size4drs)%size4drs);
				  }
				// else
				    // h2->Fill((evt.oldoldfirstcap[ipix][igain]-evt.firstcap[ipix][igain]+size4drs)%size4drs);
				    h2->Fill((posabs-evt.oldoldfirstcap[ipix][igain]+size4drs)%size4drs);
				//   h2->Fill((posabs - evt.oldfirstcap[ipix][igain]+size4drs)%size4drs);
			    
			  }
		      }
		  }
	      }
	  }
      evt.CorrTime(); // just to reset the times 
    }
  plik.close();
}


void fittimecorr(Int_t roi=40) // string namein="data/Ext160909_1405.dat", 
{
  // a few test histograms
  h1 = new TH1F ("h1", "", 4096, -0.5, 4096-0.5);
  h2 = new TH1F ("h2", "", 4096, -0.5, 4096-0.5);
  h3 = new TH1F ("h3", "", 4096, -0.5, 4096-0.5);
  h4 = new TH1F ("h4", "", 4096, -0.5, 4096-0.5);


  TH2F *h= new TH2F ("h", "", 100, -2.5, 3.5, 100, -100, 250);
  const int nnames=19;

  string names[nnames]= {
    // "Randome5kHz20kev_calp800caln800_run1.dat",
    // "Randome5kHz20kev_calp800caln800_run2.dat",
    // "Randome5kHz20kev_run1.dat",
    // "Randome5kHz20kev_run2.dat",
    "Randome7kHz20kev_run1.dat",
    "Randome7kHz20kev_run2.dat",
// "../data/Ped_Random1000_IC_ET_RD40_FEB0_IP211.dat",
// "../data/Ped_Random1000_IC_ET_RD40_FEB1_IP217.dat",
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

// string names[nnames]= {  
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

  // for (int i=0; i<nnames; i++)
  // int i=2;
  for (int i=0; i<2; i++)
    makeone (h, names[i], roi);

  // TF1 *f = new TF1("f", "117* pow(x,-0.2255) -47.8 + [0]", 0.02, 200.);
  TF1 *f0 = new TF1("f0", "[1]* pow(pow(10,x),-[2])+ [0]", -2.4, 3.4);
  f0->SetParameters(-11.95, 28.98, 0.2405); // I guess some old curve
  f0->SetLineColor(kBlack);

  TF1 *f1 = new TF1("f1", "[1]* pow(pow(10,x),-[2])+ [0]", -2.4, 3.4);
  f1->SetParameters(-8, 23, 0.2405); // curve currently in event.h
  f1->SetLineColor(kMagenta);


  TF1 *f2 = new TF1("f2", "[1]* pow(pow(10,x),-[2])+ [0]", -2.4, 3.4);
  f2->SetParameters(-8*0.5, 23*0.5, 0.2405); // half of the curve currently in event.h
  f2->SetLineColor(kMagenta);
  f2->SetLineStyle(kDashed);

  TF1 *f = new TF1("f", "[1]* pow(pow(10,x),-[2])+ [0]", -2.4, 3.4);

  TCanvas *c = new TCanvas ("c", "", 640, 480);
  TH1F *hosie = c->DrawFrame(0.01, -50, 200, 150);
  c->cd();
  c->SetLogz();
  h->GetYaxis()->SetRangeUser(-50,150);
  h->Draw("colz");
  TH1F *hprf = (TH1F*)h->ProfileX();
  hprf->Draw("same");
  // cout<<"fitting ..."<<endl;
  // hprf->Fit(f, "R");
  // f0->Draw("same");
  f1->Draw("same");
  f2->Draw("same");

  TCanvas *c2 = new TCanvas ("c2", "", 800, 600);
  c2->Divide(2,2,0.001, 0.001);
  c2->cd(1);
  h1->Draw();
  c2->cd(2);
  h2->Draw();
  c2->cd(3);
  h3->Draw();
  c2->cd(4);
  h4->Draw();
}
