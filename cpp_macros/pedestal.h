#ifndef pedestal_H
#define pedestal_H

#include "event.h"

// pedestal as a function of an absolute position in DRS4
class PedestalSimple 
{
public: 
  double meanped[nch][2][size4drs];
  double rmsped[nch][2][size4drs];
  int numped[nch][2][size4drs];

  PedestalSimple()
  {
    for (int i=0; i<nch; i++)
      for (int j=0; j<2; j++)
	{
	  for (int k=0; k<size4drs; k++)
	    {
	      meanped[i][j][k]=0;
	      rmsped[i][j][k]=0;
	      numped[i][j][k]=0;
	    }
	}
  }
  
  void fillPedEvent(Event &evt)
  {
    const int nbits=9; 
    for (int i=0; i<nch; i++)
      for (int j=0; j<2; j++)
	{
	  int fc=evt.firstcap[i][j];
	  for (int k=2; k<evt.roisize-2; k++)  // skipped first and last 2 slices
	    {
	      int posabs = (k+fc)%size4drs;
	      unsigned short val = evt.samples[i][j][k];
	      
	      meanped[i][j][posabs]+=val;
	      rmsped[i][j][posabs]+=val*val;
	      numped[i][j][posabs]++;
	    }
	}
  }

  void finalizePed()
  {
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      {
	for (int k=0; k<size4drs; k++)
	  if (numped[i][j][k])
	    {
	      meanped[i][j][k]/=numped[i][j][k];
	      rmsped[i][j][k]/=numped[i][j][k];
	      rmsped[i][j][k]=sqrt(rmsped[i][j][k]-meanped[i][j][k]*meanped[i][j][k]);
	      // cout<<i<<" "<<j<<" "<<k<<" "<<meanped[i][j][k]<<endl;
#ifdef USEINTEGER
	      meanped[i][j][k]=round(meanped[i][j][k]);
#endif
	    }
	  else
	    cout<<"ERROR !!! Pix"<<i<<" gain: "<<j<<", capacitor "<<k<<" not enough events"<<endl;
      }
  }
};


// pedestal as a function of absolute and relative position in DRS4
class PedestalUltimate 
{
public: 
  double *meanped[nch][2][size4drs];
  double *rmsped[nch][2][size4drs];
  int *numped[nch][2][size4drs];
  int roisize;

  PedestalUltimate(int roi): roisize(roi)
  {
    for (int i=0; i<nch; i++)
      for (int j=0; j<2; j++)
	for (int k=0; k<size4drs; k++)
	  {
	    meanped[i][j][k]=new double[roisize];
	    rmsped[i][j][k]=new double[roisize];
	    numped[i][j][k]=new int[roisize];
	    for (int l=0; l<roisize; l++)
	      {
		meanped[i][j][k][l]=0;
		rmsped[i][j][k][l]=0;
		numped[i][j][k][l]=0;
	      }
	  }
  }

  ~PedestalUltimate()
  {
    for (int i=0; i<nch; i++)
      for (int j=0; j<2; j++)
	for (int k=0; k<size4drs; k++)
	  {
	    delete []meanped[i][j][k];
	    delete []rmsped[i][j][k];
	    delete []numped[i][j][k];
	  }
  }

  void fillPedEvent(Event &evt)
  {
    for (int i=0; i<nch; i++)
      for (int j=0; j<2; j++)
	{
	  int fc=evt.firstcap[i][j];
	  for (int k=0; k<evt.roisize; k++) 
	    {
	      int posabs = (k+fc)%size4drs;
	      unsigned short val = evt.samples[i][j][k];
	      
	      meanped[i][j][posabs][k]+=val;
	      rmsped[i][j][posabs][k]+=val*val;
	      numped[i][j][posabs][k]++;
	    }
	}
  }

  void finalizePed()
  {
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<size4drs; k++)
	for (int l=0; l<roisize; l++)
	  if (numped[i][j][k][l])
	    {
	      meanped[i][j][k][l]/=numped[i][j][k][l];
	      rmsped[i][j][k][l]/=numped[i][j][k][l];
	      rmsped[i][j][k][l]=sqrt(rmsped[i][j][k][l]-meanped[i][j][k][l]*meanped[i][j][k][l]);
	    }
	  else
	    cout<<"ERROR !!! Pix"<<i<<" gain: "<<j<<", capacitor "<<k<<"(absolute)"
		<<l<<"(relative) not enough events"<<endl;
  }

};

void removePed(Event &evt, PedestalSimple &ped)
{
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      {
	int fc=evt.firstcap[i][j];
	for (int k=0; k<evt.roisize; k++)  
	    evt.samples[i][j][k]-=ped.meanped[i][j][(k+fc)%size4drs];
      }
}

void removePed(Event &evt, PedestalUltimate &ped)
{
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      {
	int fc=evt.firstcap[i][j];
	for (int k=0; k<evt.roisize; k++)  
	  evt.samples[i][j][k]-=ped.meanped[i][j][(k+fc)%size4drs][k];
      }
}

PedestalSimple* prepPedSimple(string nameped, int roi, int numev=10000)
{
  cout<<"opening pedestal"<<nameped<<endl;
  ifstream plik1(nameped.data(), ios::binary);
  if (!plik1.is_open())
    {
      cout<<"file: "<<nameped<<" is not open, exiting"<<endl;
      return 0;
    }
  Event evt1(roi);

  PedestalSimple *pedsim = new PedestalSimple();
  cout<<"preparing pedestal from "<<nameped<<endl;
  for (int ev=0; ev<numev; ev++)
    {
      if (ev % 10000 == 0)
	cout<<"Event "<<ev<<endl;
      if (!evt1.read(plik1))
	{
	  cout<<"event "<<ev<<"end of file"<<endl;
	  break;
	}
      evt1.CorrTime();
      pedsim->fillPedEvent(evt1);
    }
  plik1.close();
  pedsim->finalizePed();
  return pedsim;
}

class PedestalTsutomu
{
public:
  double tsutomujump[nch][2];
  double tsutomubase[nch][2];
  int ntsutomujump[nch][2];
  int ntsutomubase[nch][2];

  PedestalTsutomu()
  {
    for (int i=0; i<nch; i++)
      for (int j=0; j<2; j++)
	{
	  tsutomujump[i][j]=0;
	  tsutomubase[i][j]=0;
	  ntsutomujump[i][j]=0;
	  ntsutomubase[i][j]=0;	  
	}
  }

  void fillPedEvent(Event &evt)
  {
    const int nbits=9; 
    for (int i=0; i<nch; i++)
      for (int j=0; j<2; j++)
	{
	  int fc=evt.firstcap[i][j];
	  for (int k=0; k<nbits; k++)  
	    {
	      int stepsize=1024/(1<<(k+2));
	      fc%=(4*stepsize);
	      if ((fc>=stepsize) && (fc<2*stepsize))
		{
		  tsutomujump[i][j]+=evt.samples[i][j][k];
		  ntsutomujump[i][j]++;		      
		}
	      else
		{
		  tsutomubase[i][j]+=evt.samples[i][j][k];
		  ntsutomubase[i][j]++;		      
		}
	    }	  
	}
  }
  void finalizePed()
  {
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      {
	if (ntsutomujump[i][j])
	  tsutomujump[i][j]/=ntsutomujump[i][j];
	if (ntsutomubase[i][j])
	  tsutomubase[i][j]/=ntsutomubase[i][j];
	tsutomujump[i][j]-=tsutomubase[i][j];
#ifdef USEINTEGER
	tsutomujump[i][j]=round(tsutomujump[i][j]);
#endif

	cout<<"pix "<<i<<" in gain "<<j<<", jump = "<<tsutomujump[i][j]<<endl;
      }
  }

  // correction of the Tsutomu pattern: shift by ~10 counts depending on the binary patern of the first capacitor
  void CorrEvt(Event &evt)
  {
    const int nbits=9; 
    for (int i=0; i<nch; i++)
      for (int j=0; j<2; j++)
	{
	  float shift = tsutomujump[i][j];
	  int fc=evt.firstcap[i][j];
	  for (int k=0; k<nbits; k++)  
	    {
	      int stepsize=1024/(1<<(k+2));
	      fc%=(4*stepsize);
	      if ((fc>=stepsize) && (fc<2*stepsize))
		evt.samples[i][j][k]-=shift; 
	    }
      }
}

};

#endif
