#ifndef ffactor_H
#define ffactor_H

#include <cmath>

#include "event.h"
#include "pixstat.h"
#include "pedestal.h"

class FFactorStat{
public:
  float f2fact; // assumed F-factor
  float meansig; // mean signal in nphe phe
  float rmssig; // rms of the above signal
  bool isfit;  // flag if the fit was used (not used) 
  bool isgood; // flag if the pixel is good
  float pedrms; // pedestal RMS
  float conv; // cnts to phe conversion
  float nphe; // number of phe
  float peak; // peak signal of 1 phe, not necessarily filled

  FFactorStat() : f2fact(1.2), meansig(0), rmssig(0), isfit(true), isgood(true), pedrms(true), conv(0), nphe(0), peak(0) {}
  bool doFFactor()
  {
    isgood=true;
    if (meansig<5*pedrms)
      isgood=false;
    nphe = f2fact*meansig*meansig /(rmssig*rmssig - pedrms*pedrms);
    conv = nphe/meansig;
    return isgood;
  }
  void print()
  {
    cout<<(isgood?"GOOD":"BAD")<<"<signal>="<<meansig<<", Nphe="<<nphe<<", conv = "<<conv<<", noise="<<pedrms<<"cnt="<<pedrms*conv<<"phe"<<endl;
  }
};

PedStat* getrms (string nameped, int roi, PedestalSimple* pedsim, int windowsize=6, int nevmax=1000)
{
  PedStat *pedstat = new PedStat();
  ifstream plik(nameped.data(), ios::binary);
  if (!plik.is_open())
    {
      cout<<"file: "<<nameped<<" is not open, exiting"<<endl;
      return 0;
    }
  Event evt(roi);

  const int nn=nch-1;

  int nev=0;
  for (int i=0; i<nevmax; i++)
    {
      if (!evt.read(plik))
	break;

      removePed(evt, *pedsim);
      evt.CorrTime();
      evt.InterpolatePseudoPulses();

      for (int i=0; i<nn; i++)
	for (int j=0; j<2; j++)
	  {
	    float sum = evt.Sliding(i, j, windowsize);
	    pedstat->fmean[i][j]+=sum;
	    pedstat->frms[i][j]+=sum*sum;
	  }
	nev++;
    }
  if (nev>0)
    for (int i=0; i<nn; i++)
      for (int j=0; j<2; j++)
	{
	  pedstat->fmean[i][j]/=nev;
	  pedstat->frms[i][j]/=nev;
	  pedstat->frms[i][j]=sqrt(pedstat->frms[i][j] - pedstat->fmean[i][j]*pedstat->fmean[i][j]);
	}
  plik.close();
  return pedstat;
}

#endif 
