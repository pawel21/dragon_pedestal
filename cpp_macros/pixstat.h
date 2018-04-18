#ifndef pixstat_H
#define pixstat_H

#include "event.h"
// container to store pedestal statistics 
class PedStat
{
public:
  double fmean[nch-1][2]; // mean pedestal for all pixels for high and low gain
  double frms[nch-1][2]; // rms of the pedestal for all pixels for high and low gain
  void print()
  {
    cout<<"Mean (high gain): ";
    for (int i=0; i<nch-1; i++) cout<<fmean[i][0]<<", ";
    cout<<endl;
    cout<<"Mean (low gain): ";
    for (int i=0; i<nch-1; i++) cout<<fmean[i][1]<<", ";
    cout<<endl;
    cout<<"RMS (high gain): ";
    for (int i=0; i<nch-1; i++) cout<<frms[i][0]<<", ";
    cout<<endl;
    cout<<"RMS (low gain): ";
    for (int i=0; i<nch-1; i++) cout<<frms[i][1]<<", ";
    cout<<endl;
  }
};


#endif 
