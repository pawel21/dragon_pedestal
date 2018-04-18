// methods for the calibration of the time response
#ifndef timecal_H
#define timecal_H
#include "event.h"
#include "pedestal.h"

#include <iostream>
#include <cmath>

// class for performind and holding DRS4 time calibration for a single cluster
class TimeCalClu
{
private:
  int fNumCap;          // number of capacitors (1024 or 4096)
  int fNumCombine;      // how many capacitors are combines in a single point in above graphs
  int fNumHarmonics;    // how many harmonics are computed, total number of coefficients is 1+fNumHarmonics*2
  int fNumPoints;       // how many points will be in final graph
  float *fMeanVal[nch-1][2]; // mean values per bin
  int *fNumMean[nch-1][2]; // entries per bin

public: 
  float *fan[nch-1][2]; // cos coeff
  float *fbn[nch-1][2]; // sin coeff
  bool isgood[nch-1][2]; // 

  TimeCalClu(int ncap, int ncomb, int nharm):fNumCap(ncap), fNumCombine(ncomb), fNumHarmonics(nharm)
  {
    fNumPoints=(fNumCap-1)/fNumCombine + 1;
    for (int i=0; i<nch-1; i++)
      for (int j=0; j<2; j++)
	{
	  fan[i][j]=new float[fNumHarmonics+1];
	  fbn[i][j]=new float[fNumHarmonics+1];
	  fMeanVal[i][j]=new float[fNumPoints];
	  fNumMean[i][j]=new int[fNumPoints];
	  for (int k=0; k<fNumPoints; k++)
	    {
	      fMeanVal[i][j][k]=0; 
	      fNumMean[i][j][k]=0; 
	    }
	}
  }
  
  ~TimeCalClu()
  {
    for (int i=0; i<nch-1; i++)
      for (int j=0; j<2; j++)
	{
	  delete []fan[i][j];
	  delete []fbn[i][j];
	  delete []fMeanVal[i][j];
	  delete []fNumMean[i][j];
	}
  }

  void Fill(int pixid, int gain, int firstcap, float arrslice)
  {
    int bin=firstcap/fNumCombine;
    if (bin>=fNumPoints)
      std::cout<<"Time calib ERROR, firstcap="<<firstcap<<" beyond the table range !" << std::endl;
    else
      {
	fMeanVal[pixid][gain][bin]+=arrslice;
	fNumMean[pixid][gain][bin]++;
      }    
  }

  bool FillFile(string namecal, PedestalSimple &pedsim, int roi, int windowsize=6)
  {
    const int nn=nch-1;
    ifstream plik2(namecal.data(), ios::binary);
    if (!plik2.is_open())
    {
      cout<<"file: "<<namecal<<" is not open, exiting"<<endl;
      return false;
    }

    Event evt2(roi);
    int ev=0; 
    while (evt2.read(plik2))
      {
	if (ev % 10000 == 0)
	  cout<<"Event "<<ev<<endl;
	ev++;
	
	removePed(evt2, pedsim);
	evt2.CorrTime();
	evt2.InterpolatePseudoPulses();
	float times[nn][2];
	for (int i=0; i<nn; i++)
	  for (int j=0; j<2; j++)//gain
	    {
	      float sumhi = evt2.Sliding(i, j, windowsize, &times[i][j]);
	      Fill(i, j, evt2.firstcap[i][j]%1024, times[i][j]);
	    }
      }  
    cout<<ev<<" events read"<<endl;
    plik2.close();
    return true;    
  }

  void IntegrateWithTrig (float *x, float *y, int n, float &an, float &bn)
  {
    float suma=0, sumb=0;
    for (int i=0; i<fNumPoints; i++)
    {
      suma+=y[i]*(float)fNumCombine*cos(2*M_PI*n*(x[i]/(float)fNumCap));
      sumb+=y[i]*(float)fNumCombine*sin(2*M_PI*n*(x[i]/(float)fNumCap));
    }
    an=suma*(2./(fNumPoints*fNumCombine));
    bn=sumb*(2./(fNumPoints*fNumCombine));
    //   cout<<"n="<<n<<", an="<<an<<", bn="<<bn<<", a0="<<sqrt(an*an+bn*bn)<<endl;
  }

  void Finalize()
  {
    // X table
    float *pos = new float [fNumPoints];
    for (int i=0; i<fNumPoints; i++)
      pos[i]=(i+0.5)*fNumCombine;

    // now we are checking if there are any empty bins
    for (int pixid=0; pixid<nch-1; pixid++)
      for (int gain=0; gain<2; gain++)
	{
	  isgood[pixid][gain]=true;
	  for (int i=0; i<fNumPoints; i++)
	    if (fNumMean[pixid][gain][i]==0)
	      {
		int prevbin=(i==0)?fNumPoints-1:i-1;
		int nextbin=(i==fNumPoints-1)?0:i+1;
		
		// if there are events in both of the surrounding bins
		if (fNumMean[pixid][gain][prevbin] * fNumMean[pixid][gain][nextbin] > 0)
		  {
		    std::cout<<"Pixel "<<pixid<<", gain:"<<gain<<", no entries in time calibration bin "<<i<<", will interpolate from bins "<<prevbin<<" and "<<nextbin<<endl;
		    fMeanVal[pixid][gain][i] = fMeanVal[pixid][gain][prevbin] + fMeanVal[pixid][gain][nextbin];
		    fNumMean[pixid][gain][i] = fNumMean[pixid][gain][prevbin] + fNumMean[pixid][gain][nextbin];
		  } 
		else
		  {
		    isgood[pixid][gain]=false;
		    std::cout<<"Pixel "<<pixid<<", gain: "<<gain<<", interpolation failed in bin "<<i<<" ! "<<endl;
		  }
	      }    
	}
  
    for (int pixid=0; pixid<nch-1; pixid++)
      for (int gain=0; gain<2; gain++)
	if (isgood[pixid][gain])
	  {
	    for (int i=0; i<fNumPoints; i++)
	      // if (fNumMean[pixid][i]) // already checked if pixel is good
	      fMeanVal[pixid][gain][i]/=fNumMean[pixid][gain][i]; // mean
	    
	    //expanding into Fourier series	  
	    for (int in=0; in<=fNumHarmonics; in++)
	      IntegrateWithTrig(pos, fMeanVal[pixid][gain], in, fan[pixid][gain][in], fbn[pixid][gain][in]);	    
	  }	  

    delete []pos;
  }
  void Print()
  {
    for (int pixid=0; pixid<nch-1; pixid++)
      for (int gain=0; gain<2; gain++)
	{
	  cout<<"Pix "<<pixid<<((gain==0)?"Hi":"Lo")<<endl;
	  for (int i=0; i<=fNumHarmonics; i++)
	    cout<<"a"<<i<<"="<<fan[pixid][gain][i]<<", b"<<i<<"="<<fbn[pixid][gain][i]<<endl;
	}
  }
  
  float GetCorrTime(int pixid, int gain, int firstcap)
  {
    float time=fan[pixid][gain][0]/2.;
    for (int in=1; in<=fNumHarmonics; in++)
      {
	time+=fan[pixid][gain][in] * cos( firstcap * (float)in*2*M_PI/(float)fNumCap );
	time+=fbn[pixid][gain][in] * sin( firstcap * (float)in*2*M_PI/(float)fNumCap );
      }
    return time; 
  }

  float GetMeanTime(int pixid, int gain)
  {
    return fan[pixid][gain][0]/2.;
  }

// returns fourier expansion amplitudes added in quadrature
  float GetSumAmp(int pixid, int gain)
  {
    float sum=0;
    for (int in=1; in<=fNumHarmonics; in++)
      sum+=fan[pixid][gain][0]*fan[pixid][gain][0] + fbn[pixid][gain][0]*fbn[pixid][gain][0];
    return sqrt(sum);
  }
};

//=====================================================================
// a class to calibrate timing of individual capacitors using Taka's method,
// i.e. use a very long calibration run and check how often the peak falls
// in a given capacitor
class TimeTakaClu
{
private:
  int fNumCap; // number of capacitors (1024 or 4096)
public:
  int nhit[nch][4096]; // number of hits in a given capacitor
  float width[nch][4096]; // sample width (in "average" samples)
  float dwidth[nch][4096]; // error on sample width (in "average" samples)
  TimeTakaClu(int nc)
  {
    fNumCap=nc;
    for (int ich=0; ich<nch; ich++)
      for (int i=0; i<fNumCap; i++)
	{
	  nhit[ich][i]=0;
	  width[ich][i]=1;
	  dwidth[ich][i]=-1;
	}
  }

  bool FillFile(string namecal, PedestalSimple &pedsim, int roi)
  {
    ifstream plik2(namecal.data(), ios::binary);
    if (!plik2.is_open())
    {
      cout<<"file: "<<namecal<<" is not open, exiting"<<endl;
      return false;
    }

    Event evt2(roi);
    int ev=0; 
    while (evt2.read(plik2))
      {
	if (ev % 10000 == 0)
	  cout<<"Event "<<ev<<endl;
	ev++;
	
	removePed(evt2, pedsim);
	evt2.CorrTime();
	// evt2.InterpolatePseudoPulses();
	int skip=2;
	for (int ich=0; ich<nch-1; ich++)
	  {
	    int bestslice=-1;
	    float bestsig=-1.e7;
	    for (int i=skip; i<roi-skip; i++)
	      if (evt2.samples[ich][0][i]>bestsig)
		{
		  bestsig=evt2.samples[ich][0][i];
		  bestslice=i;
		}
	    if (bestsig>0)
	      {
		int abspos= (bestslice + evt2.firstcap[ich][0])%fNumCap;
		nhit[ich][abspos]++;
	      }
	  }
      }  
    cout<<ev<<" events read"<<endl;
    plik2.close();
    return true;
  }

  bool Finalize()
  {
    // first check if there are any fields without hits:
    bool isgood=true;
    int nhitall[nch]={0};
    for (int ich=0; ich<nch-1; ich++)
      for (int i=0; i<fNumCap; i++)
	{
	  nhitall[ich]+=nhit[ich][i];
	  if (nhit[ich][i]<1)
	    {
	      cout<<"chanel "<<ich<<" capacitor "<<i<<"nhit="<<nhit[ich][i]<<endl;
	      isgood=false;
	    }
	}
    if (!isgood)
      {
	cout<<"found capacitors without sufficient number of hits, exiting";
	 return false;
      }
    
    for (int ich=0; ich<nch-1; ich++)
      for (int i=0; i<fNumCap; i++)
	{
	  width[ich][i]=1.*nhit[ich][i]*fNumCap/nhitall[ich];    
	  dwidth[ich][i]=sqrt(1.*nhit[ich][i])*fNumCap/nhitall[ich];    
	}
   }
};
#endif
