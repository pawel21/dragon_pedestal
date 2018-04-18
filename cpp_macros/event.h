#ifndef event_H
#define event_H

const int size4drs=4*1024;
const int nch=8;

//#define USEINTEGER // for test of checking how bad the corrections works if only integers can be used
inline float pedtimecorr(float timediff)
{
  // return 0; // switch off the correction !!!
  // if (timediff>0)
  if ((timediff>0) && (timediff<45))
    {
      // return  28.98 * pow(timediff,-0.2405) -11.95;
      // return  23. * pow(timediff,-0.2405) -8;
#ifdef USEINTEGER
      return  round(29.3 * pow(timediff,-0.2262) -12.4);
#else
      return  29.3 * pow(timediff,-0.2262) -12.4;
#endif
    }
  else return 0;
}

class Event
{
public:
  unsigned short roisize;
  unsigned short cntpps;
  unsigned int cnt10mhz;
  unsigned int evtcnt;
  unsigned int trigcnt;
  unsigned long cnt133mhz;
  unsigned short flags[nch];
  unsigned short firstcap[nch][2];
  unsigned short oldfirstcap[nch][2];
  unsigned short oldoldfirstcap[nch][2]; // first cap in 2 events ago
  unsigned short oldoldoldfirstcap[nch][2]; // first cap in 2 events ago
  // unsigned short *samples[nch][2];
  float *samples[nch][2];
  int samplestag[nch][2][size4drs];

  double lasttime[nch][2][size4drs];
  bool read(ifstream &plik);
  bool convertendian;
  Event(short r);
  ~Event();
  void Print(bool alsosamples=false);
  double GetDRSTime(){return cnt133mhz*(1./133.e6)*1.e3;}// [ms] from 133MHz clock
  double GetDRSTime2(){return cnt10mhz*(1./10.e6)*1.e3;}// [ms] from 10 MHz clock (much less bits)
  void ResetTime();
  void CorrTime();
  void InterpolatePseudoPulse(int pix, int gain, int abspos);
  void InterpolatePseudoPulses();
  float SumSlices(int pix, int gain, int pos0, int num);
  float Sliding(int pix, int gain, int num=6, float *time=0, int nfirst=-1, int nlast=-1, bool usehalfs=false);
};

Event::Event(short r)
{
  roisize=r;
  convertendian=true;
  // convertendian=false;
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      // samples[i][j] = new unsigned short [roisize];
      samples[i][j] = new float [roisize];
  ResetTime();
}

Event::~Event()
{
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      delete []samples[i][j];
}

bool Event::read(ifstream &plik)
{
  //cout<<"Hello here event"<<endl;
  char header[64];
  plik.read(header, 64);
  //cout<<"header: "<< *(unsigned short*)header <<endl;
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      {
	oldoldoldfirstcap[i][j]=oldoldfirstcap[i][j];
	oldoldfirstcap[i][j]=oldfirstcap[i][j];
	oldfirstcap[i][j]=firstcap[i][j];
      }
  // big endian to small endian conversion
  if (convertendian)
    {
      swap (header[2], header[3]);
      swap (header[4], header[7]);
      swap (header[5], header[6]);
      swap (header[8], header[11]);
      swap (header[9], header[10]);
      swap (header[12], header[15]);
      swap (header[13], header[14]);
      swap (header[16], header[23]);
      swap (header[17], header[22]);
      swap (header[18], header[21]);
      swap (header[19], header[20]);
      for (int i=0; i<2*nch; i++)
	swap (header[32+i*2], header[32+i*2+1]);
    }

  unsigned short dms=*(unsigned short*)header;
  //cout<<"dms = "<<dms<<endl;
  if (dms!=0xAAAA)
    {
      cout<<"header error, found "<<dms<<" instead of "<<0xAAAA<<endl;
      return false;
    }
  cntpps=*(unsigned short*)(header+2);
  //cout<<"cntpps: "<<cntpps<<endl;
  cnt10mhz=*(unsigned int*)(header+4);
  //cout<<"cnt10mhz: "<<cnt10mhz<<endl;
  evtcnt=*(unsigned int*)(header+8);
  //cout<<"evtcnt: "<<evtcnt<<endl;
  trigcnt=*(unsigned int*)(header+12);
  //cout<<"trigcnt: "<<trigcnt<<endl;
  cnt133mhz=*(unsigned long*)(header+16);
  //cout<<"cnt133mhz: "<<cnt133mhz<<endl;
  unsigned long dml = *(unsigned long*)(header+24);
  //cout<<"dml: "<<dml<<endl;

  if (dml!=0xDDDDDDDDDDDDDDDD)
    {
      cout<<"header error, found "<<dml<<" instead of "<<0xDDDDDDDDDDDDDDDD<<endl;
      return false;
    }
  for (int i=0; i<nch; i++) {
    flags[i] = *(unsigned short*)(header+32+i*2);
    //cout<<"flags "<<i<<" = "<<flags[i]<<endl;
   }
  for (int i=0; i<nch; i++)
    {
      int ich=(i/2)*2;
      int gain=i%2;
      firstcap[ich][gain]   = *(unsigned short*)(header+48+i*2); // same chip
      firstcap[ich+1][gain] = *(unsigned short*)(header+48+i*2); // is used for both
      cout<<"first cap = " << firstcap[ich][gain] << endl;
    }

  // now read the actual event
  for (int k=0; k<roisize*2; k++)
    {
      char data[16];
      plik.read(data, 16);
      if (convertendian)
	for (int i=0; i<8; i++)
	  swap (data[i*2], data[i*2+1]);
      for (int i=0; i<4; i++)
	{
	  int pix=(k<roisize)?i*2:i*2+1;
	  //cout<<"pix = " << pix << endl;
	  samples[pix][0][k%roisize]=*(unsigned short*)(data+i*4); // high gain
      //cout<<"samples high gain: "<<*(unsigned short*)(data+i*4)<<endl;
	  samples[pix][1][k%roisize]=*(unsigned short*)(data+i*4+2); // low gain
      //cout<<"samples low gain: "<<*(unsigned short*)(data+i*4+2)<<endl;
	}
    }
  return true;
}

void Event::Print(bool alsosamples)
{
  if (alsosamples) cout<<" ========= "<<endl;
  cout<<"PPS="<<cntpps<<", cnt10MHz="<<cnt10mhz<<", evt="<<evtcnt<<", trig="<<trigcnt<<", cnt133MHz="<<cnt133mhz<<endl;
  cout<<"flags: ";
  for (int i=0; i<nch; i++)
    cout<<flags[i]<<" ";
  cout<<endl<<"firstcap[high]=";
  for (int i=0; i<nch; i++)
    cout<<firstcap[i][0]<<" ";
  cout<<endl<<"firstcap[low]=";
  for (int i=0; i<nch; i++)
    cout<<firstcap[i][1]<<" ";
  cout<<endl;
  if (alsosamples)
    for (int i=0; i<nch; i++)
      for (int j=0; j<2; j++)
	{
	  cout<<"Channel "<<i<<((j==0)?" high":" low")<<" gain: ";
	  for (int k=0; k<roisize; k++)
	    cout<<samples[i][j][k]<<" ";
	  cout<<endl;
	}
}

void Event::ResetTime()
{
  // cout<<"resetting time"<<endl;
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      {
	for (int k=0; k<size4drs; k++)
	  lasttime[i][j][k]=-1;
	firstcap[i][j]=0;
      }
}

void Event::CorrTime()
{
  double timenow=GetDRSTime();
  cout<<"Corr Time !!!"<<endl;
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
        {
	    int fc=firstcap[i][j];
	    for (int k=0; k<roisize; k++)
	        {
	        int posabs = (k+fc)%size4drs;
		   //cout<<"posabs = "<<posabs<<endl;
	  // if (i==0 && j==0)
	  // cout<<lasttime[i][j][posabs]<<endl;
	        if (lasttime[i][j][posabs]>0)
	        {
	        double timediff=timenow - lasttime[i][j][posabs]; //ms
	        double timecorr=pedtimecorr(timediff);
		   //cout<<"timecorr = "<<timecorr<<endl;
	      // if (posabs==(oldfirstcap[i][j]+roisize-1))
	      // 	timecorr/=2; // if the slice was the last slice of previous ROI the correction is smaller
	        samples[i][j][k]-=timecorr;
	      // cout<<samples[i][j][k]<<" "<<pedtimecorr(timediff)<<" "<<timediff<<endl;
	        }
	        if (k<roisize-1) // added by JS
	        lasttime[i][j][posabs]=timenow;      // update the times
	        samplestag[i][j][posabs]=0;
	        }
	        samplestag[i][j][(roisize-1+fc)%size4drs]=1;

	// now the magic of Dragon, if the ROI is in the last quarter of each DRS4
	// for even channel numbers extra 12 slices are read in a different place
	// code from Takayuki
	        if (i%2 == 0)
	        {
	            if (fc%1024>766 && fc%1024 < 1012)
	                for(int kk=fc+1024-1;kk<fc+1024+11;kk++) // changed from 12 to 11 here after seeing Problem #1035 !
		                {
		                lasttime[i][j][kk%size4drs] = timenow;
		  // samplestag[i][j][kk%4096]=100+fc;
		                  }

	            else if(fc%1024 >= 1012) // it was '>'
	            {
		            int channel = fc/1024;
		    for(int kk=fc+1024;kk<(channel+2)*1024;kk++) // corrected a bug from +1 to +2
		        {
		        lasttime[i][j][kk%4096] = timenow;
		    // samplestag[i][j][kk%4096]=1;
		        }
	            }
	        }

	  // figured out that this correction is not needed and actually
	  // makes things worse
	  // if(fc%1024 >=1024-roisize) // before was in if (i%2 == 0 loop)
	  //   {
	  //     int channel = fc/1024;
	  //     lasttime[i][j][1024*channel] = timenow;
	  //     // samplestag[i][j][1024*channel]=1;
	  //   }
      }
}


void Event::InterpolatePseudoPulse(int pix, int gain, int abspos) // abspos and abspos+1 - absolute position in the DRS4
{
  // first check where it is in ROI
  int pos = (abspos - firstcap[pix][gain] + size4drs) % size4drs;
  if ( ((pos>=0)&&(pos<=roisize-1)) || (pos==size4drs - 1) )
    {
      if (pos==0) // XX...
	{
	  samples[pix][gain][0]=samples[pix][gain][2];
	  samples[pix][gain][1]=samples[pix][gain][2];
	}
      else if (pos==size4drs - 1) // X....
	{
	  samples[pix][gain][0]=samples[pix][gain][1];
	}
      else if (pos==roisize-2) // ...XX
	{
	  samples[pix][gain][roisize-2]=samples[pix][gain][roisize-3];
	  samples[pix][gain][roisize-1]=samples[pix][gain][roisize-3];
	}
      else if (pos==roisize-1) // ....X
	{
	  samples[pix][gain][roisize-1]=samples[pix][gain][roisize-2];
	}
      else //...XX...
	{
	  samples[pix][gain][pos]=samples[pix][gain][pos-1] + (samples[pix][gain][pos+2]-samples[pix][gain][pos-1])*0.33;
	  samples[pix][gain][pos+1]=samples[pix][gain][pos-1] + (samples[pix][gain][pos+2]-samples[pix][gain][pos-1])*0.66;

#ifdef USEINTEGER
	  samples[pix][gain][pos]=round(samples[pix][gain][pos]);
	  samples[pix][gain][pos+1]=round(samples[pix][gain][pos+1]);
#endif
	}
      // cout<<"interpolated at pos = "<<pos<<", pix="<<pix<<", gain="<<gain<<endl;
    }

}

void Event::InterpolatePseudoPulses()
{
  // interpolation of spike type A
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      for (int k=0; k<4; k++) // for each DRS4 in the chain
	{
	  InterpolatePseudoPulse(i, j, (roisize-2 + oldfirstcap[i][j]+k*1024));
	  InterpolatePseudoPulse(i, j, (1024 - roisize-2 - oldfirstcap[i][j]+k*1024+size4drs));
	}

  // interpolation of spike type B
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      {
  	int spikepos = ((short)oldfirstcap[i][j]-1 - (short)firstcap[i][j]+2*size4drs)%size4drs;

  	if (spikepos==0) // X...
  	  samples[i][j][0]=samples[i][j][1];
  	else if (spikepos==roisize-1)  // ...X
  	  samples[i][j][roisize-1]=samples[i][j][roisize-2];
  	else if (spikepos<roisize-1)
	  {
	    samples[i][j][spikepos]=0.5*(samples[i][j][spikepos-1]+samples[i][j][spikepos+1]);
#ifdef USEINTEGER
	    samples[i][j][spikepos]=round(samples[i][j][spikepos]);
#endif
	  }
      }

  // interpolation of every 32nd slice
  // just crude interpolation
  for (int i=0; i<nch; i++)
    for (int j=0; j<2; j++)
      {
	if (firstcap[i][j] % 32 == 31) // X...
	  samples[i][j][0]=samples[i][j][1];
	if ((firstcap[i][j]+roisize-1)% 32 == 31) // ...X
	  samples[i][j][roisize-1]=samples[i][j][roisize-2];
	// FIX ME! instead of for loop do a while loop with +32 increment
	for (int k=1; k<roisize-1; k++) // ..X..
	  if ((firstcap[i][j]+k) % 32 == 31)
	    {
	      samples[i][j][k]=0.5*(samples[i][j][k-1]+samples[i][j][k+1]);
#ifdef USEINTEGER
	      samples[i][j][k]=round(samples[i][j][k]);
#endif
	    }
      }

}


int maxpos(unsigned short *tab, int roi)
{
  short max=-1000;
  int maxi=-1;
  for (int i=0; i<roi; i++)
    {
      if (tab[i]> max)
	{
	  max=tab[i];
	  maxi=i;
	}
    }
  return maxi;
}

float Event::SumSlices(int pix, int gain, int pos0, int num)
{
  float sum=0;
  for (int i=pos0; i<pos0+num; i++)
    sum+=samples[pix][gain][i];
  return sum;
}

// sliding window of size num, skipping 2 first and 2 last slices
float Event::Sliding(int pix, int gain, int num, float *time, int nfirst, int nlast, bool usehalfs)
{
  const int skip=2;
  if (nfirst<0)
    nfirst=skip;
  if (nlast<0)
    nlast=roisize-1-skip;
  float sum=0;
  float *sam = samples[pix][gain];
  for (int i=nfirst; i<num+nfirst; i++)
    sum+=sam[i]; //first sum
  float bestsum=sum;
  int bestpos=nfirst;
  // for (int i=skip; i<roisize-num-skip; i++)
  for (int i=nfirst; i<nlast+1-num; i++)
    {
      sum+=(sam[i+num]-sam[i]);
      if (sum>bestsum)
	{
	  bestsum=sum;
	  bestpos=i+1; // first slice of integration
	}
    }
  float peak=-1111;
  if (time)
    {
      *time = 0;
      for (int i=bestpos; i<bestpos+num; i++)
	{
	  if (sam[i]>peak)
	    peak=sam[i];
	  *time+=sam[i]*i;
	}

      // simple tweak to reduce the problem of the partial discretization of the arrival times
      // after finding the biggest sum of slices we check if we can shift the window by half a
      // slice to get a better sum
      // does not work either way :-(
      // if (usehalfs)
      // 	{
      // 	  int bestlast = bestpos+num - 1; // last integrated slice
      // 	  float nexttoleft=-1, nexttoright=-1;
      // 	  if (bestpos>nfirst)
      // 	    nexttoleft=0.25*(sam[bestpos]+sam[bestpos-1])-0.5*sam[bestlast];
      // 	  if (bestlast<nlast)
      // 	    nexttoright=0.25*(sam[bestlast]+sam[bestlast+1])-0.5*sam[bestpos];
      // 	  bool corrtoleft=false, corrtoright=false;
      // 	  if (nexttoleft>0 && nexttoleft>nexttoright)
      // 	    corrtoleft=true;
      // 	  else if (nexttoright>0)
      // 	    corrtoright=true;

      // 	  if (corrtoleft)
      // 	    {
      // 	      *time+=(bestpos-0.5)*0.25*(sam[bestpos]+sam[bestpos-1])-0.5*bestlast*sam[bestlast];
      // 	      bestsum+=corrtoleft;
      // 	    }
      // 	  else if(corrtoright)
      // 	    {
      // 	      *time+=(bestlast+0.5)*0.25*(sam[bestlast]+sam[bestlast+1])-0.5*bestpos*sam[bestpos];
      // 	      bestsum+=corrtoright;
      // 	    }

      // 	}
      *time/=bestsum;
      // int bestlast = bestpos+num - 1; // last integrated slice
      // *time=(*time - 0.5*(bestpos*sam[bestpos] + bestlast*sam[bestlast])) /(bestsum - 0.5*(sam[bestpos] + sam[bestlast])) ;
      // cout<<bestpos<<" "<<bestsum<<" "<<sow<<" "<<*time<<endl;
    }
  // return peak; // for tests, only returns the peak
  return bestsum;
}



// merge events from different files and read next such event
// pretty buggy, but at least it should return an error if missing some event
// evts - array of pointers to events
// returns false if one of the files finishes
bool constructevent(Event **evts, ifstream **files, int nfiles)
{
  // read the first event
  if (!evts[0]->read(*files[0]))
    return false;
  unsigned int trigcnt0=evts[0]->trigcnt;
  // cout<<"read in file 0 with trigcnt="<<trigcnt0<<endl;
  for (int i=1; i<nfiles; i++)
    {
      do
	{
	  if (!evts[i]->read(*files[i]))
	    return false;

	  // cout<<"read in file "<<i<<" with trigcnt="<<evts[i]->trigcnt<<endl;

	}
      while (evts[i]->trigcnt!=trigcnt0);
    }
  return true;
}


#endif
