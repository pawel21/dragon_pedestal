#include <stdio.h>

int translatezFits(unsigned short * eventWfSrcHG,
		   unsigned short * eventWfSrcLG,
		   unsigned short * eventDRSSrcHG,
		   unsigned short * eventDRSSrcLG,
		   unsigned short * eventWfDestHG,
		   unsigned short * eventWfDestLG,
		   unsigned nbModules)
{
  int i, j, k = 0;
  for(k=0; k<nbModules; k++)
    {
      printf("k= %i \n", k);
      unsigned short * wfSrcHG = &(eventWfSrcHG[k*7*40]);
      unsigned short * wfSrcLG = &(eventWfSrcLG[k*7*40]);
      
      unsigned short * wfDestHG = &(eventWfDestHG[k*8*40]);
      unsigned short * wfDestLG = &(eventWfDestLG[k*8*40]);

      unsigned short * DRSSrcHG  = &(eventDRSSrcHG[k*40]);
      unsigned short * DRSSrcLG  = &(eventDRSSrcLG[k*40]);

      unsigned short * currentSrc;
      int line = 0;
      
      for(j=0; j<8; j++)
	{
	  if(j%2==0)
	    {
	      currentSrc = &(wfSrcHG[line*40]);
	    }
	  else
	    {
	      currentSrc = &(wfSrcLG[line*40]);
	      line++;
	    }
	  
	  for(i=0; i<5; i++)
	    {
	      wfDestHG[0+i+j*5] = currentSrc[i*8+0]; /* Pix 0 HG */
	      wfDestLG[0+i+j*5] = currentSrc[i*8+1]; /* Pix 0 LG */
	      wfDestHG[80+i+j*5] = currentSrc[i*8+2]; /* Pix 2 HG */
	      wfDestLG[80+i+j*5] = currentSrc[i*8+3]; /* Pix 2 LG */
	      wfDestHG[160+i+j*5] = currentSrc[i*8+4]; /* Pix 4 HG */
	      wfDestLG[160+i+j*5] = currentSrc[i*8+5]; /* Pix 4 LG */
	      wfDestHG[240+i+j*5] = currentSrc[i*8+6]; /* Pix 6 HG */
	      wfDestLG[240+i+j*5] = currentSrc[i*8+7]; /* Pix 6 LG */
	    }
	}
      line = 4;
      for(j=0; j<6; j++)
      	{
	  if(j%2==0)
	    {
	      currentSrc = &(wfSrcHG[line*40]);
	    }
	  else
	    {
	      currentSrc = &(wfSrcLG[line*40]);
	      line++;
	    }
      	  for(i=0; i<5; i++)
      	    {
      	      wfDestHG[40+i+j*5] = currentSrc[i*8+0]; /* Pix 1 HG */
      	      wfDestLG[40+i+j*5] = currentSrc[i*8+1]; /* Pix 1 LG */
      	      wfDestHG[120+i+j*5] = currentSrc[i*8+2]; /* Pix 3 HG */
      	      wfDestLG[120+i+j*5] = currentSrc[i*8+3]; /* Pix 3 LG */
      	      wfDestHG[200+i+j*5] = currentSrc[i*8+4]; /* Pix 5 HG */
      	      wfDestLG[200+i+j*5] = currentSrc[i*8+5]; /* Pix 5 LG */
      	      wfDestHG[280+i+j*5] = currentSrc[i*8+6]; /* Pix DRS HG */
      	      wfDestLG[280+i+j*5] = currentSrc[i*8+7]; /* Pix DRS LG */
      	    }
      	}
      
      for(j=0; j<2; j++)
      	{
      	  currentSrc = j%2==0 ? DRSSrcHG:DRSSrcLG;
      	  for(i=0; i<5; i++)
      	    {
      	      wfDestHG[40+i+j*5+30] = currentSrc[i*8+0]; /* Pix 1 HG */
      	      wfDestLG[40+i+j*5+30] = currentSrc[i*8+1]; /* Pix 1 LG */
      	      wfDestHG[120+i+j*5+30] = currentSrc[i*8+2]; /* Pix 3 HG */
      	      wfDestLG[120+i+j*5+30] = currentSrc[i*8+3]; /* Pix 3 LG */
      	      wfDestHG[200+i+j*5+30] = currentSrc[i*8+4]; /* Pix 5 HG */
      	      wfDestLG[200+i+j*5+30] = currentSrc[i*8+5]; /* Pix 5 LG */
      	      wfDestHG[280+i+j*5+30] = currentSrc[i*8+6]; /* Pix DRS HG */
      	      wfDestLG[280+i+j*5+30] = currentSrc[i*8+7]; /* Pix DRS LG */
      	    }
      	}
    }
  
  return 0;
}
