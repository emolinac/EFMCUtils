#ifndef __EFMCUTILS_H
#define __EFMCUTILS_H

#include "TH1F.h"
#include "TGraphErrors.h"

/*______________________________General_________________________________*/

Bool_t EmptyHisto(TH1F* h);

/*____________________________Pt2 Analysis_______________________________*/
Double_t GetMaxPt2(TH1F* h);

Int_t GetMaxPt2Bin(TH1F* h);

/*___________________________TGraph Customs___________________________*/
void SetXErrNull(TGraphErrors* g);

void SetXShift(TGraphErrors* g, Double_t shift);

/*____________________________Calculations_____________________________*/
Double_t BiggestContentAvge(TGraphErrors* g1 , TGraphErrors* g2);
  
#endif
