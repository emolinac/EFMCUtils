#ifndef __EFMCUTILS_H
#define __EFMCUTILS_H

#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TNtuple.h"
#include "TCut.h"
#include "TMath.h"
/*______________________________General_________________________________*/

Bool_t EmptyHisto(TH1F* h);

void GetQ2NuCentroidLiquid(Int_t Q2_bin, Int_t Nu_bin, TNtuple* limits_tuple, TNtuple* ntuple_data, TH2F* h, Double_t* pairQ2Nu);

void GetQ2NuCentroidSolid(Int_t Q2_bin, Int_t Nu_bin, TNtuple* limits_tuple, TNtuple* ntuple_data, TH2F* h, Double_t* pairQ2Nu);

/*____________________________kt2 Analysis_______________________________*/
Double_t GetMeanPt2Signori(Double_t* x, Double_t* pars);

Double_t GetMeanPt2Beta(Double_t* x, Double_t* pars);

/*____________________________Pt2 Analysis_______________________________*/
Double_t GetMaxPt2(TH1F* h);

Int_t GetMaxPt2Bin(TH1F* h);

/*___________________________TGraph Customs___________________________*/
void SetXErrNull(TGraphErrors* g);

void SetXShift(TGraphErrors* g, Double_t shift);

/*____________________________Calculations_____________________________*/
Double_t BiggestContentAvge(TGraphErrors* g1 , TGraphErrors* g2);
  
#endif
