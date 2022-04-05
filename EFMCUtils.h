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
Double_t GetMeanpt2Signori(Double_t Zh, Double_t zhat, Double_t norm, Double_t beta, Double_t delta, Double_t gamma);

Double_t GetMeanpt2Beta(Double_t Zh, Double_t norm, Double_t alpha, Double_t beta);

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

Double_t BiggestContentAvge_CustomStartPoint(TGraphErrors* g1 , TGraphErrors* g2, Int_t start_point);

#endif
