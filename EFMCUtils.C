#include "EFMCUtils.h"

#include <iostream>

/*____________________________General_________________________________*/

Bool_t EmptyHisto(TH1F* h){
  // Return kTRUE if histo is empty. Otherwise, returns kFALSE.
  for(Int_t bin = 1; bin <= h->GetNbinsX() ; bin++){
    if(h->GetBinContent(bin)!=0){return kFALSE;}
  }
  return kTRUE;
}

void GetQ2NuCentroidLiquid(Int_t Q2_bin, Int_t Nu_bin, TNtuple* limits_tuple, TNtuple* ntuple_data, TH2F* h, Double_t* pairQ2Nu){
  // Returns (Q2_centroid, Nu_centroid) for the liquid target
  
  Float_t Q2_min=0,Q2_max=0,Nu_min=0,Nu_max=0;
  Float_t Q2_min_local=0, Q2_max_local=0, Nu_min_local=0, Nu_max_local=0, Q2_bin_local=0, Nu_bin_local=0;
  
  limits_tuple->SetBranchAddress("Q2_min",&Q2_min_local);
  limits_tuple->SetBranchAddress("Q2_max",&Q2_max_local);
  limits_tuple->SetBranchAddress("Nu_min",&Nu_min_local);
  limits_tuple->SetBranchAddress("Nu_max",&Nu_max_local);
  limits_tuple->SetBranchAddress("Q2_bin",&Q2_bin_local);
  limits_tuple->SetBranchAddress("Nu_bin",&Nu_bin_local);

  for(Int_t tuple_entry = 0 ; tuple_entry < limits_tuple->GetEntries() ; tuple_entry++)
    {
      limits_tuple->GetEntry(tuple_entry);
      if(Q2_bin==Q2_bin_local&&Nu_bin==Nu_bin_local)
	{
	  Q2_min = Q2_min_local;
	  Q2_max = Q2_max_local;
	  Nu_min = Nu_min_local;
	  Nu_max = Nu_max_local;
	  break;
	}
    }

  TCut cut = Form("VC_TM==1&&%f<Q2&&Q2<%f&&%f<Nu&&Nu<%f",Q2_min,Q2_max,Nu_min,Nu_max);
      
  ntuple_data->Project("h","Nu:Q2",cut);

  pairQ2Nu[0] = h->GetMean(1);
  pairQ2Nu[1] = h->GetMean(2);
  std::cout<<"Q2_centroid = "<<pairQ2Nu[0]<<"   Nu_centroid = "<<pairQ2Nu[1]<<std::endl;

  h->Reset();
}

void GetQ2NuCentroidSolid(Int_t Q2_bin, Int_t Nu_bin, TNtuple* limits_tuple, TNtuple* ntuple_data, TH2F* h, Double_t* pairQ2Nu){
  // Returns (Q2_centroid, Nu_centroid) for the solid target
  
  Float_t Q2_min=0,Q2_max=0,Nu_min=0,Nu_max=0;
  Float_t Q2_min_local=0, Q2_max_local=0, Nu_min_local=0, Nu_max_local=0, Q2_bin_local=0, Nu_bin_local=0;
  
  limits_tuple->SetBranchAddress("Q2_min",&Q2_min_local);
  limits_tuple->SetBranchAddress("Q2_max",&Q2_max_local);
  limits_tuple->SetBranchAddress("Nu_min",&Nu_min_local);
  limits_tuple->SetBranchAddress("Nu_max",&Nu_max_local);
  limits_tuple->SetBranchAddress("Q2_bin",&Q2_bin_local);
  limits_tuple->SetBranchAddress("Nu_bin",&Nu_bin_local);

  for(Int_t tuple_entry = 0 ; tuple_entry < limits_tuple->GetEntries() ; tuple_entry++)
    {
      limits_tuple->GetEntry(tuple_entry);
      if(Q2_bin==Q2_bin_local&&Nu_bin==Nu_bin_local)
	{
	  Q2_min = Q2_min_local;
	  Q2_max = Q2_max_local;
	  Nu_min = Nu_min_local;
	  Nu_max = Nu_max_local;
	  break;
	}
    }

  TCut cut = Form("VC_TM==2&&%f<Q2&&Q2<%f&&%f<Nu&&Nu<%f",Q2_min,Q2_max,Nu_min,Nu_max);
      
  ntuple_data->Project("h","Nu:Q2",cut);

  pairQ2Nu[0] = h->GetMean(1);
  pairQ2Nu[1] = h->GetMean(2);
  std::cout<<"Q2_centroid = "<<pairQ2Nu[0]<<"   Nu_centroid = "<<pairQ2Nu[1]<<std::endl;

  h->Reset();
}

/*_________________________kt2 Analysis_______________________________*/
Double_t GetMeanPt2Signori(Double_t* x, Double_t* pars){
  Double_t Zh = x[0];
  // meanpt2 parameters
  Double_t zhat  = pars[0];
  Double_t norm  = pars[1];
  Double_t beta  = pars[2];
  Double_t delta = pars[3];
  Double_t gamma = pars[4];

  Double_t meanpt2 = norm*(TMath::Power(Zh,beta)+delta)*TMath::Power(1-Zh,gamma)/(TMath::Power(zhat,beta)+delta)*TMath::Power(1-zhat,gamma);
  
  // meankt2 parameter
  Double_t meankt2 = pars[5];

  return meanpt2 + Zh*Zh*meankt2;
}

Double_t GetMeanPt2Beta(Double_t* x, Double_t* pars){
  Double_t Zh = x[0];
  // meanpt2 parameters
  Double_t norm   = pars[0];
  Double_t alpha  = pars[1];
  Double_t beta   = pars[2];

  Double_t meanpt2 = norm*TMath::BetaDist(Zh,alpha,beta);
  
  // meankt2 parameter
  Double_t meankt2 = pars[3];

  return meanpt2 + Zh*Zh*meankt2;
}

/*_________________________Pt2 Analysis_______________________________*/

Double_t GetMaxPt2(TH1F* h){
  // Returns the MAX value of Pt2 given a certain Pt2 distribution.
  // IMPORTANT: This functions assumes that tribution does not have empty bins
  
  if(h==NULL){std::cout<<"Histo is null!"<<std::endl; return 0;}
  for(Int_t bin = 1 ; bin <= h->GetNbinsX() ; bin++){
    Double_t value = h->GetBinContent(bin);
    if(value==0){
      return h->GetBinLowEdge(bin);
    }
  }
  return 3.;
}

Int_t GetMaxPt2Bin(TH1F* h){
  // Returns the bin where the MAX value of a Pt2 distribution is located.
  // IMPORTANT: This functions assumes that tribution does not have empty bins

  if(h==NULL){std::cout<<"Histo is null!"<<std::endl; return 0;}  
  for(Int_t bin = 1 ; bin <= h->GetNbinsX() ; bin++){
    Double_t value = h->GetBinContent(bin);
    if(value==0){
      return bin-1;
    }
  }
  return h->GetNbinsX();
}

/*___________________________TGraph Customs___________________________*/

void SetXErrNull(TGraphErrors* g){
  // Function that removes the errors in the X axis for a TGraph
  Double_t* errors_Y = g->GetEY();
  for(Int_t point = 0 ; point < g->GetN() ; point++)
    {
      g->SetPointError(point,0,errors_Y[point]);  
    }
}

void SetXShift(TGraphErrors* g, Double_t shift){
  // Function that sets a uniform in the X axis for all points in a TGraph
  Double_t* content_Y = g->GetY();
  Double_t* content_X = g->GetX();
  for(Int_t point = 0; point < g->GetN() ; point++)
    {
      g->SetPoint(point, content_X[point] + shift , content_Y[point]);
    }
}

/*____________________________Calculations_____________________________*/
Double_t BiggestContentAvge(TGraphErrors* g1 , TGraphErrors* g2){
  // Return the average of a sample containing the biggest contents of g1 and g2
  // This method is useful to calculate the average systematic error between two variations
  // of the Nominal
  
  Int_t g1_N = g1->GetN();
  Int_t g2_N = g2->GetN();

  if(g1_N!=g2_N){std::cout<<"Number of points in TGraphs is not equal!"<<std::endl; return 0;}

  Double_t* g1_content = g1->GetY();
  Double_t* g2_content = g2->GetY();

  Double_t sum = 0; Int_t N = 0;
  for(Int_t point = 0 ; point < g1_N ; point++){
    if(TMath::Abs(g1_content[point])>TMath::Abs(g2_content[point])){
      sum += TMath::Abs(g1_content[point]);
      N++;
    }
    else{
      sum += TMath::Abs(g2_content[point]);
      N++;
    }
  }
  if(N==0){std::cout<<"N==0!. Invalid operation."<<std::endl; return 0;}
  return sum/N;
}
