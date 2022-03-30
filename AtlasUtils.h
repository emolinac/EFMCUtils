//
//   @file    AtlasUtils.h         
//   
//
//   @author M.Sutton
// 
//   Copyright (C) 2010 Atlas Collaboration
//
//   $Id: AtlasUtils.h, v0.0   Thu 25 Mar 2010 10:34:20 CET $


#ifndef __ATLASUTILS_H
#define __ATLASUTILS_H

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"

void ATLAS_LABEL(Double_t x,Double_t y,Color_t color=1); 

TGraphErrors* myTGraphErrorsDivide(TGraphErrors* g1,TGraphErrors* g2);

TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);

TGraphAsymmErrors* myMakeBand(TGraphErrors* g0, TGraphErrors* g1,TGraphErrors* g2);

void myAddtoBand(TGraphErrors* g1, TGraphAsymmErrors* g2);

TGraphErrors* TH1TOTGraph(TH1 *h1);

void myText(Double_t x,Double_t y,Color_t color,const char *text,Double_t tsize=22,Double_t langle=0);

void myBoxText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor,const char *text, Int_t mlinecolor=1, Int_t mlinestyle=1, Int_t fstyle=1001, Double_t tsize=22);

void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle,const char *text,Float_t msize=2.,Double_t tsize=22);

void myLineText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor,const char *text,Int_t lstyle, const Double_t tsize=22);

TCanvas* MakeCanvas(Double_t wx=480, Double_t wy=480, const char* name="c", const char* title="title") {return new TCanvas(name,title,wx,wy);}

#endif // __ATLASUTILS_H
