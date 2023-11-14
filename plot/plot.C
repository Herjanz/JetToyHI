//Quick dirty plotting macro for jet response from trees of JetToyHI framework

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH3F.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

using namespace std;

int markerColor[3] = {1,kRed+1,4};
int markerStyle[3] = {20,25,24};

TH1F *DrawFrame(double xmin = 0., double xmax = 1., double ymin = 0., double ymax = 1., TString xTitle = "x", TString yTitle = "y", bool setMargins = true);

TTree *getTree(TString str = "JetToyHIResultSimpleJetAnalysis.root", TString tree = "HerjansBoompje");

void phi_eta(TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();


  auto C = new TCanvas();

  // partons
  tr->Draw("partonsPhi:partonsEta>>partons");
  TH1F *hDiv = (TH1F*)gDirectory->Get("partons");
  hDiv->SetMarkerStyle(kFullCircle);
  hDiv->SetMarkerColor(1);
  // jets
  tr->Draw("sigJetPhi:sigJetEta>>jets");
  TH1F *hDiv2 = (TH1F*)gDirectory->Get("jets");
  hDiv2->SetMarkerStyle(kFullCircle);
  hDiv2->SetMarkerColor(2);
  // daughterpartons
  tr->Draw("daughterPartonsPhi:daughterPartonsEta>>daughterPartons");
  TH1F *hDiv3 = (TH1F*)gDirectory->Get("daughterPartons");
  hDiv3->SetMarkerStyle(kFullCircle);
  hDiv3->SetMarkerColor(3);


  hDiv->GetXaxis()->SetTitle("#eta");
  hDiv->GetYaxis()->SetTitle("#phi");


  TH1 *h1 = hDiv->DrawCopy();
  TH1 *h2 = hDiv2->DrawCopy("same");
  TH1 *h3 = hDiv3->DrawCopy("same");
  
  TLegend *legend = new TLegend(0.55,0.65,0.76,0.82);
  legend->AddEntry(h1,"partons");
  legend->AddEntry(h3,"daughterpartons");
  legend->AddEntry(h2,"jets");
  legend->Draw();
}

void dp2p_pt_r(TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();


  auto C = new TCanvas();

  // partons
  tr->Draw("DP2P_PtFraction:DP2P_Dr>>dp2p_pt_r");
  TH1F *hDiv = (TH1F*)gDirectory->Get("dp2p_pt_r");
  hDiv->SetMarkerStyle(kFullCircle);
  hDiv->SetMarkerColor(1);

  hDiv->GetYaxis()->SetTitle("ptFrac");
  hDiv->GetXaxis()->SetTitle("#delta_R");

  TH1 *h1 = hDiv->DrawCopy();
}

void j2jp_pt_r(TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();


  auto C = new TCanvas();

  // partons
  tr->Draw("J2JP_PtFraction:J2JP_Dr>>j2jp_pt_r");
  TH1F *hDiv = (TH1F*)gDirectory->Get("j2jp_pt_r");
  hDiv->SetMarkerStyle(kFullCircle);
  hDiv->SetMarkerColor(1);

  hDiv->GetYaxis()->SetTitle("ptFrac");
  hDiv->GetXaxis()->SetTitle("#delta_R");

  TH1 *h1 = hDiv->DrawCopy();
}



void dp2p_r(TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();


  auto C = new TCanvas();

  // partons
  tr->Draw("DP2P_Dr>>dp2p_r","DP2P_Dr<1");
  TH1F *hDiv = (TH1F*)gDirectory->Get("dp2p_r");
  hDiv->SetMarkerStyle(kFullCircle);
  hDiv->SetMarkerColor(1);

  hDiv->GetYaxis()->SetTitle("Count");
  hDiv->GetXaxis()->SetTitle("#delta_R");

  TH1 *h1 = hDiv->DrawCopy();
}

void dp2dp_r(TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();


  auto C = new TCanvas();

  // partons
  tr->Draw("DP2DP_Dr>>dp2dp_r","DP2DP_Dr<1");
  TH1F *hDiv = (TH1F*)gDirectory->Get("dp2dp_r");
  hDiv->SetMarkerStyle(kFullCircle);
  hDiv->SetMarkerColor(1);

  hDiv->GetYaxis()->SetTitle("Count");
  hDiv->GetXaxis()->SetTitle("#delta_R");

  TH1 *h1 = hDiv->DrawCopy();
}


void FF(TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();


  auto C = new TCanvas();

  // partons
  tr->Draw("jetsWithQuarkParentsConstPt/jetsWithQuarkParentsPt>>quark");
  TH1F *hDiv = (TH1F*)gDirectory->Get("quark");
  hDiv->SetLineColor(1);
  hDiv->SetLineWidth(3);
  hDiv->SetTitle("FF by jetParent");
  // jets
  tr->Draw("jetsWithGluonParentsConstPt/jetsWithGluonParentsPt>>gluon");
  TH1F *hDiv2 = (TH1F*)gDirectory->Get("gluon");
  hDiv2->SetLineColor(2);
  hDiv2->SetLineWidth(2);
  // daughterpartons
  tr->Draw("otherJetsConstPt/otherJetsPt>>other");
  TH1F *hDiv3 = (TH1F*)gDirectory->Get("other");
  hDiv3->SetLineColor(4);
  hDiv3->SetLineWidth(2);
  hDiv3->SetLineStyle(kDotted);


  hDiv->GetXaxis()->SetTitle("ConstituentPt/JetPt");
  hDiv->GetYaxis()->SetTitle("#");
  gPad->SetLogy();


  TH1 *h1 = hDiv->DrawCopy();
  TH1 *h2 = hDiv2->DrawCopy("same");
  TH1 *h3 = hDiv3->DrawCopy("same");
  
  TLegend *legend = new TLegend(0.55,0.65,0.76,0.82);
  legend->AddEntry(h1,"quarks");
  legend->AddEntry(h2,"gluons");
  legend->AddEntry(h3,"others");
  legend->Draw();
}

void FFbyPt(int minPt=25, int maxPt=30, TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();


  auto C = new TCanvas();

  // partons
  tr->Draw("jetsWithQuarkParentsConstPt/jetsWithQuarkParentsPt>>quark", (TString) ("jetsWithQuarkParentsPt > " + to_string(minPt) + " && jetsWithQuarkParentsPt < " + to_string(maxPt)));
  TH1F *hDiv = (TH1F*)gDirectory->Get("quark");
  hDiv->SetLineColor(1);
  hDiv->SetLineWidth(3);
  hDiv->SetTitle((TString) ("FF by jetParent for " + to_string(minPt) + " < Z < " + to_string(maxPt)));
  // jets
  tr->Draw("jetsWithGluonParentsConstPt/jetsWithGluonParentsPt>>gluon", (TString) ("jetsWithGluonParentsPt > " + to_string(minPt) + " && jetsWithGluonParentsPt < " + to_string(maxPt)));
  TH1F *hDiv2 = (TH1F*)gDirectory->Get("gluon");
  hDiv2->SetLineColor(2);
  hDiv2->SetLineWidth(2);
  // daughterpartons
  tr->Draw("otherJetsConstPt/otherJetsPt>>other", (TString) ("otherJetsPt > " + to_string(minPt) + " && otherJetsPt < " + to_string(maxPt)));
  TH1F *hDiv3 = (TH1F*)gDirectory->Get("other");
  hDiv3->SetLineColor(4);
  hDiv3->SetLineWidth(2);
  hDiv3->SetLineStyle(kDotted);


  hDiv->GetXaxis()->SetTitle("ConstituentPt/JetPt");
  hDiv->GetYaxis()->SetTitle("#");
  gPad->SetLogy();


  TH1 *h1 = hDiv->DrawCopy();
  TH1 *h2 = hDiv2->DrawCopy("same");
  TH1 *h3 = hDiv3->DrawCopy("same");
  
  TLegend *legend = new TLegend(0.55,0.65,0.76,0.82);
  legend->AddEntry(h1,"quarks");
  legend->AddEntry(h2,"gluons");
  legend->AddEntry(h3,"others");
  legend->Draw();
}


// void FF(TString str = "JetToyHIResultSimpleJetAnalysis.root")
// {
//   TTree *tr = getTree();

//   int nEvt = tr->GetEntriesFast();


//   auto C = new TCanvas();

//   // partons
//   tr->Draw("FF_quark>>quark");
//   TH1F *hDiv = (TH1F*)gDirectory->Get("quark");
//   hDiv->SetLineColor(1);
//   hDiv->SetLineWidth(3);
//   hDiv->SetTitle("FF by jetParent");
//   // jets
//   tr->Draw("FF_gluon>>gluon");
//   TH1F *hDiv2 = (TH1F*)gDirectory->Get("gluon");
//   hDiv2->SetLineColor(2);
//   hDiv2->SetLineWidth(2);
//   // daughterpartons
//   tr->Draw("FF_other>>other");
//   TH1F *hDiv3 = (TH1F*)gDirectory->Get("other");
//   hDiv3->SetLineColor(4);
//   hDiv3->SetLineWidth(2);
//   hDiv3->SetLineStyle(kDotted);


//   hDiv->GetXaxis()->SetTitle("ConstituentPt/JetPt");
//   hDiv->GetYaxis()->SetTitle("#");


//   TH1 *h1 = hDiv->DrawCopy();
//   TH1 *h2 = hDiv2->DrawCopy("same");
//   TH1 *h3 = hDiv3->DrawCopy("same");
  
//   TLegend *legend = new TLegend(0.55,0.65,0.76,0.82);
//   legend->AddEntry(h1,"quarks");
//   legend->AddEntry(h2,"gluons");
//   legend->AddEntry(h3,"others");
//   legend->Draw();
// }






void dp_pt_phi(TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();


  auto C = new TCanvas();

  // partons
  tr->Draw("DP2P_PtFraction:daughterPartonsPhi>>dpptphi");
  TH1F *hDiv = (TH1F*)gDirectory->Get("dpptphi");
  hDiv->SetMarkerStyle(kFullCircle);
  hDiv->SetMarkerColor(1);

  hDiv->GetYaxis()->SetTitle("ptFrac");
  hDiv->GetXaxis()->SetTitle("#phi");

  TH1 *h1 = hDiv->DrawCopy();
}

void pt_phi(TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();


  auto C = new TCanvas();

  // partons
  tr->Draw("minPtFracDPsPt/partonsPt:minPtFracDPsPhi>>dps");
  TH1F *hDiv = (TH1F*)gDirectory->Get("dps");
  hDiv->SetMarkerStyle(kFullCircle);
  hDiv->SetMarkerColor(1);

  hDiv->GetYaxis()->SetTitle("min(pt_{dp1}, pt_{dp2})/pt_{m}");
  hDiv->GetXaxis()->SetTitle("#phi");

  TH1 *h1 = hDiv->DrawCopy();
}

void pt_r(TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();


  auto C = new TCanvas();

  // partons
  tr->Draw("minPtFracDPsPt/partonsPt:DP2DP_Dr>>dps");
  TH1F *hDiv = (TH1F*)gDirectory->Get("dps");
  hDiv->SetMarkerStyle(kFullCircle);
  hDiv->SetMarkerColor(1);

  hDiv->GetYaxis()->SetTitle("min(pt_{dp1}, pt_{dp2})/pt_{m}");
  hDiv->GetXaxis()->SetTitle("#delta R");

  TH1 *h1 = hDiv->DrawCopy();
}


void plot(TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();

  //----------------------------------------------------------------
  // plot info from tree
  //----------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1","c1",450,400);
  TH1F *fr1 = DrawFrame(0.,50.,0.,10.,"pt","N");
  tr->Draw("partonsPt","","same");
  tr->SetLineColor(2);
  tr->Draw("daughterPartonsPt","","same");
}

TH1F *DrawFrame(double xmin, double xmax, double ymin, double ymax, TString xTitle, TString yTitle, bool setMargins)
{
  if(setMargins) {
    gPad->SetLeftMargin(0.22);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.1);//0.05);
    gPad->SetTopMargin(0.05);
  }

  TH1F *frame = gPad->DrawFrame(xmin,ymin,xmax,ymax);
  frame->SetXTitle(xTitle.Data());
  frame->SetYTitle(yTitle.Data());
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleOffset(1.3);
  frame->GetXaxis()->CenterTitle(true);
  frame->GetYaxis()->CenterTitle(true);

  gPad->SetTicks(1,1);

  return frame;
}

void histpalettecolor()
{
   auto C = new TCanvas();
 
   gStyle->SetOptTitle(kFALSE);
   gStyle->SetOptStat(0);
 
   auto h1 = new TH1F("h1","Histogram drawn with full circles",100,-4,4);
   auto h2 = new TH1F("h2","Histogram drawn with full squares",100,-4,4);
   auto h3 = new TH1F("h3","Histogram drawn with full triangles up",100,-4,4);
   auto h4 = new TH1F("h4","Histogram drawn with full triangles down",100,-4,4);
   auto h5 = new TH1F("h5","Histogram drawn with empty circles",100,-4,4);
 
   TRandom3 rng;
   Double_t px,py;
   for (Int_t i = 0; i < 25000; i++) {
      rng.Rannor(px,py);
      h1->Fill(px,10.);
      h2->Fill(px, 8.);
      h3->Fill(px, 6.);
      h4->Fill(px, 4.);
      h5->Fill(px, 2.);
   }
 
   h1->SetMarkerStyle(kFullCircle);
   h2->SetMarkerStyle(kFullSquare);
   h3->SetMarkerStyle(kFullTriangleUp);
   h4->SetMarkerStyle(kFullTriangleDown);
   h5->SetMarkerStyle(kOpenCircle);
 
   h1->Draw("PLC PMC");
   h2->Draw("SAME PLC PMC");
   h3->Draw("SAME PLC PMC");
   h4->Draw("SAME PLC PMC");
   h5->Draw("SAME PLC PMC");
 
   gPad->BuildLegend();
}





TTree *getTree(TString str = "JetToyHIResultSimpleJetAnalysis.root", TString tree = "HerjansBoompje")
{
  TFile *f = new TFile(str.Data());
  TTree *tr = dynamic_cast<TTree*>(f->Get(tree.Data()));

  tr->SetLineWidth(2);
  return tr;
}