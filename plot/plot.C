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

int getNumberOfValuesInRange(TTree *tr, TString field, double start, double end, TString conditions = "")
{
  tr->Draw((TString) (field + ">>f(100, 0., 100.)"), conditions);
  TH1F *hist = (TH1F*)gDirectory->Get("f");
  int binLow = hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->FindBin(start));
  int binHigh = hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->FindBin(end))+1;
  int number = hist->Integral(binLow, binHigh);
  //cout << number << endl;
  return number;
}


void FF(double minPt=0., double maxPt=50., TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree();

  int nEvt = tr->GetEntriesFast();

  auto C = new TCanvas();

  TString quark_conditions = (TString) ("jetsWithQuarkParentsPt > " + to_string(minPt) + " && jetsWithQuarkParentsPt < " + to_string(maxPt));
  TString gluon_conditions = (TString) ("jetsWithGluonParentsPt > " + to_string(minPt) + " && jetsWithGluonParentsPt < " + to_string(maxPt));
  TString other_conditions = (TString) ("otherJetsPt > " + to_string(minPt) + " && otherJetsPt < " + to_string(maxPt));

  // get number of jets for scaling
  int nQuarkJets = getNumberOfValuesInRange(tr, "jetsWithQuarkParentsPt", minPt, maxPt, quark_conditions);
  int nGluonJets = getNumberOfValuesInRange(tr, "jetsWithGluonParentsPt", minPt, maxPt, gluon_conditions);
  int nOtherJets = getNumberOfValuesInRange(tr, "otherJetsPt", minPt, maxPt, other_conditions);

  // partons
  tr->Draw("jetsWithQuarkParentsConstPt/jetsWithQuarkParentsPt>>quark(20, 0., 1.)", quark_conditions);
  TH1F *hDiv = (TH1F*)gDirectory->Get("quark");
  hDiv->SetLineColor(1);
  hDiv->SetLineWidth(3);
  hDiv->SetTitle((TString) ("FF for " + to_string(int(minPt)) + " < jetPt < " + to_string(int(maxPt))));
  hDiv->Scale(1./nQuarkJets, "width");
  //hDiv->Sumw2(); already created due to Scale()
  // jets
  tr->Draw("jetsWithGluonParentsConstPt/jetsWithGluonParentsPt>>gluon(20, 0., 1.)", gluon_conditions);
  TH1F *hDiv2 = (TH1F*)gDirectory->Get("gluon");
  hDiv2->SetLineColor(2);
  hDiv2->SetLineWidth(2);
  hDiv2->Scale(1./nGluonJets, "width");
  // daughterpartons
  tr->Draw("otherJetsConstPt/otherJetsPt>>other(20, 0., 1.)", other_conditions );
  TH1F *hDiv3 = (TH1F*)gDirectory->Get("other");
  hDiv3->SetLineColor(4);
  hDiv3->SetLineWidth(2);
  hDiv3->SetLineStyle(kDotted);
  hDiv3->Scale(1./nOtherJets, "width");


  hDiv->GetXaxis()->SetTitle("ConstituentPt/JetPt");
  hDiv->GetYaxis()->SetTitle("1/Njets dN/dz");
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

void FF_HEP(double minPt=4.5, double maxPt=45., TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree(str);
  //int nEvt = tr->GetEntriesFast();

  auto C = new TCanvas();

  string conditions = "jetsWithQuarkParentsPt > " + to_string(minPt) + " && jetsWithQuarkParentsPt < " + to_string(maxPt);
  conditions += " && abs(jetsWithQuarkParentsEta) <= 1.15";
  TString quark_conditions = (TString) conditions;
  conditions = "jetsWithGluonParentsPt > " + to_string(minPt) + " && jetsWithGluonParentsPt < " + to_string(maxPt);
  conditions += " && abs(jetsWithGluonParentsEta) <= 1.15";
  TString gluon_conditions = (TString) conditions;

  // get number of jets for scaling
  int nQuarkJets = getNumberOfValuesInRange(tr, "jetsWithQuarkParentsPt", minPt, maxPt, quark_conditions);
  int nGluonJets = getNumberOfValuesInRange(tr, "jetsWithGluonParentsPt", minPt, maxPt, gluon_conditions);
  
  const Int_t nBins = 8;
  const Double_t binEdges[nBins+1] = {0, 0.05, 0.1, 0.15, 0.25, 0.35, 0.55, 0.8, 1.2};

  // quark
  tr->Draw("jetsWithQuarkParentsConstPt/jetsWithQuarkParentsPt>>quark(16, 0., 0.8)", quark_conditions);
  TH1F *hDiv = (TH1F*)gDirectory->Get("quark");
  hDiv->SetBins(nBins, binEdges);
  hDiv->SetLineColor(1);
  hDiv->SetLineWidth(2);
  hDiv->SetTitle((TString) ("FF for " + to_string(int(minPt)) + " < jetPt < " + to_string(int(maxPt))));
  hDiv->Scale(1./nQuarkJets, "width");
  //hDiv->Sumw2(); already created due to Scale()

  // gluon
  tr->Draw("jetsWithGluonParentsConstPt/jetsWithGluonParentsPt>>gluon(16, 0., 0.8)", gluon_conditions);
  TH1F *hDiv2 = (TH1F*)gDirectory->Get("gluon");
  hDiv2->SetBins(nBins, binEdges);
  hDiv2->SetLineColor(2);
  hDiv2->SetLineWidth(2);
  hDiv2->Scale(1./nGluonJets, "width");


  // File from research paper to compare with
  TFile f("HEPData-ins467225.root");
  const Int_t nBinsHEP = 7;

  TH1F *hDiv3 = (TH1F*)f.Get("Table 6/Hist1D_y1");
  hDiv3->SetLineColor(1);
  hDiv3->SetLineWidth(2);
  hDiv3->SetLineStyle(7);
  TH1F *hDiv3_e1_stat = (TH1F*)f.Get("Table 6/Hist1D_y1_e1");
  TH1F *hDiv3_e2_stat = (TH1F*)f.Get("Table 6/Hist1D_y1_e2");
  for(int i = 1; i <= nBinsHEP; i++)
  {
    hDiv3->SetBinError(i, sqrt(hDiv3_e1_stat->GetBinContent(i)*hDiv3_e1_stat->GetBinContent(i) + hDiv3_e2_stat->GetBinContent(i)*hDiv3_e2_stat->GetBinContent(i)));
  }
  //hDiv3->SetError(hDiv3_e_stat->GetEntries());

  TH1F *hDiv4 = (TH1F*)f.Get("Table 6/Hist1D_y2");
  hDiv4->SetLineColor(2);
  hDiv4->SetLineWidth(2);
  hDiv4->SetLineStyle(7);
  TH1F *hDiv4_e1_stat = (TH1F*)f.Get("Table 6/Hist1D_y2_e1");
  TH1F *hDiv4_e2_stat = (TH1F*)f.Get("Table 6/Hist1D_y2_e2");
  for(int i = 1; i <= nBinsHEP; i++)
  {
    hDiv4->SetBinError(i, sqrt(hDiv4_e1_stat->GetBinContent(i)*hDiv4_e1_stat->GetBinContent(i) + hDiv4_e2_stat->GetBinContent(i)*hDiv4_e2_stat->GetBinContent(i)));
  }
  // END HEPdata

  // scale pythia data with HEP data
  if(false)
  {
    const int bin = 3;
    int numberPytQuark = hDiv->Integral(bin, bin);
    int numberHEPQuark = hDiv3->Integral(bin, bin);
    hDiv->Scale(1.0 * numberHEPQuark/numberPytQuark);


    int numberPytGluon = hDiv2->Integral(bin, bin);
    int numberHEPGluon = hDiv4->Integral(bin, bin);
    hDiv2->Scale(1.0 * numberHEPGluon/numberPytGluon);
  }
  // end of scaling

  hDiv->GetXaxis()->SetTitle("ConstituentPt/JetPt");
  hDiv->GetYaxis()->SetTitle("1/Njets dN/dz");
  hDiv->GetYaxis()->SetRangeUser(0.01, 600);
  gPad->SetLogy();


  TH1 *h1 = hDiv->DrawCopy("");
  TH1 *h2 = hDiv2->DrawCopy("same");
  TH1 *h3 = hDiv3->DrawCopy("same");
  TH1 *h4 = hDiv4->DrawCopy("same");
  
  TLegend *legend = new TLegend(0.55,0.65,0.76,0.82);
  legend->AddEntry(h1,"quarks");
  legend->AddEntry(h2,"gluons");
  legend->AddEntry(h3,"HEP quarks");
  legend->AddEntry(h4,"HEP gluons");
  legend->Draw();
}




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
  return tr;
}