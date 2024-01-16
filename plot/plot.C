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
  TH1F *hepDiv = (TH1F*)gDirectory->Get("daughterPartons");
  hepDiv->SetMarkerStyle(kFullCircle);
  hepDiv->SetMarkerColor(3);


  hDiv->GetXaxis()->SetTitle("#eta");
  hDiv->GetYaxis()->SetTitle("#phi");


  TH1 *h1 = hDiv->DrawCopy();
  TH1 *h2 = hDiv2->DrawCopy("same");
  TH1 *h3 = hepDiv->DrawCopy("same");
  
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

void FF(bool HEP = true, bool isQuark = false, TString str = "JetToyHIResultSimpleJetAnalysis.root", double Q = 20, double minPt=4.5, double maxPt=45.)
{
  bool FIT = HEP; // remove ! (and keep HEP parameter true) to draw FIT vs ALEPH data

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
  
  const Int_t nBins = 7;
  const Double_t binEdges[nBins+1] = {0, 0.05, 0.1, 0.15, 0.25, 0.35, 0.55, 0.8};

  // quark
  tr->Draw("jetsWithQuarkParentsConstPt/jetsWithQuarkParentsPt>>quark(16, 0., 0.8)", quark_conditions);
  TH1F *simDiv = (TH1F*)gDirectory->Get("quark");
  if(HEP)
    simDiv->SetBins(nBins, binEdges);
  //simDiv->SetTitle((TString) ("FF for " + to_string(int(minPt)) + " < jetPt < " + to_string(int(maxPt))));
  simDiv->Scale(1./nQuarkJets, "width");
  //simDiv->Sumw2(); already created due to Scale()

  // gluon
  tr->Draw("jetsWithGluonParentsConstPt/jetsWithGluonParentsPt>>gluon(16, 0., 0.8)", gluon_conditions);
  TH1F *simDiv2 = (TH1F*)gDirectory->Get("gluon");
  if(HEP)
    simDiv2->SetBins(nBins, binEdges);
  simDiv2->Scale(1./nGluonJets, "width");

  if(!isQuark)
    simDiv = simDiv2;

  simDiv->SetLineColor(1);
  simDiv->SetLineWidth(3);
  // END Pythia


  // HEP
  TH1F *hepDiv, *hepDiv2;
  // File from research paper to compare with
  TFile f("HEPData-ins467225.root");
  const Int_t nBinsHEP = 7;

  hepDiv = (TH1F*)f.Get("Table 6/Hist1D_y1");
  TH1F *hepDiv_e1_stat = (TH1F*)f.Get("Table 6/Hist1D_y1_e1");
  TH1F *hepDiv_e2_stat = (TH1F*)f.Get("Table 6/Hist1D_y1_e2");
  for(int i = 1; i <= nBinsHEP; i++)
  {
    hepDiv->SetBinError(i, sqrt(hepDiv_e1_stat->GetBinContent(i)*hepDiv_e1_stat->GetBinContent(i) + hepDiv_e2_stat->GetBinContent(i)*hepDiv_e2_stat->GetBinContent(i)));
  }
  //hepDiv->SetError(hepDiv_e_stat->GetEntries());

  hepDiv2 = (TH1F*)f.Get("Table 6/Hist1D_y2");
  TH1F *hepDiv2_e1_stat = (TH1F*)f.Get("Table 6/Hist1D_y2_e1");
  TH1F *hepDiv2_e2_stat = (TH1F*)f.Get("Table 6/Hist1D_y2_e2");
  for(int i = 1; i <= nBinsHEP; i++)
  {
    hepDiv2->SetBinError(i, sqrt(hepDiv2_e1_stat->GetBinContent(i)*hepDiv2_e1_stat->GetBinContent(i) + hepDiv2_e2_stat->GetBinContent(i)*hepDiv2_e2_stat->GetBinContent(i)));
  }

  if(!isQuark)
    hepDiv = hepDiv2;

  hepDiv->SetLineColor(2);
  hepDiv->SetLineWidth(3);
  hepDiv->SetLineStyle(7);
  // END HEP


  // FIT
  string hadron = "pi";

  bool QUARK = true;
  TF1 *fitDiv = loadKKP(hadron, QUARK, Q);
  //TH1F *q = new TH1F();//tf_q->DoCreateHistogram(0, 1);
  //q->GetListOfFunctions()->Add(tf_q);
  //q->SetBins(nBins, binEdges);

  TF1 *fitDiv2 = loadKKP(hadron, !QUARK, Q);

  if(!isQuark)
    fitDiv = fitDiv2;

  fitDiv->SetLineColor(3);
  fitDiv->SetLineWidth(3);
  fitDiv->SetLineStyle(2);
  // END FIT

  // scale pythia data with HEP data
  // if(false)
  // {
  //   const int bin = 1;
  //   int numberPytQuark = simDiv->Integral(bin, bin);
  //   int numberHEPQuark = hepDiv->Integral(bin, bin);
  //   hDiv->Scale(1.0 * numberHEPQuark/numberPytQuark);


  //   int numberPytGluon = simDiv2->Integral(bin, bin);
  //   int numberHEPGluon = hepDiv2->Integral(bin, bin);
  //   hDiv2->Scale(1.0 * numberHEPGluon/numberPytGluon);
  // }
  // end of scaling

  if(HEP && FIT)
    simDiv = hepDiv;

  simDiv->GetXaxis()->SetTitle("p_{T,Constituent}/p_{T,Jet}");
  simDiv->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN}{dz}");
  simDiv->GetXaxis()->SetRangeUser(0, 0.8);
  simDiv->GetYaxis()->SetRangeUser(0.01, 1000);
  simDiv->GetXaxis()->SetLabelSize(0.05);
  simDiv->GetYaxis()->SetLabelSize(0.05);
  simDiv->GetXaxis()->SetTitleSize(0.05);
  simDiv->GetYaxis()->SetTitleSize(0.05);
  simDiv->GetXaxis()->SetTitleOffset(1.1);
  simDiv->GetYaxis()->SetTitleOffset(1.1);
  simDiv->SetTitle("");
  simDiv->SetStats(false);
  gPad->SetLogy();

  C->SetLeftMargin(0.13);
  C->SetBottomMargin(0.13);
  C->SetRightMargin(0.05);
  C->SetTopMargin(0.05);



  TLegend *legend = new TLegend(0.4, 0.8, 0.95, 0.95);
  gStyle->SetLegendTextSize(0.05);
  gStyle->SetLegendBorderSize(2);
  
  string type = "quark-initiated jets";
  if(!isQuark)
    type = "gluon-initiated jets";

  if(!(HEP && FIT))
  {
    TH1 *h1 = simDiv->DrawCopy("");
    legend->AddEntry(h1,(TString) ("PYTHIA " + type));
  }

  if(HEP){
    TH1 *h2 = hepDiv->DrawCopy(((HEP && FIT) ? "" : "same"));
    legend->AddEntry(h2,(TString) ("ALEPH " + type));
  }
  if(FIT)
  {
    auto *h3 = fitDiv->DrawCopy("same");
    legend->AddEntry(h3,(TString) ("KKP " + type));
  }

  legend->Draw();
  string pl = "";
  if(HEP)
    pl += "_HEP";
  if(FIT)
    pl += "_FIT";
  string title = type + pl + "_Q=" + to_string(int(Q)) + "_minPt=" + to_string(int(minPt));
  C->SaveAs((TString) ("img/ff/" + (string)(str.Data()) + "_" + title + ".jpg"));
}


void pt_jets(bool quark = true, TString str = "JetToyHIResultSimpleJetAnalysis.root")
{
  TTree *tr = getTree(str);

  int nEvt = tr->GetEntriesFast();

  auto C = new TCanvas();

  string type = "Quark";
  if(!quark)
    type = "Gluon";

  // partons
  tr->Draw((TString) ("jetsWith" + type + "ParentsPt>>pts"));
  TH1F *hDiv = (TH1F*)gDirectory->Get("pts");
  //hDiv->SetMarkerStyle(kFullCircle);
  //hDiv->SetMarkerColor(1);
  hDiv->SetLineWidth(3);

  hDiv->GetYaxis()->SetTitle("N_{jets}");
  hDiv->GetXaxis()->SetTitle("p_{T} GeV");
  
  hDiv->GetXaxis()->SetRangeUser(0., 50.);
  hDiv->GetYaxis()->SetRangeUser(0, (quark ? 400 : 80));
  hDiv->GetXaxis()->SetLabelSize(0.05);
  hDiv->GetYaxis()->SetLabelSize(0.05);
  hDiv->GetXaxis()->SetTitleSize(0.05);
  hDiv->GetYaxis()->SetTitleSize(0.05);
  hDiv->GetXaxis()->SetTitleOffset(1.1);
  hDiv->GetYaxis()->SetTitleOffset(1);
  hDiv->SetTitle("");
  hDiv->SetStats(false);

  C->SetLeftMargin(0.11);
  C->SetBottomMargin(0.13);
  C->SetRightMargin(0.05);
  C->SetTopMargin(0.05);

  // const Int_t nBins = 100;
  // double binEdges[nBins+1];
  // for(int i = 0; i <= nBins; i++)
  //   binEdges[i] = i;
  // hDiv->SetBins(nBins, binEdges);

  TH1 *h1 = hDiv->DrawCopy();

  string title = type + "_" + (string)(str.Data());
  C->SaveAs((TString) ("img/jetpts/" + title + ".jpg"));
}

void plot_all()
{ // makes use of the root files generated for each algorithm, name your files as below to execute this
  std::string roots[5] = {"dot4R.root", "dot8R.root", "1dot2R.root","1dot6R.root", "Durham.root"};
  for(string r : roots)
  {
    pt_jets(true, (TString) r);
    pt_jets(false, (TString) r);
  }

  for(string r : roots)
  {
    FF(true, true, (TString) r);
    FF(true, false, (TString) r);
    FF(false, true, (TString) r);
    FF(false, false, (TString) r);
  }
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