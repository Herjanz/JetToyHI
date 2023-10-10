#include <iostream>
#include <chrono>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "PU14/EventMixer.hh"
#include "PU14/CmdLine.hh"
#include "PU14/PU14.hh"

#include "include/extraInfo.hh"
#include "include/jetCollection.hh"
#include "include/softDropGroomer.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/Angularity.hh"

#include "include/Fragmentation.hh"

using namespace std;
using namespace fastjet;

// ./runSimpleJetAnalysis -hard /Users/mverweij/mnt/eos/project/j/jetquenching/JetWorkshop2017/samples/pythia8/dijet120/PythiaEventsTune14PtHat120_0.pu14 -nev 10

double getDr(PseudoJet a, PseudoJet b)
{
  double dphi = a.phi() - b.phi();
  if(dphi > pi)
    dphi -= 2*pi;
  if(dphi < -pi)
    dphi += 2*pi;
  double deta = a.eta() - b.eta();
  return sqrt(dphi*dphi + deta*deta);
}

int main (int argc, char ** argv) {

  auto start_time = chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

  cout << "will run on " << nEvent << " events" << endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trw("jetTree");

  //Jet definition
  double R                   = 0.4;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  AreaDefinition area_def = AreaDefinition(active_area,ghost_spec);
  JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 3.0;
  Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  Angularity width(1.,1.,R);
  Angularity pTD(0.,2.,R);

  Fragmentation PtZ;
    
  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  //loop over events
  int iev = 0;
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  while ( mixer.next_event() && iev < nEvent )
  {
    // increment event number    
    iev++;

    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    vector<PseudoJet> particlesMergedAll = mixer.particles();

    vector<double> eventWeight;
    eventWeight.push_back(mixer.hard_weight());
    eventWeight.push_back(mixer.pu_weight());

    // extract hard partons that initiated the jets
    fastjet::Selector parton_selector = SelectorVertexNumber(-1);
    vector<PseudoJet> partons = parton_selector(particlesMergedAll);
    //cout << partons.size() << endl;

    // select final state particles from hard event only
    fastjet::Selector sig_selector = SelectorVertexNumber(0);
    vector<PseudoJet> particlesSig = sig_selector(particlesMergedAll);

    // select final state particles from background event only
    fastjet::Selector bkg_selector = SelectorVertexNumber(1);
    vector<PseudoJet> particlesBkg = bkg_selector(particlesMergedAll);

    vector<PseudoJet> particlesMerged = particlesBkg;
    particlesMerged.insert( particlesMerged.end(), particlesSig.begin(), particlesSig.end() );
    
    
    //vector<PseudoJet> particlesBkg, particlesSig;
    //SelectorIsHard().sift(particlesMerged, particlesSig, particlesBkg); // this sifts the full event into two vectors of PseudoJet, one for the hard event, one for the underlying event

    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    //calculate some angularities & fragmentations
    vector<double> widthSig; widthSig.reserve(jetCollectionSig.getJet().size());
    vector<double> pTDSig;   pTDSig.reserve(jetCollectionSig.getJet().size());
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      widthSig.push_back(width.result(jet));
      pTDSig.push_back(pTD.result(jet));
    }
    jetCollectionSig.addVector("widthSig", widthSig);
    jetCollectionSig.addVector("pTDSig", pTDSig);


    // ------ START Herjans code ------ \\

    // --- Find splitted particles from partons that can be viewed as originators for jets --- \\
    
    // originators are particles after first split of partons that pass checks. Thus allowed to be parents for jets
    vector<PseudoJet> originators;

    // Originator to parent parton info
    vector<double> OriginatorPartonDr;// OriginatorPartonDr.reserve(jetCollectionSig.getJet().size());
    vector<double> OriginatorPartonPtFraction;// OriginatorPartonPtFraction.reserve(jetCollectionSig.getJet().size());
    vector<int> OriginatorPartonPDG;// OriginatorPartonPDG.reserve(jetCollectionSig.getJet().size());
    
    // get first splitted gluons and handle them like jet parents if they pass conditions [minimum momentum/angle]
    // extract hard partons from first splitting
    fastjet::Selector parton_selector_split = SelectorVertexNumber(-2);
    vector<PseudoJet> firstSplittedParticles = parton_selector_split(particlesMergedAll);
    for(PseudoJet jet : firstSplittedParticles) {
      // get distance, to find to which parton jet belongs
      PseudoJet closestParton;
      double shortestDr = 9;
      double PtFraction = 0;
      for (PseudoJet parton : partons)
      {
        double dr = getDr(parton, jet);
        double PtFrac = jet.pt() / parton.pt(); // child pt can not be higher than mother pt
        if(dr < shortestDr && PtFrac <= 1)
        {
          shortestDr = dr;
          closestParton = parton;
          PtFraction = PtFrac;
        }
      }

      // only attach jets to parton if they have an obvious parent, if dr is too big, its uncertain
      if(shortestDr < pi && PtFraction > 0.1)
      {
        // checks passed, add splitted particle as possible parent for jets
        originators.push_back(jet);
      }

      // add originator to parton info
      OriginatorPartonDr.push_back(shortestDr);
      OriginatorPartonPtFraction.push_back(PtFraction);
      OriginatorPartonPDG.push_back(jet.user_info<PU14>().pdg_id());
    }

    // originator to parton info
    trw.addDoubleCollection("OriginatorPartonDr", OriginatorPartonDr);
    trw.addDoubleCollection("OriginatorPartonPtFraction", OriginatorPartonPtFraction);
    trw.addIntCollection("OriginatorPartonPDG", OriginatorPartonPDG);

    // originator info
    trw.addPartonCollection("originators", originators);

    // --- Find origin particles for jets --- \\

    // jet to originator info
    vector<double> JetOriginatorDr;// JetOriginatorDr.reserve(jetCollectionSig.getJet().size());
    vector<double> JetOriginatorPtFraction;
    vector<int> jetOriginatorTypeCounter;// jetParentTypeCounter.reserve(jetCollectionSig.getJet().size());
    // int i = 1; // for individual fragmentation functions [only for individual FF]

    cout << "Originators: " << originators.size() << " Jets: " << jetCollectionSig.getJet().size() << endl;
    
    for(PseudoJet jet : jetCollectionSig.getJet())
    {

      if (!jet.has_constituents())
        continue; // nothing to do with an empty jet

      vector<double> FF = PtZ.getFF(jet);
      
      // individual fragmentation function for this jet [disabled]
      //trw.addCollection("FF" + to_string(iev) + "jet" + to_string(i++), FF);

      trw.addDoubleCollection("FF", FF); // avg in physics

      // get distance, to find to which parton jet belongs
      PseudoJet closestOriginator;
      double shortestDr = 999;
      double PtFraction = 0;
      for (PseudoJet originator : originators)
      {
        double dr = getDr(originator, jet);
        double PtFrac = jet.pt() / originator.pt(); // child pt can not be higher than mother pt

        if(dr < shortestDr && PtFrac <= 1)
        {
          shortestDr = dr;
          closestOriginator = originator;
          PtFraction = PtFrac;
        }
      }

      JetOriginatorDr.push_back(shortestDr);
      JetOriginatorPtFraction.push_back(PtFraction);

      // only attach jets to parton if they have an obvious parent, if dr is too big, its uncertain
      if(shortestDr < 3)
      {
        const int &pdgid = closestOriginator.user_info<PU14>().pdg_id();
        jetOriginatorTypeCounter.push_back(pdgid);
      }

    }
    
    // jet to originator info
    trw.addCollection("JetOriginatorDr", JetOriginatorDr);
    trw.addCollection("JetOriginatorPtFraction", JetOriginatorPtFraction);
    trw.addIntCollection("jetOriginatorTypeCounter", jetOriginatorTypeCounter);

    // ------ END Herjans code ------ \\
    

    //---------------------------------------------------------------------------
    //   Groom the jets
    //---------------------------------------------------------------------------

    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    // softDropGroomer sdgSigBeta00Z01(0.1, 0.0, R);
    // jetCollection jetCollectionSigSDBeta00Z01(sdgSigBeta00Z01.doGrooming(jetCollectionSig));
    // jetCollectionSigSDBeta00Z01.addVector("zgSigSDBeta00Z01",    sdgSigBeta00Z01.getZgs());
    // jetCollectionSigSDBeta00Z01.addVector("ndropSigSDBeta00Z01", sdgSigBeta00Z01.getNDroppedSubjets());
    // jetCollectionSigSDBeta00Z01.addVector("dr12SigSDBeta00Z01",  sdgSigBeta00Z01.getDR12());
  

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'PseudoJet' are supported

    //trw.addCollection("eventWeight",   eventWeight);
    trw.addPartonCollection("partons",       partons);

    trw.addCollection("sigJet",        jetCollectionSig, true);
    //trw.addCollection("sigJetSDBeta00Z01",      jetCollectionSigSDBeta00Z01);
    
    //cout << partons.size() << endl;
    //cout << partons[0].pt() << endl;
    //cout << partons[1].pt() << endl; // both are the same

    //cout << jetCollectionSig.getVectorDouble("FF1jet" + to_string(iev)) << endl;
    trw.fillTree();
  }//event loop

  trw.setTreeName("HerjansBoompje");

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile(cmdline.value<string>("-output", "JetToyHIResultSimpleJetAnalysis.root").c_str(), "RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = chrono::duration_cast<chrono::milliseconds>
    (chrono::steady_clock::now() - start_time).count() / 1000.0;
  cout << "runFromFile: " << time_in_seconds << endl;
}