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
#include "include/ParticleToParticle.hh"

using namespace std;
using namespace fastjet;

// ./runSimpleJetAnalysis -hard /Users/mverweij/mnt/eos/project/j/jetquenching/JetWorkshop2017/samples/pythia8/dijet120/PythiaEventsTune14PtHat120_0.pu14 -nev 10

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
  ParticleToParticle PtoP;
    
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




    // ------ START Herjans code ------

    // --- Find splitted particles from partons that can be viewed as originators for jets --- \\
    
    // jetParents are partons or daughterPartons. Allowed to be parents for jets
    vector<PseudoJet> jetParents;

    // get first splitted gluons and handle them like jet parents if they pass conditions [minimum momentum/angle]
    // extract hard partons from first splitting
    fastjet::Selector parton_selector_split = SelectorVertexNumber(-2);
    vector<PseudoJet> daughterPartons = parton_selector_split(particlesMergedAll);

    trw.addPartonCollection("daughterPartons", daughterPartons);
    
    // daughterparton[0/1] belong to parton[0] and dp[2/3] belong to p[1]
    // get dr between daughter partons
    vector<double> DPdr;
    DPdr.push_back(daughterPartons[0].delta_R(daughterPartons[1]));
    DPdr.push_back(daughterPartons[2].delta_R(daughterPartons[3]));
    // get ptFraction between daughterPartons and partons
    vector<double> DPptFrac;
    DPptFrac.push_back(daughterPartons[0].pt()/partons[0].pt());
    DPptFrac.push_back(daughterPartons[1].pt()/partons[0].pt());
    DPptFrac.push_back(daughterPartons[2].pt()/partons[1].pt());
    DPptFrac.push_back(daughterPartons[3].pt()/partons[1].pt());
    
    // daughter parton to daughter parton graphs
    trw.addDoubleCollection("DP2DP_Dr", DPdr);
    trw.addDoubleCollection("DP2P_PtFraction", DPptFrac);

    vector<PseudoJet> minPtFracDPs; // DaughterPartons that min(Pt_d1, Pt_d2)
    if(DPptFrac[0] < DPptFrac[1])
      minPtFracDPs.push_back(daughterPartons[0]);
    else
      minPtFracDPs.push_back(daughterPartons[1]);
    if(DPptFrac[2] < DPptFrac[3])
      minPtFracDPs.push_back(daughterPartons[2]);
    else
      minPtFracDPs.push_back(daughterPartons[3]);
    
    trw.addPartonCollection("minPtFracDPs", minPtFracDPs);

    // ----- select jetParents -----

    // if dr between DPs is very small, ignore them because they belong to one jet
    int validDPs = 0;
    for (int i = 0; i < DPdr.size(); ++i)
    {
      jetParents.push_back(partons[i]);
      if(DPdr[i] > 10.2)
      {
        jetParents.push_back(daughterPartons[i*2]);
        jetParents.push_back(daughterPartons[i*2+1]);
        validDPs += 2;
      }
    }

    // ----- END select jetParents -----

    cout << "validPartons: " << partons.size() << " validDPs: " << validDPs << " jetParents: " << jetParents.size() << endl;

    // Jet to Jet parents
    vector<PtoPInfo> JtoJPmatches = PtoP.findMatches(jetParents, jetCollectionSig.getJet(), false);
    // jet to jetParent graphs
    trw.addDoubleCollection("jet2Parent_Dr", getDrVector(JtoJPmatches));
    trw.addDoubleCollection("jet2Parent_PtFraction", getPtFracVector(JtoJPmatches));
    
    trw.addPartonCollection("jetParents", jetParents);

    // ------ END Herjans code ------
    

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