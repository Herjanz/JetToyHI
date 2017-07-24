#include <iostream>
#include <chrono>

#include "TFile.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "include/ProgressBar.h"

#include "PU14/EventMixer.hh"
#include "PU14/CmdLine.hh"
#include "PU14/PU14.hh"

#include "include/jetCollection.hh"
#include "include/thermalEvent.hh"
#include "include/pythiaEvent.hh"
#include "include/csSubtractor.hh"
#include "include/skSubtractor.hh"
#include "include/softDropGroomer.hh"
#include "include/treeWriter.hh"
#include "include/jetMatcher.hh"
#include "include/randomCones.hh"

using namespace std;
using namespace fastjet;

int main (int argc, char ** argv) {

  auto start_time = std::chrono::steady_clock::now();
  
  CmdLine cmdline(argc,argv);
  // inputs read from command line
  int nEvent = cmdline.value<int>("-nev",1);  // first argument: command line option; second argument: default value
  //bool verbose = cmdline.present("-verbose");

  std::cout << "will run on " << nEvent << " events" << std::endl;

  // Uncomment to silence fastjet banner
  ClusterSequence::set_fastjet_banner_stream(NULL);

  //to write info to root tree
  treeWriter trw("jetTree");

  //Jet definition
  double R                   = 0.4;
  double ghostRapMax         = 6.0;
  double ghost_area          = 0.005;
  int    active_area_repeats = 1;
  fastjet::GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
  fastjet::JetDefinition jet_def(antikt_algorithm, R);

  double jetRapMax = 3.0;
  fastjet::Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  ProgressBar Bar(cout, nEvent);
  Bar.SetStyle(-1);

  EventMixer mixer(&cmdline);  //the mixing machinery from PU14 workshop

  // loop over events
  int iev = 0;
  unsigned int entryDiv = (nEvent > 200) ? nEvent / 200 : 1;
  while ( mixer.next_event() && iev < nEvent ) {
    // increment event number    
    iev++;

    Bar.Update(iev);
    Bar.PrintWithMod(entryDiv);

    std::vector<fastjet::PseudoJet> particlesMerged = mixer.particles() ;

    // cluster hard event only
    std::vector<fastjet::PseudoJet> particlesBkg, particlesSig;
    //particlesBkg = particlesMerged;
    //particlesSig = particlesMerged;
    SelectorIsHard().sift(particlesMerged, particlesSig, particlesBkg); // this sifts the full event into two vectors of PseudoJet, one for the hard event, one for the underlying event
    
    //---------------------------------------------------------------------------
    //   jet clustering
    //---------------------------------------------------------------------------
    
    // run the clustering, extract the jets
    fastjet::ClusterSequenceArea csMerged(particlesMerged, jet_def, area_def);
    jetCollection jetCollectionMerged(sorted_by_pt(jet_selector(csMerged.inclusive_jets())));

    // fastjet::ClusterSequenceArea csBkg(particlesBkg, jet_def, area_def);
    // jetCollection jetCollectionBkg(sorted_by_pt(csBkg.inclusive_jets()));

    fastjet::ClusterSequenceArea csSig(particlesSig, jet_def, area_def);
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(csSig.inclusive_jets(10.))));

    //---------------------------------------------------------------------------
    //   background subtraction
    //---------------------------------------------------------------------------
    
    //run constituent subtraction on mixed (hard+UE) event
    csSubtractor csSub(R, 1., -1, 0.005,ghostRapMax,jetRapMax);
    //csSub.setInputParticles(particlesMerged);
    csSub.setInputJets(jetCollectionMerged.getJet());
    jetCollection jetCollectionCS(csSub.doSubtraction());

    //Background densities used by constituent subtraction
    std::vector<double> rho;    rho.push_back(csSub.getRho());
    std::vector<double> rhom;   rhom.push_back(csSub.getRhoM());

    //Uncomment if youw ant to study random cones
    // randomCones rc(4,R,2.3,rho[0]);
    // rc.setInputParticles(particlesMerged);
    // jetCollection jetCollectionRC(rc.run());

    //run soft killer on mixed event
    skSubtractor skSub(0.4, 3.0);
    skSub.setInputParticles(particlesMerged);
    std::vector<fastjet::PseudoJet> skEvent = skSub.doSubtraction();
    std::vector<double> skPtThreshold;
    skPtThreshold.push_back(skSub.getPtThreshold()); //SoftKiller pT threshold

    fastjet::ClusterSequenceArea csSK(skEvent, jet_def, area_def);
    jetCollection jetCollectionSK(sorted_by_pt(jet_selector(csSK.inclusive_jets())));

    //std::cout << "njets CS: " << jetsCS.size() << " SK: " << jetsSK.size() << std::endl;
    
    //---------------------------------------------------------------------------
    //   Groom the jets
    //---------------------------------------------------------------------------
    
    //SoftDrop grooming classic for CS jets (zcut=0.1, beta=0)
    softDropGroomer sdgCS(0.1, 0.0, R);
    jetCollection jetCollectionCSSD(sdgCS.doGrooming(jetCollectionCS));
    jetCollectionCSSD.addVector("zgCSSD",    sdgCS.getZgs());
    jetCollectionCSSD.addVector("ndropCSSD", sdgCS.getNDroppedBranches());

    //SoftDrop grooming classic for signal jets (zcut=0.1, beta=0)
    softDropGroomer sdgSig(0.1, 0.0, R);
    jetCollection jetCollectionSigSD(sdgSig.doGrooming(jetCollectionSig));
    jetCollectionSigSD.addVector("zgSigSD",    sdgSig.getZgs());
    jetCollectionSigSD.addVector("ndropSigSD", sdgSig.getNDroppedBranches());

    //match the CS jets to signal jets
    jetMatcher jmCS(R);
    jmCS.setBaseJets(jetCollectionCS);
    jmCS.setTagJets(jetCollectionSig);
    jmCS.matchJets();

    jmCS.reorderedToTag(jetCollectionCS);
    jmCS.reorderedToTag(jetCollectionCSSD);

    //match the SK jets to signal jets
    jetMatcher jmSK(R);
    jmSK.setBaseJets(jetCollectionSK);
    jmSK.setTagJets(jetCollectionSig);
    jmSK.matchJets();
    
    jmSK.reorderedToTag(jetCollectionSK);

    //match the unsubtracted jets to signal jets
    jetMatcher jmUnSub(R);
    jmUnSub.setBaseJets(jetCollectionMerged);
    jmUnSub.setTagJets(jetCollectionSig);
    jmUnSub.matchJets();

    jmUnSub.reorderedToTag(jetCollectionMerged);

    //---------------------------------------------------------------------------
    //   write tree
    //---------------------------------------------------------------------------
    
    //Give variable we want to write out to treeWriter.
    //Only vectors of the types 'jetCollection', and 'double', 'int', 'fastjet::PseudoJet' are supported

    trw.addCollection("sigJet",        jetCollectionSig);
    trw.addCollection("csJet",         jetCollectionCS);
    trw.addCollection("sigJetSD",      jetCollectionSigSD);
    trw.addCollection("csJetSD",       jetCollectionCSSD);
    trw.addCollection("skJet",         jetCollectionSK);
    //trw.addCollection("randomCones",   jetCollectionRC);

    trw.addCollection("csRho",         rho);
    trw.addCollection("csRhom",        rhom);
    trw.addCollection("skPtThreshold", skPtThreshold);

    trw.addCollection("unsubJet",      jetCollectionMerged);
        
    trw.fillTree();

  }//event loop

  Bar.Update(nEvent);
  Bar.Print();
  Bar.PrintLine();

  TTree *trOut = trw.getTree();

  TFile *fout = new TFile("JetToyHIResultFromFile.root","RECREATE");
  trOut->Write();
  fout->Write();
  fout->Close();

  double time_in_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time).count() / 1000.0;
  std::cout << "runFromFile: " << time_in_seconds << std::endl;
}