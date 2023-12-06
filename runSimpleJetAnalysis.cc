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

#include "include/ParticleToParticle.hh"

using namespace std;
using namespace fastjet;

// ./runSimpleJetAnalysis -hard /Users/mverweij/mnt/eos/project/j/jetquenching/JetWorkshop2017/samples/pythia8/dijet120/PythiaEventsTune14PtHat120_0.pu14 -nev 10

bool containsHadrons(vector<PseudoJet> particles);

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
  // double R                   = 0.4;
  // double ghostRapMax         = 6.0;
  // double ghost_area          = 0.005;
  // int    active_area_repeats = 1;
  // GhostedAreaSpec ghost_spec(ghostRapMax, active_area_repeats, ghost_area);
  // AreaDefinition area_def = AreaDefinition(active_area,ghost_spec);
  double ycut = 0.02;
  JetDefinition jet_def(ee_kt_algorithm);

  double jetRapMax = 3.0;
  Selector jet_selector = SelectorAbsRapMax(jetRapMax);

  //Angularity width(1.,1.,R);
  //Angularity pTD(0.,2.,R);

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

    fastjet::ClusterSequence clust_seq(particlesSig, jet_def);
    int n = clust_seq.n_exclusive_jets_ycut(ycut); // get 3 exclusive jets
    jetCollection jetCollectionSig(sorted_by_pt(jet_selector(clust_seq.exclusive_jets(n))));

    //calculate some angularities & fragmentations
    /*vector<double> widthSig; widthSig.reserve(jetCollectionSig.getJet().size());
    vector<double> pTDSig;   pTDSig.reserve(jetCollectionSig.getJet().size());
    for(PseudoJet jet : jetCollectionSig.getJet()) {
      widthSig.push_back(width.result(jet));
      pTDSig.push_back(pTD.result(jet));
    }
    jetCollectionSig.addVector("widthSig", widthSig);
    jetCollectionSig.addVector("pTDSig", pTDSig);
    */ // only works with old clustering algorithm based on R



    // ------ START Herjans code ------

    // --- Find splitted particles from partons that can be viewed as originators for jets --- \\
    
    // jetParents are partons or daughterPartons. Allowed to be parents for jets
    vector<PseudoJet> jetParents;

    // get first splitted gluons and handle them like jet parents if they pass conditions [minimum momentum/angle]
    // extract hard partons from first splitting
    fastjet::Selector parton_selector_split = SelectorVertexNumber(-2);
    vector<PseudoJet> daughterPartons = parton_selector_split(particlesMergedAll);

    trw.addPartonCollection("daughterPartons", daughterPartons);
    // skip this event if daughters are already hadronized
    if(containsHadrons(daughterPartons))
      continue;
    
    // daughterparton[0/1] belong to parton[0] and dp[2/3] belong to p[1]
    // get dr between daughter partons
    vector<double> DP2DP_Dr;
    DP2DP_Dr.push_back(daughterPartons[0].delta_R(daughterPartons[1]));
    DP2DP_Dr.push_back(daughterPartons[2].delta_R(daughterPartons[3]));
    // get ptFraction between daughterPartons and partons
    vector<double> DP2P_PtFraction;
    DP2P_PtFraction.push_back(daughterPartons[0].pt()/partons[0].pt());
    DP2P_PtFraction.push_back(daughterPartons[1].pt()/partons[0].pt());
    DP2P_PtFraction.push_back(daughterPartons[2].pt()/partons[1].pt());
    DP2P_PtFraction.push_back(daughterPartons[3].pt()/partons[1].pt());

    vector<double> DP2P_Dr;
    DP2P_Dr.push_back(daughterPartons[0].delta_R(partons[0]));
    DP2P_Dr.push_back(daughterPartons[1].delta_R(partons[0]));
    DP2P_Dr.push_back(daughterPartons[2].delta_R(partons[1]));
    DP2P_Dr.push_back(daughterPartons[3].delta_R(partons[1]));
    
    // daughter parton to daughter parton graphs
    trw.addDoubleCollection("DP2DP_Dr", DP2DP_Dr);
    trw.addDoubleCollection("DP2P_PtFraction", DP2P_PtFraction);
    trw.addDoubleCollection("DP2P_Dr", DP2P_Dr);

    vector<PseudoJet> minPtFracDPs; // DaughterPartons that min(Pt_d1, Pt_d2)
    if(daughterPartons[0].pt() < daughterPartons[1].pt())
      minPtFracDPs.push_back(daughterPartons[0]);
    else
      minPtFracDPs.push_back(daughterPartons[1]);
    if(daughterPartons[2].pt() < daughterPartons[3].pt())
      minPtFracDPs.push_back(daughterPartons[2]);
    else
      minPtFracDPs.push_back(daughterPartons[3]);
    
    trw.addPartonCollection("minPtFracDPs", minPtFracDPs);

    // ----- select jetParents -----

    // if dr between DPs is very small, ignore them because they belong to one jet
    int validDPs = 0;
    for (unsigned int i = 0; i < DP2DP_Dr.size(); ++i)
    {
      int dpsAdded = 0;
      double minDr = 0.1;
      if(DP2DP_Dr[i] > minDr)
      {
        if(DP2P_Dr[i*2] > minDr && daughterPartons[i*2].pt() > 4)
        {
          jetParents.push_back(daughterPartons[i*2]);
          dpsAdded++;
        }
        if(DP2P_Dr[i*2+1] > minDr && daughterPartons[i*2+1].pt() > 4)
        {
          jetParents.push_back(daughterPartons[i*2+1]);
          dpsAdded++;
        }
      }
      if(dpsAdded < 2)
      {
        jetParents.push_back(partons[i]);
      }
      validDPs += dpsAdded;
    }

    // ----- END select jetParents -----

    
    // Jet to Jet parents
    vector<PtoPInfo> JtoJPmatches = PtoP.findMatches(jetParents, jetCollectionSig.getJet(), true);
    // jet to jetParent graphs
    trw.addDoubleCollection("J2JP_Dr", getDrVector(JtoJPmatches));
    trw.addDoubleCollection("J2JP_PtFraction", getPtFracVector(JtoJPmatches));
    
    trw.addPartonCollection("jetParents", jetParents);

    // sorting jets by parent PDG
    vector<PseudoJet> jetsWithQuarkParents, jetsWithGluonParents, otherJets;
    for(PtoPInfo match : JtoJPmatches)
    {
      int parentPDG = getPDG(match.in);

      vector<PseudoJet> jet;
      jet.push_back(match.out);

      if(parentPDG > -7 && parentPDG < 7)
      {
        jetsWithQuarkParents.push_back(match.out);
      }
      else if(parentPDG == 21)
      {
        jetsWithGluonParents.push_back(match.out);
      }
      else
      {
        cout << parentPDG << endl;
        otherJets.push_back(match.out);
      }
    }

    trw.addJetCollection("jetsWithQuarkParents", jetsWithQuarkParents, true);
    trw.addJetCollection("jetsWithGluonParents", jetsWithGluonParents, true);
    trw.addJetCollection("otherJets", otherJets, true);

    // ---- ANALYZE ----
    bool doPrintInfo = false;
    // take exceptions separately to analyze them.
    // in this case jets with very high energy compared to parent
    for(double frac : getPtFracVector(JtoJPmatches))
    {
      if(frac > 2)
      {
        cout << "validPartons: " << partons.size() << " validDPs: " << validDPs << " jetParents: " << jetParents.size() << endl;

        for(unsigned int i = 0; i < getPtFracVector(JtoJPmatches).size(); i++)
          cout << "JtoJPfrac: " << getPtFracVector(JtoJPmatches)[i] << " JtoJPdr: " << getDrVector(JtoJPmatches)[i] << endl;

        cout << "DP2DP 0 to 1: "<< DP2DP_Dr[0] << endl;
        cout << "DP2DP 2 to 3: "<< DP2DP_Dr[1] << endl;

        for(int i = 0; i < 4; i++)
        {
          cout << "DP2P dr: " << DP2P_Dr[i] << " frac: " << DP2P_PtFraction[i] << endl;
        }
        doPrintInfo = true;
        break;
      }
    }


    if(doPrintInfo)
    {
      for(PtoPInfo i : JtoJPmatches)
      {
        printInfo(i);
      }
    }

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

bool containsHadrons(vector<PseudoJet> particles)
{
  for(unsigned int i = 0; i < particles.size(); i++)
  {
    int daughterPDG = getPDG(particles[i]);
    if(daughterPDG > 99 || daughterPDG < -99)
      return true;
  }
  return false;
}