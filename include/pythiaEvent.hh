#ifndef pythiaEvent_h
#define pythiaEvent_h

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <ctime>

#include "Pythia8/Pythia.h"

#include "extraInfo.hh"

//using namespace std;

//---------------------------------------------------------------
// Description
// This class generates a pythia8 event
// Author: M. Verweij
//---------------------------------------------------------------

class pythiaEvent {

private :
  Pythia8::Pythia pythia;
  double pthat_;
  unsigned int tune_;
  double rapMin_;
  double rapMax_;
  bool   partonLevel_;
  bool   vinciaShower_;

  std::vector<fastjet::PseudoJet> partons;

public :
  pythiaEvent(double pthat = 120., unsigned int tune = 14, double rapMin = -3., double rapMax = 3., bool partonLevel = false, bool vinciaShower = false, bool flatPtHat = false);
  std::vector<fastjet::PseudoJet> createPythiaEvent();
  
  std::vector<fastjet::PseudoJet> getPartonList() const { return partons; }

  void getStat() {pythia.stat();}
  double getWeight() {return pythia.info.weight();}
  double getPtHat()  {return pythia.info.pTHat();}

};
  
pythiaEvent::pythiaEvent(double pthat, unsigned int tune, double rapMin, double rapMax, bool partonLevel, bool vinciaShower, bool flatPtHat) :
  pthat_(pthat), tune_(tune), rapMin_(rapMin), rapMax_(rapMax), partonLevel_(partonLevel), vinciaShower_(vinciaShower)
{
    
  // Generator. LHC process and output selection. Initialization.
  // tunes: http://home.thep.lu.se/~torbjorn/pythia82html/Tunes.html
  pythia.readString("Beams:idA = 11.");
  pythia.readString("Beams:idB = -11");
  //pythia.readString("Beams:eCM = 5020.");
  double mZ = pythia.particleData.m0(23);
  pythia.settings.parm("Beams:eCM", mZ);

  // Allow no substructure in e+- beams: normal for corrected LEP data.
  pythia.readString("PDF:lepton = off");
  // Process selection.
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  // Switch off all Z0 decays and then switch back on those to quarks.
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5 6 21");

  //pythia.readString("HardQCD:all = on");
  //pythia.readString(Form("PhaseSpace:pTHatMin = %.1f",pthat_));
  pythia.readString("PhaseSpace:pTHatMin = 0");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  if(partonLevel_) {
    pythia.readString("HadronLevel:all = off");
  }
  if(vinciaShower_)
    pythia.readString("PartonShowers:Model = 2"); //activate the VINCIA parton shower
  else
    pythia.readString(Form("Tune:ee = %d",7));
      //pythia.readString(Form("Tune:pp = %d",tune_));

  if(false) { // flatPtHat
    //flat pthard
    pythia.readString("PhaseSpace:bias2Selection = on");
    pythia.readString("PhaseSpace:bias2SelectionPow = 4.5");
    pythia.readString("PhaseSpace:bias2SelectionRef = 15");
  }
  
  pythia.init();

}

std::vector<fastjet::PseudoJet> pythiaEvent::createPythiaEvent() {

  pythia.next(); //generate next event
    
  std::vector<fastjet::PseudoJet> particles;
  partons.clear(); //empty list before storing partons of new event

  //int iprint = 0;
  
  for (int i = 0; i < pythia.event.size(); ++i) {
    if (pythia.event[i].isFinal()) { //all final state particles
      fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      p.set_user_info(new extraInfo(pythia.event[i].id(), 0)); 
      if(p.rap()>rapMin_ && p.rap()<rapMax_)
        particles.push_back(p);
    } else if (pythia.event[i].status()==-23) { //outgoing partons from the hard scattering
      fastjet::PseudoJet p(pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e());
      p.set_user_info(new extraInfo(pythia.event[i].id(), -1)); 
      partons.push_back(p);
      
      //find the case where the splitting to two separate daughters happens
      int d1 = pythia.event[i].daughter1();
      int d2 = pythia.event[i].daughter2();
      while(d1==d2) {
        d1 = pythia.event[d1].daughter1();
        d2 = pythia.event[d2].daughter2();
      }

      if(vinciaShower_) { //in vincia need to pass the recoil type of particles
        int st1 = pythia.event[d1].status();
        int st2 = pythia.event[d2].status();
        //std::cout << "d1: " << d1 << " st1: " << st1 << " d2: " << d2 << " st2: " << st2 << std::endl;
        
        while(abs(st1)<51 && abs(st2)<51) {
          if(st1==-44) { //44 : outgoing shifted by a branching
            d1 = pythia.event[d1].daughter1();
            d2 = pythia.event[d1].daughter2();
          } else if(st2==-44) {
            d1 = pythia.event[d2].daughter1();
            d2 = pythia.event[d2].daughter2();
          } else {
            //std::cout << "mom: " << i << " d1: " << d1 << " st1: " << pythia.event[d1].status() << " d2: " << d2 << " st2: " << pythia.event[d2].status() << std::endl;
            double pt1 = pythia.event[d1].px()*pythia.event[d1].px() + pythia.event[d1].py()*pythia.event[d1].py();
            double pt2 = pythia.event[d2].px()*pythia.event[d2].px() + pythia.event[d2].py()*pythia.event[d2].py();
            if(pt1 > pt2) {
              d1 = pythia.event[d1].daughter1();
              d2 = pythia.event[d1].daughter2();
            } else {
              d1 = pythia.event[d2].daughter1();
              d2 = pythia.event[d2].daughter2();
            }
            //iprint = 1;
            //std::cout << "mom: " << i << " d1: " << d1 << " st1: " << pythia.event[d1].status() << " d2: " << d2 << " st2: " << pythia.event[d2].status() << std::endl;
          }
          st1 = pythia.event[d1].status();
          st2 = pythia.event[d2].status();

          //std::cout << "d1: " << d1 << " st1: " << st1 << " d2: " << d2 << " st2: " << st2 << std::endl;
        }
      }

      //std::cout << "mom: " << i << " d1: " << d1 << " d2: " << d2 << std::endl; 
      
      fastjet::PseudoJet pd1(pythia.event[d1].px(),pythia.event[d1].py(),pythia.event[d1].pz(),pythia.event[d1].e());
      pd1.set_user_info(new extraInfo(pythia.event[d1].id(), -2)); 
      partons.push_back(pd1);

      fastjet::PseudoJet pd2(pythia.event[d2].px(),pythia.event[d2].py(),pythia.event[d2].pz(),pythia.event[d2].e());
      pd2.set_user_info(new extraInfo(pythia.event[d2].id(), -2)); 
      partons.push_back(pd2);
      
    }

    if(i == 0)
    {
      // current date/time based on current system
      time_t now = time(0);
      // convert now to string form
      const char* dt = ctime(&now);
      string path = "../event_";
      path += dt;
      path += ".txt";
      
      FILE *fp = freopen(path.c_str(),"w",stdout);
      pythia.event.list();
      fclose(fp);
    }
  }

  //if(iprint==1) pythia.event.list();
  //pythia.event.list();


  return particles;
}


#endif
