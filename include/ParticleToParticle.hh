#ifndef __ParticleToParticle_HH__
#define __ParticleToParticle_HH__

//------------------------------------------------------------------------
/// ParticleToParticle Function
// calculates distances from outgoing particles/jets to nearest incoming particles

struct PtoPInfo {
  fastjet::PseudoJet in;
  fastjet::PseudoJet out;
  double dr = -1;
  double ptFraction = -1;
};

vector<double> getDrVector(vector<PtoPInfo> matches)
{
  vector<double> dr;
  for(PtoPInfo i : matches)
  {
    dr.push_back(i.dr);
  }
  return dr;
}

vector<double> getPtFracVector(vector<PtoPInfo> matches)
{
  vector<double> ptFrac;
  for(PtoPInfo i : matches)
  {
    ptFrac.push_back(i.ptFraction);
  }
  return ptFrac;
}

int getPDG(fastjet::PseudoJet p)
{
  return p.user_info<PU14>().pdg_id(); 
}

void printInfo(PtoPInfo i)
{
  const int & PDGin = getPDG(i.in);
  //const int & PDGout = i.out.user_info<PU14>().pdg_id();
  cout << "dr: " << i.dr << " frac: " << i.ptFraction << endl;
  cout << "in phi " << i.in.phi() << " eta " << i.in.eta() << " pt " << i.in.pt() << " PDG " << PDGin << endl;
  cout << "out phi " << i.out.phi() << " eta " << i.out.eta() << " pt " << i.out.pt() << endl;
}

class ParticleToParticle {
public:
  /// default ctor
  ParticleToParticle(){}

  /// compute the function
  virtual vector<PtoPInfo> findMatches(vector<fastjet::PseudoJet> incoming, vector<fastjet::PseudoJet> outgoing, bool childParent = false) const {
    vector<PtoPInfo> matches;

    for(PseudoJet outP : outgoing)
    {
      PtoPInfo info, second;
      info.out = outP;
      second.out = outP;

      for (PseudoJet inP : incoming)
      {
        if(outP == inP)
          continue;

        // temporary
        double dr = inP.delta_R(outP);
        if(info.dr == -1 || dr < info.dr)
        {
          second.dr = info.dr;
          second.ptFraction = info.ptFraction;
          second.in = info.in;

          info.dr = dr;
          info.ptFraction = outP.pt() / inP.pt();
          info.in = inP;
        }
        else if(second.dr == -1 || dr < second.dr)
        {
          second.dr = dr;
          second.ptFraction = outP.pt() / inP.pt();
          second.in = inP;
        }
      }

      if(childParent)
      {
        if(info.ptFraction > 1.1 && second.ptFraction < info.ptFraction)
          if(abs(second.dr - info.dr) < 0.1)
            info = second;
      }


      matches.push_back(info);
    }

    return matches;
  }
};

#endif
