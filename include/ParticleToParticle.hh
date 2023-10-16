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

class ParticleToParticle {
public:
  /// default ctor
  ParticleToParticle(){}

  /// compute the function
  virtual vector<PtoPInfo> findMatches(vector<fastjet::PseudoJet> incoming, vector<fastjet::PseudoJet> outgoing, bool childParent = false) const {
    vector<PtoPInfo> matches;

    for(PseudoJet outP : outgoing)
    {
      PtoPInfo info;
      info.out = outP;

      for (PseudoJet inP : incoming)
      {
        if(outP == inP)
          continue;


        if(childParent && outP.pt() / inP.pt() > 1)
          continue;
          
        // temporary
        double dr = inP.delta_R(outP);
        if(info.dr == -1 || dr < info.dr)
        {

          info.dr = dr;
          info.ptFraction = outP.pt() / inP.pt();
          info.in = inP;
        }
      }

      matches.push_back(info);
    }

    return matches;
  }
};

#endif
