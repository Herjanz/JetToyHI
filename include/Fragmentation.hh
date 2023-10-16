#ifndef __Fragmentation_HH__
#define __Fragmentation_HH__

//------------------------------------------------------------------------
/// Fragmentation Function
// pT/Z probability graph

class Fragmentation {
public:
  /// default ctor
  Fragmentation(){}

  /// compute the function
  virtual vector<double> getFF(const fastjet::PseudoJet &jet) const {
        // check the jet is appropriate for computation
    if (!jet.has_constituents()) {
      Printf("Fragmentation can only be applied on jets for which the constituents are known.");
      return vector<double>();
    }
    vector<fastjet::PseudoJet> constits = jet.constituents();
    double Z = jet.pt();
    vector<double> PtZ; PtZ.reserve(constits.size());
    for(fastjet::PseudoJet p : constits) {
      PtZ.push_back(p.pt()/Z);
    }
    return PtZ;
  }
};

#endif
