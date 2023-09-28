#ifndef __Fragmentation_HH__
#define __Fragmentation_HH__

//------------------------------------------------------------------------
/// Fragmentation Function
// pT/Z probability graph

class Fragmentation {
public:
  /// default ctor
  Fragmentation(double beta=1.0, double kappa = 1., double R0 = 0.4) :
    _beta(beta),
    _kappa(kappa),
    _R0(R0)
  {}

  virtual vector<double> test(const fastjet::PseudoJet &jet) const {
    vector<double> v;
    v.push_back(3.);
    return v;
  }

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
  
protected:
  double _beta;
  double _kappa;
  double _R0;
};

#endif
