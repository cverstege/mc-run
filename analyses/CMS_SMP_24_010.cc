// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {

  // Main analysis class for CMS_SMP_24_010
  class CMS_SMP_24_010 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_SMP_24_010);


    // === Analysis selection constants ===
    const double _muondressingDeltaR = 0.1;
    const double _minmuonptcut = 29*GeV;
    const double _maxmuoneta = 2.4;
    const double _zmass = 91.1876*GeV;
    const double _zmassdiff = 20*GeV;
    const double _minptZ = 25*GeV;
    const double _minjetpt = 20*GeV;
    const double _maxabsjetrap = 2.4;
    const double _jetMuCleaningDeltaR = 0.4;


    // === Rivet init: projections and histograms ===
    void init() {

      // Print analysis cuts
      MSG_INFO("Analysis cuts:\n"
        << "\tMin jet pT: " << _minjetpt << "\n"
        << "\tMax jet rapidity: " << _maxabsjetrap << "\n"
        << "\tMin muon pT: " << _minmuonptcut/GeV << " GeV\n"
        << "\tMax muon eta: " << _maxmuoneta << "\n"
        << "\tZ-mass window +/-: " << _zmassdiff/GeV << " GeV\n"
        << "\tMin Z pT: " << _minptZ/GeV << " GeV\n");

      const FinalState fs(Cuts::abseta < 4.9);
      // CMS GenJets include muons, cleaning is done later with deltaR with respect to Z-muons
      declare(FastJets(fs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::NONE), "jetsAK4");
      declare(FastJets(fs, JetAlg::ANTIKT, 0.8, JetMuons::ALL, JetInvisibles::NONE), "jetsAK8");

      // Dress Muons with photons and find Z boson
      const Cut muon_cuts = Cuts::abspid == PID::MUON
                            && Cuts::abseta < _maxmuoneta
                            && Cuts::pT > _minmuonptcut;
      const Cut z_cuts = Cuts::massIn(_zmass - _zmassdiff, _zmass + _zmassdiff)
                         && Cuts::pT > _minptZ;
      DileptonFinder zfinder_mm_dressed(
        fs,
        _zmass,
		    0.1,
		    muon_cuts,
		    z_cuts
      );
      declare(zfinder_mm_dressed, "zfinder_mm_dressed");

      // Book histograms for jet multiplicity and Z pT in y* and yboost bins
      vector<double> binedges_ZPt;
      for (const string _jettype : {"AK4", "AK8"}) {
        book(_h["NJets" + _jettype], "NJets" + _jettype, 10, 1, 11);
        // Bin edges and naming use lower bound
        for (double _ystar : _binedges_Ystar) {
          for (double _yboost : _binedges_Yboost) {
            if (_ystar + _yboost > 2.0) continue;
            if (_ystar >= 2.0) {
              // extra binning extreme ystar
              binedges_ZPt = {25., 30., 40., 50., 70., 90., 110., 150., 250.};
            } else if ((_ystar + _yboost >= 2.0) || (_yboost < 0.5 && _ystar >= 1.5)) {
              // edge binning
              binedges_ZPt = {25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 110., 130., 150., 170., 190., 250., 1000.};
            } else {
              // default binning
              binedges_ZPt = {25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 110., 130., 150., 170., 190., 220., 250., 400., 1000.};
            }
            const string _hist_ZPt_ident = "ZPt" + _jettype + "Ys" + to_string(_ystar) + "Yb" + to_string(_yboost);
            book(_h[_hist_ZPt_ident], _hist_ZPt_ident, binedges_ZPt);
          }
        }
      }
      MSG_INFO("Booked " << _h.size() << " histograms");
      MSG_DEBUG("Histograms:");
      for (const auto& h : _h) {
        MSG_DEBUG("\t" << h.first);
      }
    }


    // === Per-event analysis: selection, reconstruction, and histogramming ===
    void analyze(const Event& event) {

      // --- Z boson candidate selection ---
      const ZFinder& zfinder_mm_dressed = apply<ZFinder>(event, "zfinder_mm_dressed");
      if (zfinder_mm_dressed.bosons().empty()) {
        vetoEvent;
      }

      // Get the Z boson and its constituent muons from the finder
      const Particle& Z = zfinder_mm_dressed.boson();
      const Particles& z_muons = zfinder_mm_dressed.constituents();
      const double rap_Z = Z.rap();
      const double pT_Z = Z.pT();
      if (pT_Z <= _minptZ) {
        vetoEvent;
      }
      MSG_DEBUG("Z-boson mass " << Z.mass()/GeV << " GeV");

      // --- Jet selection and cleaning ---
      map<string, Jets> _jetcollections;
      _jetcollections["AK4"] = apply<FastJets>(event, "jetsAK4").jetsByPt(Cuts::absrap < _maxabsjetrap && Cuts::pT > _minjetpt);
      _jetcollections["AK8"] = apply<FastJets>(event, "jetsAK8").jetsByPt(Cuts::absrap < _maxabsjetrap && Cuts::pT > _minjetpt);

      // Clean jets
      bool _jetpass = false;
      for (auto& entry : _jetcollections) {
        idiscardIfAnyDeltaRLess(entry.second, z_muons, _jetMuCleaningDeltaR);
        MSG_DEBUG("#Jets" << entry.first << " = " << entry.second.size());
        if (!entry.second.empty()) {
          _jetpass = true;
        }
      }
      if (!_jetpass) {
        vetoEvent;
      }

      // --- Histogram filling ---
      for (const auto& entry : _jetcollections) {
        if (entry.second.empty()) {
          MSG_DEBUG("No jets in collection " << entry.first << ", skipping histogram filling.");
          continue;
        }
        _h["NJets" + entry.first]->fill(entry.second.size());
        const double rap_Jet1 = entry.second.at(0).rap();
        const double rap_star = 0.5 * abs(rap_Z - rap_Jet1);
        const double rap_boost = 0.5 * abs(rap_Z + rap_Jet1);
        const pair<double, double> bin = findYstarYboostBin(rap_star, rap_boost);
        MSG_DEBUG("y*: " << rap_star << ", y_b: " << rap_boost);
        if (bin.first != -1.0 && bin.second != -1.0) {
          const string _hist_ZPt_ident = "ZPt" + entry.first + "Ys" + to_string(bin.first) + "Yb" + to_string(bin.second);
          MSG_DEBUG("Filling histogram: " << _hist_ZPt_ident);
          _h[_hist_ZPt_ident]->fill(pT_Z);
        }
        else {
          MSG_DEBUG("No valid y* and yboost bin found for rapidity star: " << rap_star << ", rapidity boost: " << rap_boost);
        }
      }
    }


    // === Finalization: normalize histograms ===
    void finalize() {
      const double sf = crossSection() / picobarn / sumOfWeights();
      for (const auto& _hist : _h) {
        _hist.second->scaleW(sf);
      }
    }


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    /// @}


    // y* and yboost bin edges (public as per analysis style recommendation)
    vector<double> _binedges_Ystar = {0.0, 0.5, 1.0, 1.5, 2.0};
    vector<double> _binedges_Yboost = {0.0, 0.5, 1.0, 1.5, 2.0};

    // Helper: find the appropriate y* and yboost bin (public as per analysis style recommendation)
    pair<double, double> findYstarYboostBin(double rap_star, double rap_boost) const {
      // Bin edges and naming use lower bound
      for (double _ystar : _binedges_Ystar) {
        if (rap_star > _ystar + 0.5) continue;
        for (double _yboost : _binedges_Yboost) {
          if (rap_boost > _yboost + 0.5) continue;
          if (_ystar + _yboost > 2.0) break; // No valid bin if y* + y > 2.0
          return make_pair(_ystar, _yboost);
        }
      }
      return make_pair(-1.0, -1.0); // Not found
    }

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_SMP_24_010);

}
