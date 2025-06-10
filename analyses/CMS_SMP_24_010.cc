// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"

namespace Rivet {

  // Main analysis class for CMS_SMP_24_010
  class CMS_SMP_24_010 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_SMP_24_010);

    // === Analysis selection constants ===
    const double _muondressingDeltaR = 0.1;
    const double _muonisolationCut = 0.15;
    const double _minmuonptcut = 29*GeV;
    const double _maxmuoneta = 2.4;
    const double _zmassdiff = 20*GeV;
    const double _minptZ = 25*GeV;
    const double _minjetpt = 20*GeV;
    const double _maxabsjetrap = 2.4;
    const double _muonCleaningDeltaR = 0.4;

    // === Rivet init: projections and histograms ===
    void init() {
      // Print analysis cuts
      MSG_INFO("Analysis cuts:\n"
        << "\tMin jet pT: " << _minjetpt << "\n"
        << "\tMax jet rapidity: " << _maxabsjetrap << "\n"
        << "\tRelative muon isolation: " << _muonisolationCut << "\n"
        << "\tMin muon pT: " << _minmuonptcut/GeV << " GeV\n"
        << "\tMax muon eta: " << _maxmuoneta << "\n"
        << "\tZ-mass window +/-: " << _zmassdiff/GeV << " GeV\n"
        << "\tMin Z pT: " << _minptZ/GeV << " GeV\n");

      // Projections: final state, jets, prompt photons/muons, dressed muons
      declare(FinalState(Cuts::abseta < 4.9), "FS");
      // CMS GenJets include muons, cleaning is done later with deltaR with respect to Z-muons
      declare(FastJets(FinalState(Cuts::abseta < 4.9), FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::NONE), "jetsAK4");
      declare(FastJets(FinalState(Cuts::abseta < 4.9), FastJets::ANTIKT, 0.8, JetAlg::Muons::ALL, JetAlg::Invisibles::NONE), "jetsAK8");
      declare(PromptFinalState(Cuts::abspid == PID::PHOTON), "photons");
      declare(PromptFinalState(Cuts::abspid == PID::MUON), "bare_muons");
      // Dress bare muons with photons and apply selections
      Cut muon_cuts = Cuts::abseta < _maxmuoneta && Cuts::pT > _minmuonptcut;
      declare(DressedLeptons(getProjection<PromptFinalState>("photons"), getProjection<PromptFinalState>("bare_muons"), _muondressingDeltaR, muon_cuts), "muons");

      // Book histograms for jet multiplicity and Z pT in y* and yboost bins
      vector<double> binedges_ZPt;
      for (const string _jettype : {"AK4", "AK8"}) {
        book(_h["NJets" + _jettype], "NJets" + _jettype, 10, 1, 11);
        // Bin edges and naming use lower bound
        for (double _ystar : binedges_Ystar) {
          for (double _yboost : binedges_Yboost) {
            if (_ystar + _yboost > 2.0) continue;
            if (_ystar >= 2.0) {
              // extra binning extreme ystar
              binedges_ZPt = {25., 30., 40., 50., 70., 90., 110., 150., 250.};
            } else if ((_ystar + _yboost >= 2.0) || (_yboost < 0.5 && _ystar>=1.5)) {
              // edge binning
              binedges_ZPt = {25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 110., 130., 150., 170., 190., 250., 1000.};
            } else {
              // default binning
              binedges_ZPt = {25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 110., 130., 150., 170., 190., 220., 250., 400., 1000.};
            }
            string _hist_ZPt_ident = "ZPt" + _jettype + "Ys" + to_string(_ystar) + "Yb" + to_string(_yboost);
            book(_h[_hist_ZPt_ident], _hist_ZPt_ident, binedges_ZPt);
          }
        }
      }
      MSG_INFO("Booked " << _h.size() << " histograms");
      MSG_DEBUG("Histograms:");
      for (const auto& h : _h) MSG_DEBUG("\t" << h.first);
    }

    // === Per-event analysis: selection, reconstruction, and histogramming ===
    void analyze(const Event& event) {
      // --- Isolated dressed muon selection ---
      const auto& muons = apply<DressedLeptons>(event, "muons").dressedLeptons();
      const auto& fs_particles = apply<FinalState>(event, "FS").particles();
      vector<DressedLepton> isolated_muons;
      for (const auto& mu : muons) {
        double sumpt = 0.0;
        for (const auto& p : fs_particles) {
          if (p.genParticle() == mu.constituentLepton().genParticle()) continue;
          if (deltaR(mu, p) < 0.4) sumpt += p.pT();
        }
        if (sumpt / mu.pT() < _muonisolationCut) isolated_muons.push_back(mu);
      }
      MSG_DEBUG(muons.size() << " dressed muons; " << isolated_muons.size() << " isolated muons");
      if (isolated_muons.size() < 2) vetoEvent;

      // --- Z boson candidate selection ---
      bool bosoncandidateexists = false;
      double best_massdiff = _zmassdiff;
      const DressedLepton* muon = nullptr;
      const DressedLepton* antimuon = nullptr;
      for (size_t i = 0; i < isolated_muons.size(); ++i) {
        for (size_t j = i + 1; j < isolated_muons.size(); ++j) {
          if (isolated_muons[i].pid() != -isolated_muons[j].pid()) continue;
          double mass = (isolated_muons[i].mom() + isolated_muons[j].mom()).mass();
          double massdiff = fabs(mass - 91.1876*GeV);
          if (massdiff < best_massdiff) {
            bosoncandidateexists = true;
            best_massdiff = massdiff;
            if (isolated_muons[i].pid() > 0) {
              muon = &isolated_muons[i];
              antimuon = &isolated_muons[j];
            } else {
              muon = &isolated_muons[j];
              antimuon = &isolated_muons[i];
            }
          }
        }
      }
      if (!bosoncandidateexists) vetoEvent;
      const double rap_Z = (muon->mom() + antimuon->mom()).rap();
      const double pT_Z = (muon->mom() + antimuon->mom()).pT()/GeV;
      if (pT_Z <= _minptZ) vetoEvent;
      auto z_muons = make_pair(*muon, *antimuon);
      MSG_DEBUG("Z-boson mass " << (muon->mom() + antimuon->mom()).mass()/GeV << " GeV");

      // --- Jet selection and cleaning ---
      map<string,Jets> _jetcollections;
      _jetcollections["AK4"] = apply<FastJets>(event, "jetsAK4").jetsByPt(Cuts::absrap < _maxabsjetrap && Cuts::pT > _minjetpt);
      _jetcollections["AK8"] = apply<FastJets>(event, "jetsAK8").jetsByPt(Cuts::absrap < _maxabsjetrap && Cuts::pT > _minjetpt);
      MSG_DEBUG(_jetcollections["AK4"].size() << " AK4 jets");
      MSG_DEBUG(_jetcollections["AK8"].size() << " AK8 jets");
      set<string> _jetcollectionstoerase;
      bool jet1pass = false;
      // Prepare a vector of DressedLepton for cleaning
      vector<DressedLepton> z_leptons = {*muon, *antimuon};
      for (auto& jets: _jetcollections) {
        idiscardIfAnyDeltaRLess(jets.second, z_leptons, _muonCleaningDeltaR);
        MSG_DEBUG(jets.second.size() << " cleaned " << jets.first << " jets");
        // Require at least one hard jet
        if (!jets.second.empty()) {
          if (jets.second.at(0).pT() > _minjetpt) {
            MSG_DEBUG("Hardest " << jets.first << " jet pt: " << jets.second.at(0).pT());
            jet1pass = true;
          } else {
            _jetcollectionstoerase.insert(jets.first);
          }
        }
        else {
            _jetcollectionstoerase.insert(jets.first);
        }
      }
      if (!(jet1pass)) vetoEvent;
      for (string c: _jetcollectionstoerase) {
        _jetcollections.erase(c);
      }

      // --- Histogram filling ---
      for (const auto& jets: _jetcollections) {
        MSG_DEBUG("Filling Histograms for " << jets.first);
        _h["NJets"+jets.first]->fill(jets.second.size());
        double rap_Jet1 = jets.second.at(0).rap();
        double rap_star = 0.5 * abs(rap_Z - rap_Jet1);
        double rap_boost = 0.5 * abs(rap_Z + rap_Jet1);
        auto bin = findYstarYboostBin(rap_star, rap_boost);
        if (bin.first != -1.0 && bin.second != -1.0) {
          string _hist_ZPt_ident = "ZPt"+jets.first+"Ys"+to_string(bin.first)+"Yb"+to_string(bin.second);
          _h[_hist_ZPt_ident]->fill(pT_Z);
        }
      }
    }

    // === Finalization: normalize histograms ===
    void finalize() {
      const double sf = crossSection()/picobarn/sumOfWeights();
      for(auto const& _hist : _h) scale(_hist.second, sf);
    }

    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    /// @}

  private:
    // y* and yboost bin edges
    std::vector<double> binedges_Ystar = {0.0, 0.5, 1.0, 1.5, 2.0};
    std::vector<double> binedges_Yboost = {0.0, 0.5, 1.0, 1.5, 2.0};
    // Helper: find the appropriate y* and yboost bin
    pair<double, double> findYstarYboostBin(double rap_star, double rap_boost) const {
      // Bin edges and naming use lower bound
      for (double _ystar : binedges_Ystar) {
        if (rap_star < _ystar+0.5) continue;
        for (double _yboost : binedges_Yboost) {
          if (_ystar + _yboost > 2.) continue;
          if (rap_boost < _yboost+0.5) {
            double _ystar_label = _ystar;
            double _yboost_label = _yboost;
            return make_pair(_ystar_label, _yboost_label);
          }
        }
      }
      return make_pair(-1.0, -1.0); // Not found
    }
  };

  RIVET_DECLARE_PLUGIN(CMS_SMP_24_010);

}
