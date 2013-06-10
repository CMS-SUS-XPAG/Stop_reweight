/////// ##### Stop -> top neutralino reweighting and utilitites #####
/////// J. Alcaraz, B. De La Cruz (CMS, Feb 2013)
#include "TROOT.h"
#include "TLorentzVector.h"

struct SUSYGenParticle { // To be filled with status-3 genParticles
      int pdgId; // PDG identifier (with sign, please)
      int firstMother; // first mother, set to <0 if no mothers
      double energy; // energy [GeV]
      double pt; // pt [GeV]
      double eta; // eta
      double phi; // phi
};

// Get theta effective mixing for a given top polarization hypothesis
// The result also depends on the stop, top and chi0 mass hypotheses
// Also valid for an off-shell scenario
double GetThetaMixing(double topPol, double m_stop, double m_top, double m_chi0) {
      double p_chi0 = sqrt(pow(m_top*m_top+m_chi0*m_chi0-m_stop*m_stop,2)/4 - pow(m_top*m_chi0,2))/m_stop;
      double e_chi0 = sqrt(p_chi0*p_chi0+m_chi0*m_chi0);
      double sqrPol = 0.; if (fabs(topPol)<1) sqrPol = sqrt(1-topPol*topPol);
      double tanThetaEff = (p_chi0*sqrPol-m_chi0*topPol)/(topPol*e_chi0+p_chi0);
      // This is also a valid solution leading to the same top polarization value
      //double tanThetaEff = (-p_chi0*sqrPol-m_chi0*topPol)/(topPol*e_chi0+p_chi0);
      return atan(tanThetaEff);
}

// Get the top polarization given an effective theta mixing hypothesis
// The result also depends on the stop, top and chi0 mass hypotheses
// Also valid for an off-shell scenario
double GetTopPolarization(double thetaMixing, double m_stop, double m_top, double m_chi0) {
      double p_chi0 = sqrt(pow(m_top*m_top+m_chi0*m_chi0-m_stop*m_stop,2)/4 - pow(m_top*m_chi0,2))/m_stop;
      double e_chi0 = sqrt(p_chi0*p_chi0+m_chi0*m_chi0);
      double pol = p_chi0*cos(2*thetaMixing)/(e_chi0+m_chi0*sin(2*thetaMixing));
      return pol;
}

// This function estimates the average polarization in the event
// for a given thetaMixingEffective scenario. We note that
// top and (-)antitop polarizations only coincide for narrow stop, top, chi0 masses
// I.e., the function only makes full sense for an on-shell case with narrow widths
double AverageTopPolarization_Stop_to_TopChi0(double thetaMixing, std::vector<SUSYGenParticle> genParticles) {

    double pol1 = 0.;
    double pol2 = 0.;

    unsigned int ngen = genParticles.size();

    for (unsigned int ig=0; ig<ngen; ++ig) {
      // Find the stop
      const SUSYGenParticle& genStop = genParticles[ig];
      if (abs(genStop.pdgId)!=1000006) continue; // find the stop
      int igstop = ig;

      TLorentzVector stop4;
      stop4.SetPtEtaPhiE(genStop.pt, genStop.eta, genStop.phi, genStop.energy);
      TVector3 betaV(-stop4.Px()/stop4.Energy(),-stop4.Py()/stop4.Energy(),-stop4.Pz()/stop4.Energy());

      for (unsigned int ig=igstop+1; ig<ngen; ++ig) { // assume chi0 comes later
            // Move chi0 to the stop rest frame
            const SUSYGenParticle& genChi0 = genParticles[ig];
            if (genChi0.firstMother!=igstop) continue;
            if (abs(genChi0.pdgId)!=1000022) continue;

            TLorentzVector chi04;
            chi04.SetPtEtaPhiE(genChi0.pt, genChi0.eta, genChi0.phi, genChi0.energy);
            chi04.Boost(betaV);

            if (genParticles[igstop].pdgId>0) {
                  pol1 = chi04.P()*cos(2*thetaMixing)/(chi04.Energy()+chi04.M()*sin(2*thetaMixing));
            } else {
                  pol2 = chi04.P()*cos(2*thetaMixing)/(chi04.Energy()+chi04.M()*sin(2*thetaMixing));
            }
      
            break;
      }
      
    }

    return (pol1+pol2)/2;
}

// The following reweighting only makes sense for on-shell stop, top and chi0
// In the off-shell case top and anti-top may get very different polarizations
double Reweight_Stop_to_TopChi0 (std::vector<SUSYGenParticle> genParticles, double referenceTopPolarization, double requestedTopPolarization) {
    double weight = 1.;

    unsigned int ngen = genParticles.size();

    for (unsigned int ig=0; ig<ngen; ++ig) {
      const SUSYGenParticle& gen = genParticles[ig];
      if (gen.firstMother<0) continue;
      if (abs(gen.pdgId)>20) continue; // expect quarks or leptons from W decay

      // Navigate upwards in the stop->top->W->fermion decay chain
      const SUSYGenParticle& genW = genParticles[gen.firstMother];
      if (genW.firstMother<0) continue;
      if (abs(genW.pdgId)!=24) continue;
      const SUSYGenParticle& genTop = genParticles[genW.firstMother];
      if (abs(genTop.pdgId)!=6) continue;

      // We only care about the down-type fermion
      if (genTop.pdgId*gen.pdgId>0) continue;

      // We also need a stop
      if (genTop.firstMother<0) continue;
      const SUSYGenParticle& genStop = genParticles[genTop.firstMother];
      if (abs(genStop.pdgId)!=1000006) continue;

      // Move top and fermion to the stop center-of-mass frame
      TLorentzVector stop4;
      stop4.SetPtEtaPhiE(genStop.pt, genStop.eta, genStop.phi, genStop.energy);
      TVector3 betaV(-stop4.Px()/stop4.Energy(),-stop4.Py()/stop4.Energy(),-stop4.Pz()/stop4.Energy());

      TLorentzVector top4;
      top4.SetPtEtaPhiE(genTop.pt, genTop.eta, genTop.phi, genTop.energy);
      top4.Boost(betaV);

      TLorentzVector ferm4;
      ferm4.SetPtEtaPhiE(gen.pt, gen.eta, gen.phi, gen.energy);
      ferm4.Boost(betaV);

      // Do not reweight if by any reason top/fermion directions are undefined
      // This should be pathological if things are fine
      if (top4.P()<=0 || ferm4.P()<=0) {
            printf("Warning: particles at rest, no weight applied: ptop: %.3e, pf: %.3e\n", top4.P(), ferm4.P());
            continue; 
      }

      double costh = (top4.Px()*ferm4.Px()+top4.Py()*ferm4.Py()+top4.Pz()*ferm4.Pz())/top4.P()/ferm4.P();
      
      double weight_L = (top4.Energy()+top4.P())*(1-costh);
      double weight_R = (top4.Energy()-top4.P())*(1+costh);
      weight *= ((1+requestedTopPolarization)*weight_R+(1-requestedTopPolarization)*weight_L)/((1+referenceTopPolarization)*weight_R+(1-referenceTopPolarization)*weight_L);

    }

    return weight;

};

// The following reweighting makes sense in iall cases (off-shell case too)
// It assumes isotropy for the reference MC and any SUSY scenario as target 
// (via a thetaMixingTarget input)
double Reweight_Stop_to_TopChi0_with_SUSYmodel (std::vector<SUSYGenParticle> genParticles, double thetaMixingTarget) { 
    double weight = 1.;

    unsigned int ngen = genParticles.size();

    for (unsigned int ig=0; ig<ngen; ++ig) {
      const SUSYGenParticle& gen = genParticles[ig];
      if (gen.firstMother<0) continue;
      if (abs(gen.pdgId)>20) continue; // expect quarks or leptons from W decay

      // Navigate upwards in the stop->top->W->fermion decay chain
      const SUSYGenParticle& genW = genParticles[gen.firstMother];
      if (genW.firstMother<0) continue;
      if (abs(genW.pdgId)!=24) continue;
      const SUSYGenParticle& genTop = genParticles[genW.firstMother];
      if (abs(genTop.pdgId)!=6) continue;

      // We only care about the down-type fermion
      if (genTop.pdgId*gen.pdgId>0) continue;

      // We also need a stop
      if (genTop.firstMother<0) continue;
      const SUSYGenParticle& genStop = genParticles[genTop.firstMother];
      if (abs(genStop.pdgId)!=1000006) continue;

      // Move top and fermion to the stop center-of-mass frame
      TLorentzVector stop4;
      stop4.SetPtEtaPhiE(genStop.pt, genStop.eta, genStop.phi, genStop.energy);
      TVector3 betaV(-stop4.Px()/stop4.Energy(),-stop4.Py()/stop4.Energy(),-stop4.Pz()/stop4.Energy());

      int igstop = genTop.firstMother;
      int igchi0 = -1;
      TLorentzVector chi04;
      for (unsigned int ig=igstop+1; ig<ngen; ++ig) { // assume chi0 comes later
            // Move chi0 to the stop rest frame
            const SUSYGenParticle& genChi0 = genParticles[ig];
            if (genChi0.firstMother!=igstop) continue;
            if (abs(genChi0.pdgId)!=1000022) continue;
            igchi0 = ig;
            chi04.SetPtEtaPhiE(genChi0.pt, genChi0.eta, genChi0.phi, genChi0.energy);
            chi04.Boost(betaV);
            break;
      }
      if (igchi0<0) continue;

      // Determine top polarization in this decay according to model
      double topPolarization = chi04.P()*cos(2*thetaMixingTarget)/(chi04.Energy()+chi04.M()*sin(2*thetaMixingTarget));

      /*
      if (genTop.pdgId>0) {
            printf("Top polarization: %.3f\n", topPolarization);
      } else {
            printf("Antitop polarization: %.3f\n", -topPolarization);
      }
      */
      
      TLorentzVector top4;
      top4.SetPtEtaPhiE(genTop.pt, genTop.eta, genTop.phi, genTop.energy);
      top4.Boost(betaV);

      TLorentzVector ferm4;
      ferm4.SetPtEtaPhiE(gen.pt, gen.eta, gen.phi, gen.energy);
      ferm4.Boost(betaV);

      // Do not reweight if by any reason top/fermion directions are undefined
      // This should be pathological if things are fine
      if (top4.P()<=0 || ferm4.P()<=0) {
            printf("Warning: particles at rest, no weight applied: ptop: %.3e, pf: %.3e\n", top4.P(), ferm4.P());
            continue; 
      }

      double costh = (top4.Px()*ferm4.Px()+top4.Py()*ferm4.Py()+top4.Pz()*ferm4.Pz())/top4.P()/ferm4.P();
      
      double weight_L = (top4.Energy()+top4.P())*(1-costh);
      double weight_R = (top4.Energy()-top4.P())*(1+costh);
      weight *= ((1+topPolarization)*weight_R+(1-topPolarization)*weight_L)/(weight_R+weight_L);

    }

    return weight;

};

// The following reweighting also makes sense in the off-shell case
// It goes from a SUSY MC generated with thetaMixingReference to another
// SUSY MC scenario with thetaMixingTarget
double Reweight_Stop_to_TopChi0_with_SUSYmodel (std::vector<SUSYGenParticle> genParticles, double thetaMixingReference, double thetaMixingTarget) { 
    double weight = 1.;

    unsigned int ngen = genParticles.size();

    for (unsigned int ig=0; ig<ngen; ++ig) {
      const SUSYGenParticle& gen = genParticles[ig];
      if (gen.firstMother<0) continue;
      if (abs(gen.pdgId)>20) continue; // expect quarks or leptons from W decay

      // Navigate upwards in the stop->top->W->fermion decay chain
      const SUSYGenParticle& genW = genParticles[gen.firstMother];
      if (genW.firstMother<0) continue;
      if (abs(genW.pdgId)!=24) continue;
      const SUSYGenParticle& genTop = genParticles[genW.firstMother];
      if (abs(genTop.pdgId)!=6) continue;

      // We only care about the down-type fermion
      if (genTop.pdgId*gen.pdgId>0) continue;

      // We also need a stop
      if (genTop.firstMother<0) continue;
      const SUSYGenParticle& genStop = genParticles[genTop.firstMother];
      if (abs(genStop.pdgId)!=1000006) continue;

      // Move top and fermion to the stop center-of-mass frame
      TLorentzVector stop4;
      stop4.SetPtEtaPhiE(genStop.pt, genStop.eta, genStop.phi, genStop.energy);
      TVector3 betaV(-stop4.Px()/stop4.Energy(),-stop4.Py()/stop4.Energy(),-stop4.Pz()/stop4.Energy());

      int igstop = genTop.firstMother;
      int igchi0 = -1;
      TLorentzVector chi04;
      for (unsigned int ig=igstop+1; ig<ngen; ++ig) { // assume chi0 comes later
            // Move chi0 to the stop rest frame
            const SUSYGenParticle& genChi0 = genParticles[ig];
            if (genChi0.firstMother!=igstop) continue;
            if (abs(genChi0.pdgId)!=1000022) continue;
            igchi0 = ig;
            chi04.SetPtEtaPhiE(genChi0.pt, genChi0.eta, genChi0.phi, genChi0.energy);
            chi04.Boost(betaV);
            break;
      }
      if (igchi0<0) continue;

      // Determine top polarization in this decay according to models
      double topPolarizationReference = chi04.P()*cos(2*thetaMixingReference)/(chi04.Energy()+chi04.M()*sin(2*thetaMixingReference));
      double topPolarizationTarget = chi04.P()*cos(2*thetaMixingTarget)/(chi04.Energy()+chi04.M()*sin(2*thetaMixingTarget));

      /*
      if (genTop.pdgId>0) {
            printf("Top polarization: %.3f\n", topPolarization);
      } else {
            printf("Antitop polarization: %.3f\n", -topPolarization);
      }
      */
      
      TLorentzVector top4;
      top4.SetPtEtaPhiE(genTop.pt, genTop.eta, genTop.phi, genTop.energy);
      top4.Boost(betaV);

      TLorentzVector ferm4;
      ferm4.SetPtEtaPhiE(gen.pt, gen.eta, gen.phi, gen.energy);
      ferm4.Boost(betaV);

      // Do not reweight if by any reason top/fermion directions are undefined
      // This should be pathological if things are fine
      if (top4.P()<=0 || ferm4.P()<=0) {
            printf("Warning: particles at rest, no weight applied: ptop: %.3e, pf: %.3e\n", top4.P(), ferm4.P());
            continue; 
      }

      double costh = (top4.Px()*ferm4.Px()+top4.Py()*ferm4.Py()+top4.Pz()*ferm4.Pz())/top4.P()/ferm4.P();
      
      double weight_L = (top4.Energy()+top4.P())*(1-costh);
      double weight_R = (top4.Energy()-top4.P())*(1+costh);
      weight *= ((1+topPolarizationTarget)*weight_R+(1-topPolarizationTarget)*weight_L)/((1+topPolarizationReference)*weight_R+(1-topPolarizationReference)*weight_L);

    }

    return weight;

};
