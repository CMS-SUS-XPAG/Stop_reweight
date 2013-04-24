/////// ##### Stop -> b chargino reweighting #####
/////// J. Alcaraz, B. De La Cruz (CMS, Apr 2013)
/////// Mostly built from expressions in: 
///////     a) I. Low, arXiv:1304.0491
///////     b) J.F. Gunion, H.E. Haber, Phys. Rev. D37 (1988) 2515 (+erratum)

// This code assumes:
//    a) That the reference Pythia decay for the chargino is performed in two steps:
//       1) chi+ -> chi0 W*, 2) W* -> f fbar', and that both decays are done 
//       isotropically in the center-of-mass (fully flat phase space). This is 
//       never the case in reality, due to the presence of non-zero spins and the 
//       V-A coupling of the W to SM fermions.
//    b) On-shell W production. However, we think that the derived weight is (at least)
//       a good approximation for the off-shell case too (to be studied in more detail).
//       
// The weight depends on the left-right mixing in the chargino-stop-b coupling
// and on the left-right mixing in the chargino-neutralino-W coupling. Assuming no
// pathological cancellations in SUSY mixing matrices, the low tanbeta limit,
// tanbeta << mtop/mb leads to left-handed charginos (thataChi_eff=0), while in 
// the huge tanbeta limit, tanbeta >> mtop/mb (thetaChi_eff=+-pi/2), the chargino 
// is right-handed.

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

double Norm(double M1, double M2, double MV, double CL, double CR) {
      double lambda = pow(M1,4) + pow(M2,4) + pow(MV,4) - 2*pow(M1*M2,2) - 2*pow(M1*MV,2) - 2*pow(M2*MV,2);
      double norm = (CL*CL+CR*CR)*(lambda + 3*MV*MV*(M1*M1+M2*M2-MV*MV)) - 12*CL*CR*M1*M2*MV*MV;
      norm /= 3.;
      return norm;
}

void Boost_To_Stop_Rest_Frame(TLorentzVector& stop4, TLorentzVector& chargino4, TLorentzVector& b4, TLorentzVector& neutralino4, TLorentzVector& W4, TLorentzVector& up4, TLorentzVector& down4, TLorentzVector& s4){
      TVector3 betaV(-stop4.Px()/stop4.Energy(),-stop4.Py()/stop4.Energy(),-stop4.Pz()/stop4.Energy());
      stop4.Boost(betaV);
      chargino4.Boost(betaV);
      b4.Boost(betaV);
      neutralino4.Boost(betaV);
      W4.Boost(betaV);
      up4.Boost(betaV);
      down4.Boost(betaV);
      s4.SetE(chargino4.P()/chargino4.M());
      s4.SetVect(chargino4.Vect().Unit()*chargino4.Gamma());
}

double Reweight_T2bW (double thetaChi_eff, double thetaW_eff, std::vector<SUSYGenParticle> genParticles) {
    double weight = 1.;

    unsigned int ngen = genParticles.size();

    for (unsigned int i_stop=0; i_stop<ngen; ++i_stop) {
      // Look for stops
      const SUSYGenParticle& gen = genParticles[i_stop];
      if (abs(gen.pdgId)!=1000006) continue;

      // Look for stop decay products
      int i_b = -1;
      int i_chargino = -1;
      for (unsigned int ig=i_stop+1; ig<ngen; ++ig) {
            const SUSYGenParticle& gen = genParticles[ig];
            if (abs(gen.firstMother)!=i_stop) continue;
            if (abs(gen.pdgId)==5) i_b = ig;
            else if (abs(gen.pdgId)==1000024) i_chargino = ig;
            if (i_b>=0 && i_chargino>=0) break;
      }
      if (i_b<0 || i_chargino<0) continue;

      int i_neutralino = -1;
      int i_W = -1;
      for (unsigned int ig=i_chargino+1; ig<ngen; ++ig) {
            const SUSYGenParticle& gen = genParticles[ig];
            if (abs(gen.firstMother)!=i_chargino) continue;
            if (abs(gen.pdgId)==24) i_W = ig;
            else if (abs(gen.pdgId)==1000022) i_neutralino = ig;
            if (i_W>=0 && i_neutralino>=0) break;
      }
      if (i_W<0 || i_neutralino<0) continue;

      int i_up = -1;
      int i_down = -1;
      for (unsigned int ig=i_W+1; ig<ngen; ++ig) {
            const SUSYGenParticle& gen = genParticles[ig];
            if (abs(gen.firstMother)!=i_W) continue;
            if (abs(gen.pdgId)%2==0) i_up = ig;
            else if (abs(gen.pdgId)%2==1) i_down = ig;
            if (i_up>=0 && i_down>=0) break;
      }
      if (i_up<0 || i_down<0) continue;

      const SUSYGenParticle& gen_stop = genParticles[i_stop];
      const SUSYGenParticle& gen_b = genParticles[i_b];
      const SUSYGenParticle& gen_chargino = genParticles[i_chargino];
      const SUSYGenParticle& gen_W = genParticles[i_W];
      const SUSYGenParticle& gen_neutralino = genParticles[i_neutralino];
      const SUSYGenParticle& gen_up = genParticles[i_up];
      const SUSYGenParticle& gen_down = genParticles[i_down];

      // Fill Lorentz four-vectors
      TLorentzVector stop4, chargino4, b4, neutralino4, W4, up4, down4;

      stop4.SetPtEtaPhiE(gen_stop.pt, gen_stop.eta, gen_stop.phi, gen_stop.energy);
      chargino4.SetPtEtaPhiE(gen_chargino.pt, gen_chargino.eta, gen_chargino.phi, gen_chargino.energy);
      b4.SetPtEtaPhiE(gen_b.pt, gen_b.eta, gen_b.phi, gen_b.energy);
      neutralino4.SetPtEtaPhiE(gen_neutralino.pt, gen_neutralino.eta, gen_neutralino.phi, gen_neutralino.energy);
      W4.SetPtEtaPhiE(gen_W.pt, gen_W.eta, gen_W.phi, gen_W.energy);
      up4.SetPtEtaPhiE(gen_up.pt, gen_up.eta, gen_up.phi, gen_up.energy);
      down4.SetPtEtaPhiE(gen_down.pt, gen_down.eta, gen_down.phi, gen_down.energy);

      // Reference spin four-vector along the chargino direction (filled in after boost)
      TLorentzVector s4;

      // Move everything to the stop center-of-mass frame
      Boost_To_Stop_Rest_Frame(stop4, chargino4, b4, neutralino4, W4, up4, down4, s4);

      double c_L = sin(thetaW_eff);
      double c_R = cos(thetaW_eff);
      double norm_target = Norm(chargino4.M(), neutralino4.M(), W4.M(), c_L, c_R);
      double target = 0;
      for (int hel = -1; hel<2; hel += 2) {
            TLorentzVector t4 = s4*hel;
            TLorentzVector chargino4_plus = chargino4 + t4*chargino4.M();
            TLorentzVector chargino4_minus = chargino4 - t4*chargino4.M();
            target += (1. - chargino4.M()*cos(2*thetaChi_eff)*(b4*t4)/((b4*chargino4) - (b4.M()*chargino4.M())*sin(2*thetaChi_eff)))/2 *
           (8*c_L*c_L*(down4*chargino4_minus)*(up4*neutralino4)
                      + 8*c_R*c_R*(up4*chargino4_plus)*(down4*neutralino4)
                      - 4*c_L*c_R*neutralino4.M()*(pow(W4.M(),2)*chargino4.M()-2*(down4*t4)*(up4*chargino4)+2*(down4*chargino4)*(up4*t4)))/norm_target; 
      }

      weight *= target;

    }

    return weight;

};
