#include "Stop_TopChi0_Reweighting.C" 

#include "TStyle.h"
#include "TRint.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TVector3.h"
#include "TF1.h"
#include "MuTree.h"

bool IS_BATCH = false;
double POL = -1.; // -1 ==> fully left-handed top; +1 ==> fully right-handed top
int maxEventsUsed = -1; // If <0 ==> use all events in the file
//int maxEventsUsed = 100000; // If <0 ==> use all events in the file

using namespace std;
using namespace ciemat;

int main(int argc, char** argv){

  // Set it to "true" if you do not want to see the histograms interactively
  gROOT->SetBatch(IS_BATCH);
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(1111);

  // Open the input root file and set branches
  // On pcciet3a, pcciet3b, pccmscie6
  //TString sampleFile = "/data4/Fall11_WplusC_Trees_July2012/TTbar.root";
  // Just for some checks (there are selection cuts applied on the following file)
  //TString sampleFile = "/data4/Fall11_WplusC_Trees_July2012/TTbar_SSVHPNOMTNOISOreduced.root";
  // Stop file
//  TString sampleFile = "/data4/Fall11_WplusC_Trees_July2012/Stop.root";
  TString sampleFile = "dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/store/user/delacruz/STop2012/NTuplesFeb2013/v1/merge_stops_signalmc_T2tt_Mstop-225to1200_mLSP-0to1000_Pythia_new.root";

  Event ev;
  Event* pointerToEvent = &ev;
  printf("Processing sample '%s'...\n", sampleFile.Data());
  TFile input_sample(sampleFile,"READONLY");
  TTree* tree = 0;
  input_sample.GetObject("MUTREE/MUTREE",tree);
  if (!tree) input_sample.GetObject("MUTREE",tree);
  tree->SetBranchAddress("event", &pointerToEvent);

  int nentriesInTree = tree->GetEntriesFast();
  if (maxEventsUsed<0) maxEventsUsed = nentriesInTree;
  printf("\tThere are %d events in the file; running on %d events\n", nentriesInTree, maxEventsUsed);

  TH1D* hCosb = new TH1D("hCosb", "cos(#theta_{tb})", 50, -1.0, 1.0);

  printf("Input thetaEff for topPol %.3f is: %.3f\n", 1., GetThetaMixing(1., 950., 175., 425.));
  printf("Input thetaEff for topPol %.3f is: %.3f\n", 0.5, GetThetaMixing(0.5, 950., 175., 425.));
  printf("Input thetaEff for topPol %.3f is: %.3f\n", 0., GetThetaMixing(0., 950., 175., 425.));
  printf("Input thetaEff for topPol %.3f is: %.3f\n", -0.5, GetThetaMixing(-0.5, 950., 175., 425.));
  printf("Input thetaEff for topPol %.3f is: %.3f\n", -1., GetThetaMixing(-1., 950., 175., 425.));

  // Event loop
  for (int iEvent=0; iEvent<maxEventsUsed; iEvent++) {
    if (tree->LoadTree(iEvent)<0) break;
    tree->GetEntry(iEvent);

    if (ev.genInfos.size()<=0) {
      printf("This is not a MC file, EXIT!!!\n");
      return -1;
    }

    if (iEvent%1000000==0) printf("... event index %d\n", iEvent);
    
    unsigned int ngen = ev.genParticles.size();

    //double m_stop = 0.;
    //double m_chi0 = 0.;
    //double m_top = 0.;
    std::vector<SUSYGenParticle> genParticles;
    for (unsigned int ig=0; ig<ngen; ++ig) {
      GenParticle gen = ev.genParticles[ig];
      if (gen.status!=3) break;
      SUSYGenParticle part;
      part.pdgId = gen.pdgId;
      part.energy = gen.energy;
      part.pt = gen.pt;
      part.eta = gen.eta;
      part.phi = gen.phi;
      part.firstMother = -1; if (gen.mothers.size()>0) part.firstMother = gen.mothers[0];
      //if (abs(gen.pdgId)==1000006) m_stop = sqrt(pow(gen.energy,2)-pow(gen.pt*cosh(gen.eta),2));
      //if (abs(gen.pdgId)==1000022) m_chi0 = sqrt(pow(gen.energy,2)-pow(gen.pt*cosh(gen.eta),2));
      //if (abs(gen.pdgId)==6) m_top = sqrt(pow(gen.energy,2)-pow(gen.pt*cosh(gen.eta),2));

      genParticles.push_back(part);
    }
    //printf("m_stop: %.3f, m_top: %.3f, m_chi0: %.3f\n", m_stop, m_top, m_chi0);

    //double pol_new = POL;
    //double pol_new = AverageTopPolarization_Stop_to_TopChi0(-1.1, genParticles);
    //double weight = Reweight_Stop_to_TopChi0_TopOnshell (genParticles, 0., pol_new);

    // m_stop=950 GeV, m_chi0=425 GeV, m_top=175 GeV
    double thetaMixingTarget = -1.134; // Pol=-1
    //double thetaMixingTarget = -0.437; // Pol=+1
    double weight = Reweight_Stop_to_TopChi0_with_SUSYmodel (genParticles, thetaMixingTarget);

    for (unsigned int ig=0; ig<ngen; ++ig) {
      GenParticle gen = ev.genParticles[ig];
      if (gen.status!=3) break;
      if (abs(gen.pdgId)!=5) continue;
      if (gen.mothers.size()!=1) continue;
      GenParticle genTop = ev.genParticles[gen.mothers[0]];
      if (abs(genTop.pdgId)!=6) continue;

      if (genTop.pdgId*gen.pdgId<0) continue;

      double etop = genTop.energy;
      double pxtop = genTop.pt*cos(genTop.phi);
      double pytop = genTop.pt*sin(genTop.phi);
      double pztop = genTop.pt*sinh(genTop.eta);
      double ptop  = sqrt(pxtop*pxtop+pytop*pytop+pztop*pztop);
      double mtop  = sqrt(etop*etop-ptop*ptop);
      double pxb   = gen.pt*cos(gen.phi);
      double pyb   = gen.pt*sin(gen.phi);
      double pzb   = gen.pt*sinh(gen.eta);
      double pb    = sqrt(pxb*pxb+pyb*pyb+pzb*pzb);

      // We also need a stop
      if (genTop.mothers.size()==0) continue;
      GenParticle genStop = ev.genParticles[genTop.mothers[0]];
      if (abs(genStop.pdgId)!=1000006) continue;

      // Move top and fermion to the stop center-of-mass frame
      TLorentzVector stop4;
      stop4.SetPtEtaPhiE(genStop.pt, genStop.eta, genStop.phi, genStop.energy);
      TVector3 betaS(-stop4.Px()/stop4.Energy(),-stop4.Py()/stop4.Energy(),-stop4.Pz()/stop4.Energy());
      TLorentzVector topRef(pxtop,pytop,pztop,etop);
      topRef.Boost(betaS); // keept this vector to calculate costh
      TLorentzVector top4(pxtop,pytop,pztop,etop);
      top4.Boost(betaS);
      TLorentzVector b4(pxb,pyb,pzb,pb);
      b4.Boost(betaS);

      TVector3 betaV(-top4.Px()/top4.Energy(),-top4.Py()/top4.Energy(),-top4.Pz()/top4.Energy());
      top4.Boost(betaV);
      b4.Boost(betaV);

      double costh = (topRef.Px()*b4.Px()+topRef.Py()*b4.Py()+topRef.Pz()*b4.Pz())/topRef.P()/b4.P();

      hCosb->Fill(costh,weight);

    }

  }

  // To see things interactively (if IS_BATCH == false);
  TRint* app = new TRint("Wprime Analysis", &argc, argv);

  hCosb->SetMinimum(0.);
  hCosb->Draw();

  // Fitting slope
  TF1* f1 = new TF1("f1","[0]*(1+[1]*x)");
  f1->SetParName(0,"ValueAt0");
  f1->SetParName(1,"Slope");
  hCosb->Fit(f1,"","same");

  gROOT->GetListOfCanvases()->At(0)->SaveAs("costhb.jpg");

  app->Run();

  return 0;
}
