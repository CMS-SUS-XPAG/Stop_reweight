#include "T2bW_reweighting.C" 

#include "TStyle.h"
#include "TRint.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "MuTree.h"

int maxEventsUsed = -1; // If <0 ==> use all events in the file

using namespace std;
using namespace ciemat;

int main(int argc, char** argv){

  gStyle->SetOptStat(1111111);

  // stop file
  TString sampleFile = "./SMS-T2bw_Mstop-700_mChargino-570_mLSP-150_8TeV.root";

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

  TH1D* hptleptonL = new TH1D("hptleptonL", "Lepton p_{T} [GeV]", 30, 0.0, 1000.0);
  hptleptonL->SetTitle("m(#tilde{t})=700 GeV, m(#tilde{#chi}^{#pm})=570 GeV, m(#tilde{#chi}^{0})=150 GeV");
  hptleptonL->GetXaxis()->CenterTitle();
  hptleptonL->SetXTitle("Lepton p_{T} [GeV]");
  TH1D* hptleptonR = new TH1D("hptleptonR", "Lepton p_{T} [GeV]", 30, 0.0, 1000.0);
  TH1D* hptleptonPythia = new TH1D("hptleptonPythia", "Lepton p_{T} [GeV]", 30, 0.0, 1000.0);

  TH1D* hptmissL = new TH1D("hptmissL", "Missing p_{T} [GeV]", 30, 0.0, 1000.0);
  hptmissL->SetTitle("m(#tilde{t})=700 GeV, m(#tilde{#chi}^{#pm})=570 GeV, m(#tilde{#chi}^{0})=150 GeV");
  hptmissL->GetXaxis()->CenterTitle();
  hptmissL->SetXTitle("Missing p_{T} [GeV]");
  TH1D* hptmissR = new TH1D("hptmissR", "Missing p_{T} [GeV]", 30, 0.0, 1000.0);
  TH1D* hptmissPythia = new TH1D("hptmissPythia", "Missing p_{T} [GeV]", 30, 0.0, 1000.0);

  // Event loop
  for (int iEvent=0; iEvent<maxEventsUsed; iEvent++) {
    if (tree->LoadTree(iEvent)<0) break;
    tree->GetEntry(iEvent);

    if (iEvent%1000000==0) printf("... event index %d\n", iEvent);
    
    unsigned int ngen = ev.genParticles.size();

    //printf("\n>>>>>>>>>>>>>>>>>>>>>\n");
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
      
      genParticles.push_back(part);

      /*
      double mass = pow(part.energy,2);
      if (part.pt>0) mass -= pow(part.pt*cosh(part.eta),2);
      if (mass>0) mass = sqrt(mass); else mass = -sqrt(-mass);
      printf("ig: %3d, id: %10d, mother: %3d, pt: %.3f, eta: %.3f, phi: %.3f, mass: %.3e\n", ig, part.pdgId, part.firstMother, part.pt, part.eta, part.phi, mass);
      */

    }

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
      if (i_b<0 || i_chargino<0) break;

      int i_neutralino = -1;
      int i_W = -1;
      for (unsigned int ig=i_chargino+1; ig<ngen; ++ig) {
            const SUSYGenParticle& gen = genParticles[ig];
            if (abs(gen.firstMother)!=i_chargino) continue;
            if (abs(gen.pdgId)==24) i_W = ig;
            else if (abs(gen.pdgId)==1000022) i_neutralino = ig;
            if (i_W>=0 && i_neutralino>=0) break;
      }
      if (i_W<0 || i_neutralino<0) break;

      int i_up = -1;
      int i_down = -1;
      for (unsigned int ig=i_W+1; ig<ngen; ++ig) {
            const SUSYGenParticle& gen = genParticles[ig];
            if (abs(gen.firstMother)!=i_W) continue;
            if (abs(gen.pdgId)%2==0) i_up = ig;
            else if (abs(gen.pdgId)%2==1) i_down = ig;
            if (i_up>=0 && i_down>=0) break;
      }
      if (i_up<0 || i_down<0) break;

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

      double weight_L = Reweight_T2bW(M_PI/2, M_PI/2, genParticles);
      double weight_R = Reweight_T2bW(M_PI/2, 0., genParticles);

      // Pt of lepton
      hptleptonL->Fill(down4.Pt(), weight_L);
      hptleptonR->Fill(down4.Pt(), weight_R);
      hptleptonPythia->Fill(down4.Pt(), 1.);

      // Missing Pt
      int i_otherneutralino = -1;
      for (unsigned int ig=0; ig<ngen; ++ig) {
            const SUSYGenParticle& gen = genParticles[ig];
            if (ig!=i_neutralino && abs(gen.pdgId)==1000022) i_otherneutralino = ig;
            if (i_otherneutralino>=0) break;
      }
      if (i_otherneutralino<0) continue;
      const SUSYGenParticle& gen_otherneutralino = genParticles[i_otherneutralino];
      TLorentzVector otherneutralino4;
      otherneutralino4.SetPtEtaPhiE(gen_otherneutralino.pt, gen_otherneutralino.eta, gen_otherneutralino.phi, gen_otherneutralino.energy);
      TLorentzVector miss4 = up4 + neutralino4 + otherneutralino4;
      hptmissL->Fill(miss4.Pt(), weight_L);
      hptmissR->Fill(miss4.Pt(), weight_R);
      hptmissPythia->Fill(miss4.Pt(), 1.);

    }
  }

  // To see things interactively
  TRint* app = new TRint("Direct Stop Search Analysis", &argc, argv);

  TCanvas* c1 = new TCanvas("c1","My Canvas", 1024, 768);
  TPad* pad1 = new TPad("pad1","This is pad1",0.05,0.52,0.95,0.97);
  TPad* pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.47);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  hptleptonL->SetLineStyle(1); hptleptonL->SetLineWidth(3);
  hptleptonL->SetMaximum(1.2*hptleptonL->GetMaximum());
  hptleptonL->Draw("hist"); 
  hptleptonR->SetLineStyle(3); hptleptonR->SetLineWidth(3);
  hptleptonR->Draw("histsame");
  hptleptonPythia->SetLineStyle(2); hptleptonPythia->SetLineWidth(3);
  hptleptonPythia->Draw("histsame");
  TLegend* legA = new TLegend(0.52,0.3,0.9,0.6);
  legA->SetFillColor(0);
  legA->AddEntry(hptleptonL,"RH #tilde{#chi}^{#pm}, L W#tilde{#chi}^{#pm}#tilde{#chi}^{0} coupling","L");
  legA->AddEntry(hptleptonR,"RH #tilde{#chi}^{#pm}, R W#tilde{#chi}^{#pm}#tilde{#chi}^{0} coupling","L");
  legA->AddEntry(hptleptonPythia,"Our reference","L");
  legA->Draw();

  pad2->cd();
  hptmissL->SetLineStyle(1); hptmissL->SetLineWidth(3);
  hptmissL->SetMaximum(1.2*hptmissL->GetMaximum());
  hptmissL->Draw("hist");
  hptmissR->SetLineStyle(3); hptmissR->SetLineWidth(3);
  hptmissR->Draw("histsame"); 
  hptmissPythia->SetLineStyle(2); hptmissPythia->SetLineWidth(3);
  hptmissPythia->Draw("histsame");
  TLegend* legB = new TLegend(0.52,0.3,0.9,0.6);
  legB->SetFillColor(0);
  legB->AddEntry(hptmissL,"RH #tilde{#chi}^{#pm}, L W#tilde{#chi}^{#pm}#tilde{#chi}^{0} coupling","L");
  legB->AddEntry(hptmissR,"RH #tilde{#chi}^{#pm}, R W#tilde{#chi}^{#pm}#tilde{#chi}^{0} coupling","L");
  legB->AddEntry(hptmissPythia,"Our reference","L");
  legB->Draw();

  c1->SaveAs("t2bw_referenceTest.jpg");

  app->Run();

  return 0;
}
