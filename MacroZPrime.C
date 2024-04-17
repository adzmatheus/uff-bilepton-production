/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the muon invariant
mass.

root -l examples/Example1.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------

void MacroZPrime(const char *inputFile) {
    gSystem->Load("libDelphes");

    int quant = 0, i = 0;

    // Create chain of root trees
    TChain chain("Delphes");
    chain.Add(inputFile);

    // Create object of class ExRootTreeReader
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();

    // Get pointers to branches used in this analysis
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");

    float maxPT = 0;
    float maxMass = 0;
    float maxEta = 0;
    float maxPhi = 0;
    float minEta = 0;
    float minPhi = 0;
    float MaxPT1 = 0;
    float MaxPT2 = 0;
    int MaxPT1_Pos = -1;
    int MaxPT2_Pos = -1;


    //cutflow
    float m_eta;
    int total_muon = 0;
    int quant_filtro = 0;
    int quant_nao_filtro = 0;
    int filtro_eta = 0;
    int filtro_pt = 0;


    // Finding maxium values for PT and Invariant Mass
    for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
        // Load selected branches with data from specified event
        treeReader->ReadEntry(entry);

        // If event contains at least 1 muon
        if(branchMuon->GetEntries() > 0)
        {
            quant = branchMuon->GetEntries();

            total_muon += branchMuon->GetEntries();
            
            // Take first muon
            for (i=0; i<quant; i++) {
                Muon *muon = (Muon*) branchMuon->At(i);

                if (muon->PT > maxPT) {
                    maxPT = muon->PT;
                }
                if (muon->Eta > maxEta) {
                    maxEta = muon->Eta;
                }
                if (muon->Phi > maxPhi) {
                    maxPhi = muon->Phi;
                }
                if (muon->Phi < minPhi) {
                    minPhi = muon->Phi;
                }
                if (muon->Eta < minEta) {
                    minEta = muon->Eta;
                }
            }

            Muon *mu1, *mu2;
            if(branchMuon->GetEntries() > 1)
            {
                // Take first two muons
                mu1 = (Muon *) branchMuon->At(0);
                mu2 = (Muon *) branchMuon->At(1);

                // Plot their invariant mass
                if(((mu1->P4()) + (mu2->P4())).M() > maxMass) {
                    maxMass = ((mu1->P4()) + (mu2->P4())).M();
                }
            }
        }

    }

    cout<<"MaxPT = "<<maxPT<<endl;
    cout<<"MaxMass = "<<maxMass<<endl;
    cout<<"MaxPhi = " <<maxPhi<<endl;
    cout<<"MaxEta = " <<maxEta<<endl;

    // Book histograms
    TH1 *histPT = new TH1F("Pt", "Pt", 100, 0.0, maxPT + 100);
    TH1 *histMass = new TH1F("InvMass", "InvMass", 100, 0.0, maxMass + 100);
    TH1 *histPhi = new TH1F ("Phi", "Phi", 100, minPhi - 1, maxPhi + 1 );
    TH1 *histEta = new TH1F ("Eta", "Eta", 100, minEta - 1, maxEta + 1 );

    // Loop over all events
    for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
        // Load smuted branches with data from specified event
        treeReader->ReadEntry(entry);

        // If event contains at least 1 muon
        if(branchMuon->GetEntries() > 0)
        {
            quant = branchMuon->GetEntries();
            // Take first muon
            for (i=0; i<quant; i++) {
                Muon *muon = (Muon*) branchMuon->At(i);

                if (muon->PT > 30){
                    filtro_pt++;
                }
                if (m_eta < 2.5){
                    filtro_eta++;
                }

                if ((muon->PT > 30)&&(fabs(muon->Eta) < 2.5)) {
                    quant_filtro++;
                    // Plot transverse momentum
                    histPT->Fill(muon->PT);
                    // Plot Phi
                    histPhi->Fill(muon->Phi);
                    // Plot Eta
                    histEta->Fill(muon->Eta);

                }else{
                    quant_nao_filtro++;
                }
            }



            Muon *mu1, *mu2;

            // If event contains at least 2 muons
            if(branchMuon->GetEntries() > 1)
            {
                Muon *muon1;
                for (i=0; i<quant; i++) {
                    muon1 = (Muon*) branchMuon->At(i);

                    if ((muon1->PT > 30)&&(fabs(muon1->Eta) < 2.5)) {
                        if (muon1->PT > MaxPT1) {
                            MaxPT1 = muon1->PT;
                            MaxPT1_Pos = i;
                        }
                    }
                }
                

                if(MaxPT1_Pos != -1){
                muon1 = (Muon *) branchMuon->At(MaxPT1_Pos);
                }

                Muon *muon2;

                for (i=0; i<quant; i++) {
                    muon2 = (Muon*) branchMuon->At(i);

                    if ((muon2->PT > 30)&&(fabs(muon2->Eta) < 2.5)) {
                        if (i!=MaxPT1_Pos) {
                            if (muon2->PT > MaxPT2) {
                                if((muon2->Charge + muon1->Charge) == 0){
                                    MaxPT2 = muon2->PT;
                                    MaxPT2_Pos = i;
                                }
                            }
                        }
                    }
                }

                cout<<MaxPT1<<" -- "<<MaxPT2<<endl;

                if( (MaxPT1_Pos != -1) && (MaxPT2_Pos != -1) ) {
                    // Take first two muons
                    mu1 = (Muon *) branchMuon->At(MaxPT1_Pos);
                    mu2 = (Muon *) branchMuon->At(MaxPT2_Pos);

                    // Plot their invariant mass
                    histMass->Fill(((mu1->P4()) + (mu2->P4())).M());
                }
            }

            MaxPT1 = 0;
            MaxPT2 = 0;
            MaxPT1_Pos = -1;
            MaxPT2_Pos = -1;

        }

    }

    // Show resulting histograms
    TCanvas *PT1 = new TCanvas("PT1","PT1",200,100,550,550);
    histPT->Draw();
    TCanvas *Phi = new TCanvas("Phi","Phi",200,100,550,550);
    histPhi->Draw();
    TCanvas *Eta = new TCanvas("Eta","Eta",200,100,550,550);
    histEta->Draw();
    TCanvas *InvMass = new TCanvas("InvMass","InvMass",200,100,550,550);
    histMass->Draw();



    cout<<"Total of muons: "<<total_muon<<endl;
    cout<<"Total of muons with Pt > 30: "<<filtro_pt<<endl;
    cout<<"Total of muons with |Eta| < 2.5: "<<filtro_eta<<endl;
    cout<<"Total of muons within the requirements: "<<quant_filtro<<endl;
    cout<<"Total of muons out the requirements: "<<quant_nao_filtro<<endl;
}

