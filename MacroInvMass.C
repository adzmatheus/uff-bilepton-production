Trocar sua página inicial? … 
No momento, ela está definida como "Início". Você pode mudar essa opção a qualquer momento nas configurações.
MacroInvMass.C

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

void MacroInvMass(const char *inputFile,const char *outputFile)
{
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

    float MaxPT1 = 0;
    float MaxPT2 = 0;
    int MaxPT1_Pos = -1;
    int MaxPT2_Pos = -1;
    float M_Inv;
    int quant_mu = 0; //amount of events with 2 or more muons

    //cutflow
    int quant_filtro = 0;
    int quant_nao_filtro = 0;
    int total_muon = 0;
    

    TFile f(outputFile,"recreate");
    TTree t1("InvMass","MassaInvariante");
    t1.Branch("InvMass",&M_Inv, "InvMass");


    // Loop over all events
    for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    {
        // Load smuted branches with data from specified event
        treeReader->ReadEntry(entry);

        // If event contains at least 2 muon
        if(branchMuon->GetEntries() > 1)
        {
            quant = branchMuon->GetEntries();
            Muon *mu1, *mu2;
            quant_mu++;
            total_muon += branchMuon->GetEntries();

            // Finding the muon with higher Pt in the event, storing its index in MaxPT1_Pos

            Muon *muon1;
            for (i=0; i<quant; i++) {
                muon1 = (Muon*) branchMuon->At(i);

                if ((muon1->PT > 30)&&(fabs(muon1->Eta) < 2.5)) {
                    quant_filtro++;
                    if (muon1->PT > MaxPT1) {
                        MaxPT1 = muon1->PT;
                        MaxPT1_Pos = i;
                    }
                }else{
                    quant_nao_filtro++;
                }
            }

            if(MaxPT1_Pos != -1){
                muon1 = (Muon *) branchMuon->At(MaxPT1_Pos);
            }

            // Finding the muon with second higher Pt in the event, storing its index in MaxPT2_Pos

            Muon *muon2;
            for (i=0; i<quant; i++) {
                muon2 = (Muon*) branchMuon->At(i);

                if ((muon2->PT > 30)&&(fabs(muon2->Eta) < 2.5)) {
                    if (i!=MaxPT1_Pos) { // If the muon being analyzed is not the one found in the step before 
                        if (muon2->PT > MaxPT2) {
                            if((muon2->Charge + muon1->Charge) == 0){
                                MaxPT2 = muon2->PT;
                                MaxPT2_Pos = i;
                            }
                        }
                    }
                }
            }



            // If there are 2 suitable muons in the event
            if ( (MaxPT1_Pos != -1) && (MaxPT2_Pos != -1) ) {
                mu1 = (Muon *) branchMuon->At(MaxPT1_Pos);
                mu2 = (Muon *) branchMuon->At(MaxPT2_Pos);

                M_Inv = ((mu1->P4()) + (mu2->P4())).M();
                t1.Fill();
            }


            MaxPT1 = 0;
            MaxPT2 = 0;
            MaxPT1_Pos = -1;
            MaxPT2_Pos = -1;

        }

    }
    t1.Write();
    cout<<"Events with at least 2 muons: "<<quant_mu<<endl;
    cout<<"Total of events: "<<numberOfEntries<<endl;
    cout<<"Events with at least 2 muons: "<<quant_mu<<endl<<endl;
    cout<<"From the events with 2 or more muons: "<<endl;
    cout<<"- Total of muons: "<<total_muon<<endl;
    cout<<"- Total of muons within the requirements: "<<quant_filtro<<endl;
    cout<<"- Total of muons out the requirements: "<<quant_nao_filtro<<endl;
}


