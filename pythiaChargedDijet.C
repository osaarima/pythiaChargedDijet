#include "TH1D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TFile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include <Pythia8/Pythia.h>
#include <TStopwatch.h>

/*
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
*/
/*
#include "fastjet/config.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include <fastjet/SISConePlugin.hh>
#include <fastjet/CDFMidPointPlugin.hh>
#ifdef FASTJET_VERSION
#include <fastjet/Selector.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <fastjet/tools/Subtractor.hh>
*/
#include <iostream>

#include <TClonesArray.h>
#include "src/AliJCDijetHistos.h"
#include "src/AliJCDijetAna.h"
#include "src/AliJBaseTrack.h"
//#include "src/JHistos.h"
//#include "src/AliJCard.h"

#define DEBUG 0

using namespace std;
using namespace Pythia8; 
class AliJCDijetHistos;

AliJCDijetHistos *fhistos;
AliJCDijetAna *fana;


int main(int argc, char **argv) {

    if(argc<6){
        cout<<"usage: " << argv[0] << " pythia.config pTHatMin pTHatMax dijetLeadingPt <output.root> [random_seed] [tracking inefficiency]"<<endl;exit(1);
    }
    TStopwatch timer; 
    timer.Start();   

    char* pythiaconfig  = argv[1];
    double pTHatMin     = atof(argv[2]);
    double pTHatMax     = atof(argv[3]);
    double dijetLeadingPt = atof(argv[4]);
    TString outputs = argv[5];
    Int_t random_seed = argc>6 ? atoi(argv[6]) : 0;//placing the inputs into variables
    unsigned trackingInEff = argc>7 ? atoi(argv[7]) : 0; //Default: no tracking ineffciency


    TFile *fout = new TFile(outputs.Data(),"RECREATE");
    fout->cd();//opening of the output file
    TDirectoryFile *fdir = new TDirectoryFile( "JCDijetBaseTask","JCDijetBaseTask" );
    fdir->cd();

    //---------------------
    //Pythia initialization 
    //---------------------
    Pythia pythia;   // Generator.
    //Event& event      = pythia.event;
    ParticleData& pdt = pythia.particleData;

    // Read in commands from external file.
    pythia.readFile(pythiaconfig);   

    // Extract settings to be used in the main program.
    int    nEvent  = pythia.mode("Main:numberOfEvents");
    bool   showCPD = pythia.flag("Main:showChangedParticleData");
    double energy  = pythia.mode("Beams:eCM");

    pythia.readString(Form("PhaseSpace:pTHatMin ==%f",pTHatMin));
    pythia.readString(Form("PhaseSpace:pTHatMax ==%f",pTHatMax));
    cout<<"Events="<<nEvent <<" RNDM seed "<< random_seed << endl;

    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed=%02d",random_seed));

    // Initialize. Beam parameters set in .cmnd file.
    pythia.init();

    // List changed data. 
    if (showCPD) pdt.listChanged();

    //-------------------------------------------------------
    // Histograms and tools
    //-------------------------------------------------------
    fhistos = new AliJCDijetHistos();
    vector<double> centbins = {0.0, 100.0};
    fhistos->SetCentralityBinsHistos(centbins);
    fhistos->CreateEventTrackHistos();

    fhistos->fHMG->Print();

    fana = new AliJCDijetAna();

    TH1D *hCrossSectionInfo = new TH1D("hCrossSection","CrossSectionInfo",8,0,8);

    //------------------------------------------------------------------
    // Define jet reconstruction
    //------------------------------------------------------------------
    TClonesArray *inputList = new TClonesArray("AliJBaseTrack",1500);

    double partMinPtCut         = 0.15;// atlas 0.5 cms/alice 0.15
    double partMinEtaCut        = 0.8;
    double coneR                = 0.4; // atlas 0.6, cms 0.7 alice 0.4
    double ktconeR              = 0.4;
    double fusePionMassInktjets = false;
    double fuseDeltaPhiBGSubtr  = true;
    double jetConstituentCut    = 5.0;
    double dijetSubleadingPt    = 20.0;
    double dijetDeltaPhiCut     = 2.0; // Cut is pi/dijetDeltaPhiCut
    double fmatchingR           = 0.2;
    int fktScheme               = 1;
    int centBin=0;

    TString sktScheme;
    switch (fktScheme) {
        case 0:  sktScheme = "E_scheme";
                 break;
        case 1:  sktScheme = "pt_scheme";
                 break;
        case 2:  sktScheme = "pt2_scheme";
                 break;
        case 3:  sktScheme = "Et_scheme";
                 break;
        case 4:  sktScheme = "Et2_scheme";
                 break;
        case 5:  sktScheme = "BIpt_scheme";
                 break;
        case 6:  sktScheme = "BIpt2_scheme";
                 break;
        default: sktScheme = "Unknown, check macro arguments!";
                 break;
    }

    cout << endl;
    cout << "============= Settings =============" << endl;
    cout << "cent bin:                   " << centBin << endl;
    cout << "particle eta cut:           " << partMinEtaCut << endl;
    cout << "particle pt cut:            " << Form("%.2f",partMinPtCut) << endl;
    cout << "jet cone size:              " << coneR << endl;
    cout << "kt-jet cone size:           " << ktconeR << endl;
    cout << "Using pion mass in kt-jets: " << fusePionMassInktjets << endl;
    cout << "Using DeltaPhi in BG subtr: " << fuseDeltaPhiBGSubtr << endl;
    cout << "Using kt-jet scheme:        " << sktScheme.Data() << endl;
    cout << "jet costituent cut:         " << jetConstituentCut << endl;
    cout << "dijet leading pt cut:       " << dijetLeadingPt << endl;
    cout << "dijet subleading pt cut:    " << dijetSubleadingPt << endl;
    cout << "dijet DeltaPhi cut:         pi/" << dijetDeltaPhiCut << endl;
    cout << "tracking inefficiency:      " << trackingInEff << endl;
    cout << endl;

    if(fusePionMassInktjets && fktScheme!=0) {
        cout << "Warning: Using pion mass for kt-jets but not using E_scheme!" << endl;
        cout << endl;
    }

    fana->SetSettings(5,
                      partMinEtaCut,
                      partMinPtCut,
                      coneR,
                      ktconeR,
                      fktScheme,
                      fusePionMassInktjets,
                      fuseDeltaPhiBGSubtr,
                      jetConstituentCut,
                      dijetLeadingPt,
                      dijetSubleadingPt,
                      dijetDeltaPhiCut,
                      fmatchingR);
    fana->InitHistos(fhistos, true, 1);


    //--------------------------------------------------------
    //         B e g i n    e v e n t    l o o p.
    //--------------------------------------------------------
    cout<<"Let's start" <<endl; 
    int ieout = nEvent/20;
    if (ieout<1) ieout=1;
    int EventCounter = 0;
    int fCBin = 0;
    Int_t nTried = 0; 
    Int_t prev_nTried = 0;
    Int_t nTrial = 0;
    Int_t nAccepted = 0;
    Float_t sigmaGen = 0.0;
    Float_t ebeweight = 1.0;
    TRandom3 *randomGenerator = new TRandom3();

    for(int iEvent = 0; iEvent < nEvent; ++iEvent) {//begin event loop

        if (!pythia.next()) continue;
        inputList->Clear("C");
        nTried = pythia.info.nTried();
        nTrial = nTried - prev_nTried;
        prev_nTried = nTried;
        sigmaGen = pythia.info.sigmaGen();
        ebeweight = 1.0; //no event-by-event weight at all. //sigmaGen/nTrial;
        hCrossSectionInfo->Fill(7.5,ebeweight);
        if(iEvent % ieout == 0) cout << iEvent << "\t" << int(float(iEvent)/nEvent*100) << "%, nTried:" << nTried << ", nTrial:" << nTrial << ", sigma:" << sigmaGen << endl;

        for (int i = 0; i < pythia.event.size(); ++i) {//loop over all the particles in the event
            if (pythia.event[i].isFinal() && pythia.event[i].isCharged() && pythia.event[i].isHadron() ) { // Only check if it is final, charged and hadron since the acceptance is checked in the CalculateJetsDijets
                if (trackingInEff!=0 && randomGenerator->Uniform(0.0,1.0) < 0.03) continue;// Lose 3% of tracks, as in ALICE.

                TLorentzVector lParticle(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
                AliJBaseTrack track( lParticle );
                track.SetID(pythia.event[i].id());
                track.SetTrackEff(1.);
                new ((*inputList)[inputList->GetEntriesFast()]) AliJBaseTrack(track);
            }
        } // end of finalparticles

        // Here I call my function
        fana->CalculateJets(inputList, fhistos, fCBin);
        fana->FillJetsDijets(fhistos, fCBin);

        EventCounter++;
        if(iEvent == nEvent-1) cout << nEvent << "\t" << "100%, nTried:" << pythia.info.nTried() << ", sigma:" << pythia.info.sigmaGen() << endl ;
    }//event loop

    nTried = pythia.info.nTried();
    nAccepted = pythia.info.nAccepted();
    sigmaGen = pythia.info.sigmaGen();
    //double sigmaErr = pythia.info.sigmaErr();
    hCrossSectionInfo->Fill(0.5,nTried);
    cout << "nTried after loop:" << nTried << endl;// print also inside event loop and in the macro.
    hCrossSectionInfo->Fill(1.5,nAccepted);
    //cout << "nAccepted after loop:" << nAccepted << endl;
    hCrossSectionInfo->Fill(2.5,sigmaGen);
    cout << "sigma after loop:" << sigmaGen << endl;
    hCrossSectionInfo->Fill(3.5,EventCounter);
    hCrossSectionInfo->Fill(4.5,energy);
    hCrossSectionInfo->Fill(5.5,1); // for counting # of merged
    hCrossSectionInfo->Fill(6.5,pythia.info.weightSum()); // for counting # of merged

    fout->Write();
    fout->Close();
    cout << EventCounter << " events are analyzed successfully."<< endl;
    timer.Print(); 
    return 0;
}
