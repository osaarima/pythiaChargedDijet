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
AliJCDijetHistos *fhistosDet;
AliJCDijetAna *fana;
AliJCDijetAna *fanaMC;

double getEffFromHisto(TH1D* h, double pt);

int main(int argc, char **argv) {

    if(argc<6){
        cout<<"usage: " << argv[0] << " pythia.config pTHatMin pTHatMax dijetLeadingPt <output.root> [random_seed] [tracking inefficiency] [minJetPt]"<<endl;exit(1);
    }
    TStopwatch timer; 
    timer.Start();   

    char* pythiaconfig  = argv[1];
    double pTHatMin     = atof(argv[2]);
    double pTHatMax     = atof(argv[3]);
    double dijetLeadingPt = atof(argv[4]);
    TString outputs = argv[5];
    Int_t random_seed = argc>6 ? atoi(argv[6]) : 0;//placing the inputs into variables
    double trackingInEff = argc>7 ? atof(argv[7]) : 0.0; //Default: no tracking ineffciency
    double minJetPt = argc>8 ? atof(argv[8]) : 5.0; //Default: no tracking ineffciency


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
    vector<double> centbins = {0.0, 100.0};
    fhistos = new AliJCDijetHistos();
    fhistos->SetName("jcdijet");
    fhistos->SetCentralityBinsHistos(centbins);
    fhistos->CreateEventTrackHistos();
    fhistos->fHMG->Print();

    fhistosDet = new AliJCDijetHistos();
    fhistosDet->SetName("jcdijetDetMC");
    fhistosDet->SetCentralityBinsHistos(centbins);
    fhistosDet->CreateEventTrackHistos();
    fhistosDet->fHMG->Print();

    fana = new AliJCDijetAna();
    if(trackingInEff!=0.0) fanaMC = new AliJCDijetAna();

    TH1D *hCrossSectionInfo = new TH1D("hCrossSection","CrossSectionInfo",8,0,8);

    //------------------------------------------------------------------
    // Define jet reconstruction
    //------------------------------------------------------------------
    TClonesArray *inputList = new TClonesArray("AliJBaseTrack",1500);
    TClonesArray *inputListDet = new TClonesArray("AliJBaseTrack",1500);

    double partMinPtCut         = 0.15;// atlas 0.5 cms/alice 0.15
    double partMinEtaCut        = 0.9;
    double coneR                = 0.4; // atlas 0.6, cms 0.7 alice 0.4
    double ktconeR              = 0.4;
    double fusePionMassInktjets = false;
    double fuseDeltaPhiBGSubtr  = true;
    double jetConstituentCut    = 5.0;
    double dijetSubleadingPt    = 20.0;
    double dijetDeltaPhiCut     = 2.0; // Cut is pi/dijetDeltaPhiCut
    double fmatchingR           = 0.3;
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
                      fktScheme, //antikt
                      fusePionMassInktjets,
                      fuseDeltaPhiBGSubtr,
                      jetConstituentCut,
                      dijetLeadingPt,
                      dijetSubleadingPt,
                      minJetPt, //jet  pt cut
                      dijetDeltaPhiCut,
                      fmatchingR,
                      0.0);
    fana->InitHistos(fhistos, true, 1);

    if(trackingInEff!=0.0) {
        fanaMC->SetSettings(5,
                partMinEtaCut,
                partMinPtCut,
                coneR,
                ktconeR,
                fktScheme,
                fktScheme, //antikt
                fusePionMassInktjets,
                fuseDeltaPhiBGSubtr,
                jetConstituentCut,
                dijetLeadingPt,
                dijetSubleadingPt,
                minJetPt, //jet  pt cut
                dijetDeltaPhiCut,
                fmatchingR,
                0.0);//trackingInEff); //Efficiency is handled in this macro by DJ eff histo
        fanaMC->InitHistos(fhistosDet, true, 1);
    }


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

    TRandom3 *randomGenerator = new TRandom3(random_seed);
    /* // Information about the Efficiency histograms (DJ email):
    trackCuts = ["TPCOnly(128)","GlobalTightDCA(32)","GlobalDCA(16)","GlobalSDD(96)","Hybrid(768)","HybridStep1(256)”];
    histnames = [
        "gCor00","gEff00","gCon00"  # check it with ROOT file Title
    ];

    You will need Hybrid(768) and gCor00, ic=0(MB)

    gr = f.Get("Efficiency/{:s}{:02}{:02}".format(histnames[i],ic,IndexTrackCuts[k]));
    Therefore  “I” should be 0 and k = 4 , finally the graph you should take is "Efficiency/gCor000004” from the following root file.
    Production info :
    export TEST_DIR='/alice/sim/2017/LHC17f2b_fast/265343'
    export ALIEN_JDL_OUTPUTDIR='/alice/sim/2017/LHC17f2b_fast/265343'
    */
    TFile *fIn;
    TH1D* hCoeff;
    if(trackingInEff!=0.0){
        fIn= TFile::Open("Eff-LHC16q-pPb_MC_LHC17f2b_fast_1154_20201205-1425.root","read");
        hCoeff = (TH1D*)fIn->Get("Efficiency/hCor000004");
    }
    double fPtEff=0.0;

    for(int iEvent = 0; iEvent < nEvent; ++iEvent) {//begin event loop

        if (!pythia.next()) continue;
        inputList->Clear("C");
        inputListDet->Clear("C");
        nTried = pythia.info.nTried();
        nTrial = nTried - prev_nTried;
        prev_nTried = nTried;
        sigmaGen = pythia.info.sigmaGen();
        ebeweight = 1.0; //no event-by-event weight at all. //sigmaGen/nTrial;
        hCrossSectionInfo->Fill(7.5,ebeweight);
        if(iEvent % ieout == 0) cout << iEvent << "\t" << int(float(iEvent)/nEvent*100) << "%, nTried:" << nTried << ", nTrial:" << nTrial << ", sigma:" << sigmaGen << endl;

        for (int i = 0; i < pythia.event.size(); ++i) {//loop over all the particles in the event
            if (pythia.event[i].isFinal() && pythia.event[i].isCharged() && pythia.event[i].isHadron() ) { // Only check if it is final, charged and hadron since the acceptance is checked in the CalculateJetsDijets

                TLorentzVector lParticle(pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e());
                AliJBaseTrack track( lParticle );
                track.SetID(pythia.event[i].id());
                track.SetTrackEff(1.);
                new ((*inputList)[inputList->GetEntriesFast()]) AliJBaseTrack(track);

                if (trackingInEff!=0.0) {
                    fPtEff = getEffFromHisto(hCoeff, TMath::Sqrt(pythia.event[i].px()*pythia.event[i].px()+pythia.event[i].py()*pythia.event[i].py()));
                    if(fPtEff < randomGenerator->Rndm()) {
                        new ((*inputListDet)[inputListDet->GetEntriesFast()]) AliJBaseTrack(track);
                    }
                }
            }
        } // end of finalparticles

        // Here I call my function
        fana->CalculateJets(inputList, fhistos, fCBin);
        fana->FillJetsDijets(fhistos, fCBin);

        if(trackingInEff!=0.0) {
            fanaMC->CalculateJets(inputListDet, fhistosDet, fCBin);
            fanaMC->FillJetsDijets(fhistosDet, fCBin);

            // Here response matrix calculation.
            fana->CalculateResponse(fanaMC,fhistosDet);
        }

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

//Search pt bin from efficiency histo and return the content
//If pt is less than first bin then return first bin.
//If pt is more than last bin then return last bin.
double getEffFromHisto(TH1D* h, double pt) {
    int iBin;
    int iLastBin = h->GetNbinsX();

    if(pt<h->GetBinLowEdge(1)) return h->GetBinContent(1);
    if(pt>h->GetBinLowEdge(iLastBin)) return h->GetBinContent(iLastBin);
    for(int i=1; i<iLastBin; i++) {
        if(pt>h->GetBinLowEdge(i) && pt<h->GetBinLowEdge(i+1)) {
            if(h->GetBinContent(i)>1.0) cout << "Warning: efficiency>1.0!" << endl;
            if(h->GetBinContent(i)<0.0) cout << "Warning: efficiency<0.0!" << endl;
            return h->GetBinContent(i);
        }

    }
    cout << "Warning: efficiency not found!" << endl;
    return DBL_MAX;
}
