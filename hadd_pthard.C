#include "../AliJCDijet-Results/DijetHelper.h"
const int Nfiles = 10;
TString sOutFolder = "trackingEffDJRun_closure";
TString sDefName = "seed2000";
TString inputfiles[Nfiles]={
    Form("output/%s/data/11-21%s_%s.root",sOutFolder.Data(),sOutFolder.Data(),sDefName.Data()),
	Form("output/%s/data/21-36%s_%s.root",sOutFolder.Data(),sOutFolder.Data(),sDefName.Data()),
	Form("output/%s/data/36-57%s_%s.root",sOutFolder.Data(),sOutFolder.Data(),sDefName.Data()),
	Form("output/%s/data/57-84%s_%s.root",sOutFolder.Data(),sOutFolder.Data(),sDefName.Data()),
	Form("output/%s/data/84-117%s_%s.root",sOutFolder.Data(),sOutFolder.Data(),sDefName.Data()),
	Form("output/%s/data/117-152%s_%s.root",sOutFolder.Data(),sOutFolder.Data(),sDefName.Data()),
	Form("output/%s/data/152-191%s_%s.root",sOutFolder.Data(),sOutFolder.Data(),sDefName.Data()),
	Form("output/%s/data/191-234%s_%s.root",sOutFolder.Data(),sOutFolder.Data(),sDefName.Data()),
	Form("output/%s/data/234-300%s_%s.root",sOutFolder.Data(),sOutFolder.Data(),sDefName.Data()),
	Form("output/%s/data/300--1%s_%s.root",sOutFolder.Data(),sOutFolder.Data(),sDefName.Data())
};

void doweight(TString inFile, TString outfile,double scale);

void hadd_pthard() {

    for(int i=0;i<Nfiles;i++) {
        //for(int i=0;i<1;i++) {
        cout << "===================================================" << endl;
        cout << "Weight for i = " << i << endl << endl << endl << endl;
        doweight(inputfiles[i],Form("output/%s/data/Normalized_AnalysisResults_%d.root",sOutFolder.Data(),i+1),1.0);
        //}
    }
    cout << "Scaled with sigma_gen/N_evt" << endl;
}

void doweight(TString inFile, TString outfile,double scale) {
    cout << inFile <<endl;
    TFile *fin = TFile::Open(Form("%s",inFile.Data()));
    fin->Print();
    TFile *fout = new TFile(outfile,"recreate");

    //TH1D *hCrossSectionInfo = (TH1D*)fin->Get("hCrossSection");
    /*double  sigmaGen = hCrossSectionInfo->GetBinContent(3);
      int     nAccepted = hCrossSectionInfo->GetBinContent(2);
      int     nTried = hCrossSectionInfo->GetBinContent(1);
      int     nMerge = hCrossSectionInfo->GetBinContent(6);
      int     wSum = hCrossSectionInfo->GetBinContent(7);*/
    //double  normPythia = sigmaGen/nTried/nMerge; 
    double normPythia = scale;
    int const nDirs = 2;

    TString sBaseFolder = "JCDijetBaseTask";
    TString sFolder[nDirs] = {"jcdijet"
                         ,"jcdijetDetMC"};
    TDirectory *dir[nDirs];


    TString sEvents[nDirs] = {"JCDijetBaseTask/jcdijet/h_events/h_eventsCentBin00"
                              ,"JCDijetBaseTask/jcdijetDetMC/h_events/h_eventsCentBin00"
                              };

    TDirectoryFile* newBaseDir = new TDirectoryFile("JCDijetBaseTask","JCDijetBaseTask");
    for(int iFol=0; iFol<2; iFol++) {
        dir[iFol] = fin->GetDirectory(Form("%s/%s",sBaseFolder.Data(),sFolder[iFol].Data()));

        TH1D *hCrossSection = getHisto(fin,"JCDijetBaseTask/hCrossSection");
        if(hCrossSection==0) ErrorExit(Form("No hCrossSection found"));

        int Nevts_all= hCrossSection->GetBinContent(8); // # of event per segement times nsegement
        int nseq = hCrossSection->GetBinContent(6);
        double ntried = hCrossSection->GetBinContent(1);
        double sigma_all = hCrossSection->GetBinContent(3); // in mb  sigma_gen x nseq
        double sigma = sigma_all/nseq;
        double scaling = sigma/(double)Nevts_all; 

        newBaseDir->cd();
        TDirectoryFile* newHistoDir = new TDirectoryFile(sFolder[iFol].Data(),sFolder[iFol].Data());
        newHistoDir->cd();
        TString foutdir = Form("%s/%s",sBaseFolder.Data(),sFolder[iFol].Data());
        //TDirectory *newdir = fout->mkdir(Form("%s/%s",sBaseFolder.Data(),sFolder[iFol].Data()));
        if (dir[iFol]) {
            dir[iFol]->cd();
            TIter next(gDirectory->GetListOfKeys());
            TKey * key;
            TString olds;
            while( (key=(TKey*)next()) ){
                //key->Print();
                dir[iFol]->cd();
                TString sname=key->GetName();
                //cout << sname << endl;
                TH1 *h = (TH1*)key->ReadObj();
                //if(sname == "hCrossSection") continue;
                //if(sname=="hiCentr") continue;
                //if(sname =="HistManager") continue;
                //if(sname =="JCard") continue;
                //if(sname =="h2DFullEvsChEdN0") continue;
                //if(sname =="h2DFullEvsChEdNnot0") continue;

                //h->Print();
                h->Scale(scaling);
                newHistoDir->cd();
                h->Write();
                //TDirectory *subdir = fin->GetDirectory(Form("AliJDiJetTask/AliJDiJetAnalysisHistManager/%s",sname));
                TDirectory *subdir = fin->GetDirectory(Form("%s/%s",foutdir.Data(),sname.Data()));
                //subdir->Print();
                if(subdir){
                    //subdir->ls(); // Use this to print directory nicely.
                    newHistoDir->cd();
                    TDirectory *newdir = new TDirectoryFile(sname.Data(),sname.Data());
                    newdir->cd();
                    subdir->cd();
                    //cout << gDirectory->GetListOfKeys() << endl;
                    TIter next2(gDirectory->GetListOfKeys());
                    TKey *key2;
                    while( (key2 = (TKey*)next2()) ){
                        TString sname=key2->GetName();
                        //cout << sname << endl;
                        TH1 *h = (TH1*)key2->ReadObj();
                        //h->Sumw2();
                        //if(sname == "hCrossSection") continue;
                        //if(sname=="hiCentr") continue;
                        //h->Print();
                        h->Scale(scaling);

                        newdir->cd();
                        h->Write();
                        //TDirectory *subdir = fin->GetDirectory(sname);
                    }
                }
            }
            fout->cd();
            //(TH1D*)fin->Get("hCrossSection")->Write();
            //(TH1D*)fin->Get("hiCentr")->Write();
        } else {
            cout <<"No dir!" << endl;
        }
    }

    fout->Close();
    fin->Close();
}
