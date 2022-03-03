
double getEffFromHisto(TH1D* h, double pt);

int main(int argc, char **argv) {
    TFile *fIn;
    TH1D* hCoeff;
    fIn= TFile::Open("Eff-LHC16q-pPb_MC_LHC17f2b_fast_1154_20201205-1425.root","read");
    hCoeff = (TH1D*)fIn->Get("Efficiency/hCor000004");
    double fPtEff=0.0;

    TH1D* his=new TH1D("his","his",100,0,100);

    double pt = 10;

    fPtEff = getEffFromHisto(his, pt);

    if(fPtEff > randomGenerator->Rndm()) {
        //Accept particle
    }
}

//Search pt bin from efficiency histo and return the content
//If pt is less than first bin then return first bin.
//If pt is more than last bin then return last bin.
double getEffFromHisto(TH1D* h, double pt) {
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
