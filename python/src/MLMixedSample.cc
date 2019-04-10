#include "MLMixedSample.h"

#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <TIterator.h>
#include <TMath.h>

using namespace RooFit;

#include <algorithm>
#include <iostream>

ClassImp(MLMixedSample);

///Create a mixed sample from the two original samples
RooDataSet* MLMixedSample::createMixedSample(RooDataSet* data, RooDataSet* gen, RooArgSet vars) const{

        RooCategory cat("cat","cat");
    cat.defineType("A",0);
    cat.defineType("B",1);

    RooDataSet* mixed = new RooDataSet("mixed","mixed",vars,Index(cat),Import("A",*data),Import("B",*gen));
    return mixed;

}

void MLMixedSample::fillIndex(const RooDataSet& mixed, RooArgSet vars){

        //start off by getting the RMS values
        TIterator* it = vars.createIterator();
        TObject* obj = 0;
        while( (obj = it->Next()) ){
                RooRealVar* v = dynamic_cast<RooRealVar*>(obj);
                if(v){
                        (*rms_)[v->GetName()] = mixed.rmsVar(*v)->getVal();
                }
        }
        delete it;

        const unsigned int entries = mixed.numEntries();
        const unsigned int nVars = rms_->size();

        //this will contain the catchable parts of the nearest neighbor calculation - i.e. ignoring the cross-terms
        index_->resize(entries);

        //fill the data into a vector
        cat_->resize(entries);
        data_->resize(entries*nVars);

        for(unsigned int i = 0; i < entries; i++){
                const RooArgSet* row =  mixed.get(i);

                const int cat = row->getCatIndex("cat");
                cat_->at(i) = cat;

                unsigned int nameIndex = 0;
                for(rms_map::iterator name = rms_->begin(); name != rms_->end(); ++name){
                        const double val = row->getRealValue(name->first.c_str());
                        const double rms = name->second;
                        //the minimum is when \sum^D x_i/w^2 == \sum^D x_i/w^2
                        if(nameIndex == 0){
                                index_->at(i) = std::make_pair( (val/rms*rms), i);
                        }else{
                                const neighbour& n = index_->at(i);
                                index_->at(i) = std::make_pair( n.first + (val/rms*rms), i);
                        }
                        //cache the dataser to speed up access
                        data_->at(nameIndex*entries + i) = val;

                        nameIndex++;
                }
        }
        //sort by the value of \sum^D x_i/w^2
        std::sort(index_->begin(), index_->end());

}

double MLMixedSample::testStatistic(RooDataSet* data, RooDataSet* gen, RooArgSet vars){

    const RooDataSet* mixed = createMixedSample(data,gen,vars);
    fillIndex(*mixed,vars);

    const double n_a = data->numEntries();
    const double n_b = gen->numEntries();
    const double n = n_a + n_b;
    const double t = sum_I(*mixed,vars)/(n_k_*n);

    const double mu_t = (n_a*(n_a - 1) + n_b*(n_b - 1))/(n*(n-1));
    const double sigmaT = TMath::Sqrt( ( ((n_a*n_b)/(n*n)) + ((4*n_a*n_a*n_b*n_b)/(n*n*n*n)) )/(n*n_k_) );

    delete mixed;
    clear();

    return (t - mu_t)/sigmaT;
}

double MLMixedSample::confidence(const double x, const double mu, const double sigma){
        //this is not correct yet. I need to dig further
        return (TMath::Sqrt(TMath::Pi()/2.)*(1 + TMath::Erf(((x - mu)*TMath::Sqrt(TMath::Power(sigma,-2)))/
                        TMath::Sqrt(2))))/TMath::Sqrt(TMath::Power(sigma,-2));
}

std::auto_ptr<MLMixedSample::neighbours> MLMixedSample::nearest_neighbour(const RooDataSet& mixed, RooArgSet vars, const unsigned int i) const{
        std::cout << "Warning: variable vars (" << vars << ") is ignored in the nearest neighbor calculation." << std::endl;

        const unsigned int entries = mixed.numEntries();

        //double sum = 0;
        unsigned int nameIndex = 0;

        //calculate the distance metric used in the index
        double thisPoint = 0;
        for(rms_map::iterator it = rms_->begin(); it != rms_->end(); ++it){
        const double val = data_->at(nameIndex*entries + i);
        const double rms = it->second;
        thisPoint += (val/rms*rms);
        nameIndex++;
        }

        //look at the closest points in our sorted list of points
        const neighbour n = std::make_pair(thisPoint,i);
        neighbours::iterator start = index_->begin();
        neighbours::iterator finish = index_->end();
        if(diff_factor_ >= 1){
                const unsigned int diff = n_k_ * diff_factor_;
                start = std::lower_bound(index_->begin(),index_->end(),n) - diff;
                finish = std::upper_bound(index_->begin(),index_->end(),n) + diff;
                if(start < index_->begin()){
                        start = index_->begin();
                }
                if(finish > index_->end()){
                        finish = index_->end();
                }
        }

        //this final bit is O(N^2) but at least we've cut down N a lot
        const int cat_i = cat_->at(i);
        MLMixedSample::neighbours nearest;
        nearest.reserve(mixed.numEntries());

        for(neighbours::iterator it = start; it != finish; ++it){
                const unsigned int j = it->second;
                //skip this event
                if(i == j){
                        continue;
                }
                const int cat_j = cat_->at(j);
                //this loop is eqn 3.3
                double sum = 0;
                unsigned int nameIndex = 0;
                for(rms_map::iterator it = rms_->begin(); it != rms_->end(); ++it){

                    const double x_i = data_->at(nameIndex*entries + i);
                    const double x_j = data_->at(nameIndex*entries + j);

                    const double nn = (x_i - x_j)/it->second;
                    sum += (nn*nn);

                    nameIndex++;
                }
                const unsigned int same = (cat_i == cat_j) ? 1 : 0;
                nearest.push_back(std::make_pair(sum,same));
        }
        //sort to get the n_k nearest neighbors
        std::sort(nearest.begin(),nearest.end());//can be O(N*N) or better O(N Log(N))

        MLMixedSample::neighbours* result = new MLMixedSample::neighbours();
        result->reserve(n_k_);
        for(unsigned int i = 0; i < n_k_; i++){
                result->push_back(nearest.at(i));
        }
        return std::auto_ptr<MLMixedSample::neighbours>(result);

}

double MLMixedSample::sum_I(const RooDataSet& mixed, RooArgSet vars) const{

        const unsigned int entries = mixed.numEntries();
        const unsigned int printEvery = TMath::Nint(0.01*entries);

    double sum = 0;
    //has complexity O(N*N) - Rather slow
    for(unsigned int i = 0; i < entries; i++){
        //has complexity O(N)
        std::auto_ptr<MLMixedSample::neighbours> nn = MLMixedSample::nearest_neighbour(mixed, vars, i);
        for(MLMixedSample::neighbours::iterator n = nn->begin(); n != nn->end(); ++n){
                sum += n->second;
        }
        if( (i % printEvery) == 0){
                std::cout << "Evaluated events: " << i << "/" << entries << " (" << i/(1.*entries) << ")" << std::endl;
        }
    }
    return sum;

}
