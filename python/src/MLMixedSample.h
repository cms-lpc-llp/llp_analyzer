#ifndef MLMIXEDSAMPLE_H
#define MLMIXEDSAMPLE_H

class RooArgSet;
class RooDataSet;

#include <TObject.h>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

///An implementation of the mixed sample method from arXiv:1006.3019 - section 3.2
class MLMixedSample : public TObject {

protected:
        typedef std::pair<double,unsigned int> neighbour;
        typedef std::vector<neighbour> neighbours;
        typedef std::map<std::string,double> rms_map;
        typedef std::vector<double> data_vector;
        typedef std::vector<unsigned int> cat_vector;
        
public:
        
        /**
         * @param n_k :                 The number of nearest neighbours to consider - analgous to a histogram bin size, but in D dimensions
         * @param diff_factor : The full nearest neighbour problem is O(N^2). We can cut this down by using some approximations
         *                      where only a subset of points in the D dimensional parameter space are considered. 
         *                      If this number is small - i.e. close to 1, you will not get the right answer. As diff_factor*n_k -> numEntries
         *                      the answer becomes less approximate but the evaluation slows down considerably. If in doubt, increase the factor
         *                      until the result becomes stable. The default value is large, but has been seen to produce reasonable results
         *                      in an amount of time that is tolerable. If you set this parameter less than 1, you get all events considered.
         *                      For smallish datasets, this is probably best.  
         */
        MLMixedSample(const unsigned int n_k = 10, const int diff_factor = 500):
                TObject::TObject(),
                n_k_(n_k),
                diff_factor_(diff_factor),
                index_(new neighbours),
                rms_(new rms_map),
                cat_(new cat_vector),
                data_(new data_vector){
        }

        RooDataSet* createMixedSample(RooDataSet* data, RooDataSet* gen, RooArgSet vars) const;
        double testStatistic(RooDataSet* data, RooDataSet* gen, RooArgSet vars);
        
        //convert the test statistic into a confidence level
        static double confidence(const double x, const double mu = 0, const double sigma = 1);
        
protected:

        std::auto_ptr<neighbours> nearest_neighbour(const RooDataSet& mixed, RooArgSet vars, const unsigned int i) const;
        double sum_I(const RooDataSet& mixed, RooArgSet vars) const;
        
        void fillIndex(const RooDataSet& mixed, RooArgSet vars);

        ClassDef(MLMixedSample,1);
        
private:

        const unsigned int n_k_;
        const int diff_factor_;
        ///The 
        std::auto_ptr<neighbours> index_;
        std::auto_ptr<rms_map> rms_;
        std::auto_ptr<cat_vector> cat_;
        std::auto_ptr<data_vector> data_;
        
        void clear(){
                index_->clear();
                rms_->clear();
                cat_->clear();
                data_->clear();
        }

};

#endif /*MLMIXEDSAMPLE_H*/
