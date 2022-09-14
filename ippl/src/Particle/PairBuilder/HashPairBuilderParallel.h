// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * HashPairBuilderParallel follows the Hockney and Eastwood approach to efficiently
 * find particle pairs. In this version of the code a local chaining mesh per processor 
 * is used to avoid looping empty buckets.
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/



#ifndef HASH_PAIR_BUILDER_PARALLEL_H
#define HASH_PAIR_BUILDER_PARALLEL_H

#include <algorithm>
#include <array>
#include <limits>
#include <cmath>
#include <set>
#include "Message/Communicate.h"

using namespace std;

template<class PBase>
class HashPairBuilderParallel
{
public:
    enum { Dim = PBase::Dim };
    typedef typename PBase::Position_t      Position_t;

    HashPairBuilderParallel(PBase& p_r, double gammaz) : particles_mr(p_r), gammaz_m(gammaz) 
    { hr_m = p_r.get_hr(); }

    template<class Pred, class OP>
    void forEach(const Pred& pred_r, const OP& op_r)
    {
        constexpr size_t END = numeric_limits<size_t>::max();
        size_t size = particles_mr.getLocalNum()+particles_mr.getGhostNum();

        NDIndex<3> locDomain = particles_mr.getFieldLayout().getLocalNDIndex();

        rmin_m = particles_mr.getMesh().get_origin();
        //To move to the boosted frame
        rmin_m[2] *= gammaz_m;
        hr_m[2] *= gammaz_m;

        //compute local chaining mesh
        Vektor<double,3> extendLLocal, extendRLocal, domainWidthLocal;
        for (unsigned i=0; i<3; ++i) {
            extendLLocal[i] = locDomain[i].first()*hr_m[i]+rmin_m[i];
            extendRLocal[i] = rmin_m[i]+(locDomain[i].last()+1)*hr_m[i];
            domainWidthLocal[i] = extendRLocal[i]-extendLLocal[i];
        
            //make sure that the chaining mesh covers the whole domain 
            //and has a gridwidth > r_cut
            bucketsPerDim_m[i]=floor(domainWidthLocal[i]/pred_r.getRange(i));
        
            if(bucketsPerDim_m[i] == 0) {
                bucketsPerDim_m[i] = 1;
            }

            hChaining_m[i] = domainWidthLocal[i]/bucketsPerDim_m[i];
        }


        //extend the chaining mesh by one layer of chaining cells in each dimension
        rmin_m = extendLLocal-hChaining_m;
        rmax_m = extendRLocal+hChaining_m;
        bucketsPerDim_m+=2;

        size_t Nbucket = bucketsPerDim_m[0]*bucketsPerDim_m[1]*bucketsPerDim_m[2];

        //index of first particle in this bucket
        vector<size_t> buckets(Nbucket);
        //index of next particle in this bucket. END indicates last particle of bucket
        vector<size_t> next(size);
        
        fill(buckets.begin(), buckets.end(), END);
        fill(next.begin(), next.end(), END);

        //in 3D we interact with 14 neighboring cells (including self cell interaction)
        unsigned neigh = 14;

        int offset[14][3] = {{ 1, 1, 1}, { 0, 1, 1}, {-1, 1, 1},
            { 1, 0, 1}, { 0, 0, 1}, {-1, 0, 1},
            { 1,-1, 1}, { 0,-1, 1}, {-1,-1, 1},
            { 1, 1, 0}, { 0, 1, 0}, {-1, 1, 0},
            { 1, 0, 0}, { 0, 0, 0}};

        //assign all particles to a bucket
        for(size_t i = 0;i<size;++i)
        {
            size_t bucketId = getBucketId(i);
            if(bucketId >= Nbucket) {
                cout << "Bucket with id: " << bucketId << " is wrong" << endl;
                cout << "Rank: " << Ippl::myNode() << endl;
                cout << "Buckets: " << bucketsPerDim_m << endl;
                cout << "Particle coords: " << particles_mr.R[i] << endl; 
                cout << "rmin_m: " << rmin_m << "rmax_m: " << rmax_m << endl;
            }
            next[i] = buckets[bucketId];
            buckets[bucketId] = i;
        }

        //loop over all buckets
        for (int bx=0; bx<bucketsPerDim_m[0]; ++bx) {
            for (int by=0; by<bucketsPerDim_m[1]; ++by) {
                for (int bz=0; bz<bucketsPerDim_m[2]; ++bz) {
                    unsigned bucketIdSelf = bz*bucketsPerDim_m[1]*bucketsPerDim_m[0]
                                              +by*bucketsPerDim_m[0]+bx;
                    //compute index of neighboring bucket to interact with
                    for (unsigned n=0; n<neigh;++n){
                        int bxNeigh, byNeigh, bzNeigh;

                        bxNeigh = bx+offset[n][0];
                        byNeigh = by+offset[n][1];
                        bzNeigh = bz+offset[n][2];

                        if (bxNeigh >= 0 && bxNeigh<bucketsPerDim_m[0] &&
                            byNeigh >= 0 && byNeigh<bucketsPerDim_m[1] &&
                            bzNeigh >= 0 && bzNeigh<bucketsPerDim_m[2]) {

                            unsigned bucketIdNeigh =
                            bzNeigh*bucketsPerDim_m[1]*bucketsPerDim_m[0]
                            +byNeigh*bucketsPerDim_m[0]+bxNeigh;

                            //i is index of particle considered in active chaining cell, 
                            //j is index of neighbor particle considered
                            size_t i = buckets[bucketIdSelf];
                            size_t j;

                            //loop over all particles in self cell
                            //self offset avoids double counting in self cell
                            int selfOffset = 0;
                            while (i != END) {
                                j = buckets[bucketIdNeigh];
                                //increase offset by number of processed particles in self cell
                                for (int o=0;o<selfOffset;o++){
                                    j = next[j];
                                }
                                //loop over all particles in neighbor cell
                                while(j != END) {
                                    if(pred_r(particles_mr.R[i], particles_mr.R[j])) {
                                        if (i!=j) {
                                            op_r(i, j, particles_mr);
                                        }
                                    }
                                    j = next[j];
                                }
                                i = next[i];
                                //adjust selfOffset
                                if (bucketIdSelf==bucketIdNeigh) {
                                    selfOffset++;
                                }
                                else {
                                    selfOffset=0;
                                }
                            }
                        }
                    }

                }
            }
        }

        
    }
private:

    //returns the bucket id of particle i
    size_t getBucketId(size_t i)
    {

        Vektor<int,3> loc;
        bool isInside, isOutsideMin, isOutsideMax;
        int indInside;
        for (unsigned d=0; d<3; ++d) {
            indInside = (particles_mr.R[i][d]-rmin_m[d])/hChaining_m[d];
            isInside = (particles_mr.R[i][d] > rmin_m[d]) && (particles_mr.R[i][d] < rmax_m[d]);
            isOutsideMin = (particles_mr.R[i][d] <= rmin_m[d]);
            isOutsideMax = (particles_mr.R[i][d] >= rmax_m[d]);
        
            //if the particle is inside the bucket take the inside index otherwise assign it to either first or
            //last bucket in that dimension
            loc[d] = ((int)isInside * indInside) + ((int)isOutsideMin * 0) + ((int)isOutsideMax * (bucketsPerDim_m[d]-1));
        }
        
        size_t bucketId = loc[2]*bucketsPerDim_m[1]*bucketsPerDim_m[0]+loc[1]*bucketsPerDim_m[0]+loc[0];
        return bucketId;
    }

    PBase& particles_mr;
    double gammaz_m;
    Vektor<int,3> bucketsPerDim_m;
    Vektor<double,3> hChaining_m;
    Vektor<double,3> rmin_m;
    Vektor<double,3> rmax_m;
    Vektor<double,3> hr_m;
};


#endif
