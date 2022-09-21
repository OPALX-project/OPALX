// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/
 
#ifndef HASH_PAIR_BUILDER_PARALLEL_H
#define HASH_PAIR_BUILDER_PARALLEL_H

/*
 * HashPairBuilderParallel - class for finding
 * particle pairs in the P3M solver and 
 * compute particle-particle interactions 
 * between them
 *  
 * This class follows the Hockney and Eastwood approach 
 * to efficiently find particle pairs. In this version 
 * of the code a local chaining mesh per processor 
 * is used to avoid looping empty buckets.
 */


#include "Message/Communicate.h"
#include "Utility/IpplException.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>

template<class PBase>
class HashPairBuilderParallel
{
public:
    enum { Dim = PBase::Dim };
    typedef typename PBase::Position_t      Position_t;

    HashPairBuilderParallel(PBase& p_r, double gammaz) : particles_mr(p_r), 
                                                         gammaz_m(gammaz) 
    { hr_m = p_r.get_hr(); }

    template<class Pred, class OP>
    void forEach(const Pred& pred_r, const OP& op_r)
    {
        constexpr std::size_t END = std::numeric_limits<std::size_t>::max();
        std::size_t size = particles_mr.getLocalNum()+particles_mr.getGhostNum();

        NDIndex<3> locDomain = particles_mr.getFieldLayout().getLocalNDIndex();

        rmin_m = particles_mr.getMesh().get_origin();
        //To move to the boosted frame
        rmin_m[2] *= gammaz_m;
        hr_m[2] *= gammaz_m;

        //compute local chaining mesh
        Vektor<double,3> extentLLocal, extentRLocal, domainWidthLocal;
        for (unsigned i=0; i<3; ++i) {
            extentLLocal[i] = locDomain[i].first()*hr_m[i]+rmin_m[i];
            extentRLocal[i] = rmin_m[i]+(locDomain[i].last()+1)*hr_m[i];
            domainWidthLocal[i] = extentRLocal[i]-extentLLocal[i];
        
            //make sure that the chaining mesh covers the whole domain 
            //and has a gridwidth > r_cut
            bucketsPerDim_m[i] = std::floor(domainWidthLocal[i]/pred_r.getRange(i));
        
            if(bucketsPerDim_m[i] == 0) {
                bucketsPerDim_m[i] = 1;
            }

            hChaining_m[i] = domainWidthLocal[i]/bucketsPerDim_m[i];
        }


        //extend the chaining mesh by one layer of chaining cells in each dimension
        rmin_m = extentLLocal-hChaining_m;
        rmax_m = extentRLocal+hChaining_m;
        bucketsPerDim_m+=2;

        std::size_t Nbucket = bucketsPerDim_m[0]*bucketsPerDim_m[1]*bucketsPerDim_m[2];

        //index of first particle in this bucket
        std::vector<std::size_t> buckets(Nbucket);
        //index of next particle in this bucket. END indicates last particle of bucket
        std::vector<std::size_t> next(size);
        
        std::fill(buckets.begin(), buckets.end(), END);
        std::fill(next.begin(), next.end(), END);

        //As per Hockney and Eastwood ``Computer simulation using particles" (section 8.4.3) 
        //we use Newton's third law 
        //in the force calculation and hence out of the total 27 
        //interactions in 3D we interact only 
        //with 13 neighboring cells + 1 self cell interaction
        unsigned neigh = 14;

        int offset[14][3] = {{ 1, 1, 1}, { 0, 1, 1}, {-1, 1, 1},
            { 1, 0, 1}, { 0, 0, 1}, {-1, 0, 1},
            { 1,-1, 1}, { 0,-1, 1}, {-1,-1, 1},
            { 1, 1, 0}, { 0, 1, 0}, {-1, 1, 0},
            { 1, 0, 0}, { 0, 0, 0}};

        //assign all particles to a bucket
        for(std::size_t i = 0;i<size;++i)
        {
            std::size_t bucketId = getBucketId(i);
            if(bucketId >= Nbucket) {
                std::cout << "Bucket with id: " << bucketId << " is wrong" << std::endl;
                std::cout << "Rank: " << Ippl::myNode() << std::endl;
                std::cout << "Buckets: " << bucketsPerDim_m << std::endl;
                std::cout << "Particle coords: " << particles_mr.R[i] << std::endl; 
                std::cout << "rmin_m: " << rmin_m << "rmax_m: " << rmax_m << std::endl;
                throw IpplException("HashPairBuilderParallel::forEach", 
                            "Particle outside the local domain");
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
                            std::size_t i = buckets[bucketIdSelf];
                            std::size_t j;

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
    std::size_t getBucketId(std::size_t i)
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
            //last bucket in that dimension. (int)isOutsideMin * 0 is written explicitly for the ease of code
            //understanding.
            loc[d] = ((int)isInside * indInside) + ((int)isOutsideMin * 0) + ((int)isOutsideMax * (bucketsPerDim_m[d]-1));
        }
        
        std::size_t bucketId = loc[2]*bucketsPerDim_m[1]*bucketsPerDim_m[0]+loc[1]*bucketsPerDim_m[0]+loc[0];
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
