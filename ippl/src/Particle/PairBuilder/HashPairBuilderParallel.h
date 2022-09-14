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

    HashPairBuilderParallel(PBase &p, double gammaz_) : particles(p), gammaz(gammaz_) 
    { hr_m = p.get_hr(); }

    template<class Pred, class OP>
    void for_each(const Pred& pred, const OP &op)
    {
        const size_t END = numeric_limits<size_t>::max();
        size_t size = particles.getLocalNum()+particles.getGhostNum();

        NDIndex<3> locDomain = particles.getFieldLayout().getLocalNDIndex();

        rmin_m = particles.getMesh().get_origin();
        //To move to the boosted frame
        rmin_m[2] *= gammaz;
        hr_m[2] *= gammaz;

        //compute local chaining mesh
        Vektor<double,3> extend_l_local, extend_r_local, domain_width_local;
        for (unsigned i=0; i<3; ++i) {
            extend_l_local[i] = locDomain[i].first()*hr_m[i]+rmin_m[i];
            extend_r_local[i] = rmin_m[i]+(locDomain[i].last()+1)*hr_m[i];
            domain_width_local[i] = extend_r_local[i]-extend_l_local[i];
        
            //make sure that the chaining mesh covers the whole domain 
            //and has a gridwidth > r_cut
            buckets_per_dim[i]=floor(domain_width_local[i]/pred.getRange(i));
        
            if(buckets_per_dim[i] == 0)
                buckets_per_dim[i] = 1;

            h_chaining[i] = domain_width_local[i]/buckets_per_dim[i];
        }


        //extend the chaining mesh by one layer of chaining cells in each dimension
        rmin_m = extend_l_local-h_chaining;
        rmax_m = extend_r_local+h_chaining;
        buckets_per_dim+=2;

        size_t Nbucket = buckets_per_dim[0]*buckets_per_dim[1]*buckets_per_dim[2];

        //index of first particle in this bucket
        size_t *buckets = new size_t[Nbucket]; 
        //index of next particle in this bucket. END indicates last particle of bucket
        size_t *next = new size_t[size]; 
        fill(buckets, buckets+Nbucket, END);
        fill(next, next+size, END);

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
            size_t bucket_id = get_bucket_id(i);
            if(bucket_id >= Nbucket) {
                cout << "Bucket with id: " << bucket_id << " is wrong" << endl;
                cout << "Rank: " << Ippl::myNode() << endl;
                cout << "Buckets: " << buckets_per_dim << endl;
                cout << "Particle coords: " << particles.R[i] << endl; 
                cout << "rmin_m: " << rmin_m << "rmax_m: " << rmax_m << endl;
            }
            next[i] = buckets[bucket_id];
            buckets[bucket_id] = i;
        }

        //loop over all buckets
        for (int bx=0; bx<buckets_per_dim[0]; ++bx) {
            for (int by=0; by<buckets_per_dim[1]; ++by) {
                for (int bz=0; bz<buckets_per_dim[2]; ++bz) {
                    unsigned bucket_id_self = bz*buckets_per_dim[1]*buckets_per_dim[0]
                                              +by*buckets_per_dim[0]+bx;
                    //compute index of neighboring bucket to interact with
                    for (unsigned n=0; n<neigh;++n){
                        int bx_neigh, by_neigh, bz_neigh;

                        bx_neigh = bx+offset[n][0];
                        by_neigh = by+offset[n][1];
                        bz_neigh = bz+offset[n][2];

                        if (bx_neigh >= 0 && bx_neigh<buckets_per_dim[0] &&
                            by_neigh >= 0 && by_neigh<buckets_per_dim[1] &&
                            bz_neigh >= 0 && bz_neigh<buckets_per_dim[2]) {

                            unsigned bucket_id_neigh =
                            bz_neigh*buckets_per_dim[1]*buckets_per_dim[0]+by_neigh*buckets_per_dim[0]+bx_neigh;

                            //i is index of particle considered in active chaining cell, 
                            //j is index of neighbor particle considered
                            size_t i = buckets[bucket_id_self];
                            size_t j;

                            //loop over all particles in self cell
                            //self offset avoids double counting in self cell
                            int self_offset = 0;
                            while (i != END) {
                                j = buckets[bucket_id_neigh];
                                //increase offset by number of processed particles in self cell
                                for (int o=0;o<self_offset;o++){
                                    j = next[j];
                                }
                                //loop over all particles in neighbor cell
                                while(j != END) {
                                    if(pred(particles.R[i], particles.R[j])) {
                                        if (i!=j)
                                            op(i, j, particles);
                                    }
                                    j = next[j];
                                }
                                i = next[i];
                                //adjust self_offset
                                if (bucket_id_self==bucket_id_neigh)
                                    self_offset++;
                                else
                                    self_offset=0;
                            }
                        }
                    }

                }
            }
        }

        
        delete[] buckets;
        delete[] next;
    }
private:

    //returns the bucket id of particle i
    size_t get_bucket_id(size_t i)
    {

        Vektor<int,3> loc;
        bool isInside, isOutsideMin, isOutsideMax;
        int indInside;
        for (unsigned d=0; d<3; ++d) {
            indInside = (particles.R[i][d]-rmin_m[d])/h_chaining[d];
            isInside = (particles.R[i][d] > rmin_m[d]) && (particles.R[i][d] < rmax_m[d]);
            isOutsideMin = (particles.R[i][d] <= rmin_m[d]);
            isOutsideMax = (particles.R[i][d] >= rmax_m[d]);
        
            //if the particle is inside the bucket take the inside index otherwise assign it to either first or
            //last bucket in that dimension
            loc[d] = ((int)isInside * indInside) + ((int)isOutsideMin * 0) + ((int)isOutsideMax * (buckets_per_dim[d]-1));
        }
        
        size_t bucket_id = loc[2]*buckets_per_dim[1]*buckets_per_dim[0]+loc[1]*buckets_per_dim[0]+loc[0];
        return bucket_id;
    }

    PBase &particles;
    double gammaz;
    Vektor<int,3> buckets_per_dim;
    Vektor<double,3> h_chaining;
    Vektor<double,3> rmin_m;
    Vektor<double,3> rmax_m;
    Vektor<double,3> hr_m;
};


#endif
