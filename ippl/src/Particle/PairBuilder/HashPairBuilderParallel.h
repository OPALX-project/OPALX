// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * HashPairBuilderParallel follows the Hockney and Eastwood approach to efficiently
 * find particle pairs. In this version of the code a local Chaining mesh per processor is used to avoid looping
 * empty buckets.
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

template<class PBase>
class HashPairBuilderParallel
{
public:
    enum { Dim = PBase::Dim };
    typedef typename PBase::Position_t      Position_t;

    HashPairBuilderParallel(PBase &p, double gammaz_, long long timestep_) : particles(p), gammaz(gammaz_), timestep_m(timestep_) 
    { hr_m = p.get_hr(); }

    template<class Pred, class OP>
    void for_each(const Pred& pred, const OP &op)
    {
        const std::size_t END = std::numeric_limits<std::size_t>::max();
        std::size_t size = particles.getLocalNum()+particles.getGhostNum();

        Inform dmsg("debug_msg:");
        // dmsg << "We use parallel hash pair builder small chaining mesh ****************************" << endl;

        NDIndex<3> globDomain = particles.getFieldLayout().getDomain();
        NDIndex<3> locDomain = particles.getFieldLayout().getLocalNDIndex();

        rmin_m = particles.getMesh().get_origin();
        rmin_m[2] *= gammaz;
        hr_m[2] *= gammaz;

        ///////// compute local chaining mesh
        Vektor<double,3> extend_l_local, extend_r_local, domain_width_local;
        for (unsigned i=0; i<3; ++i) {
            extend_l_local[i] = locDomain[i].first()*hr_m[i]+rmin_m[i];
            extend_r_local[i] = rmin_m[i]+(locDomain[i].last()+1)*hr_m[i];
            domain_width_local[i] = extend_r_local[i]-extend_l_local[i];
        }

        //make sure that the chaining mesh covers the whole domain and has a gridwidth > r_cut
        buckets_per_dim[0]=floor(domain_width_local[0]/pred.getRange(0));
        buckets_per_dim[1]=floor(domain_width_local[1]/pred.getRange(1));
        buckets_per_dim[2]=floor(domain_width_local[2]/pred.getRange(2));
        //buckets_per_dim[0]=1;//floor(domain_width_local[0]/pred.getRange(0));
        //buckets_per_dim[1]=1;//floor(domain_width_local[1]/pred.getRange(1));
        //buckets_per_dim[2]=1;//floor(domain_width_local[2]/pred.getRange(2));

        for (unsigned dim = 0; dim<3; ++dim) {
            if(buckets_per_dim[dim] == 0)
                buckets_per_dim[dim] = 1;
        }
        

        dmsg << "local domain width = " << domain_width_local << endl;
        //dmsg << "local extends : " << extend_l_local << "\t" << extend_r_local << endl;
        //int factor = 1;
        for (unsigned dim = 0; dim<3; ++dim) {
            h_chaining[dim] = domain_width_local[dim]/buckets_per_dim[dim];
            //if(h_chaining[dim] <= pred.getRange(dim)) {
            //    factor = ceil(pred.getRange(dim)/h_chaining[dim]);
            //}
            //rmin_m[dim] = extend_l_local[dim] - factor*h_chaining[dim];
            //rmax_m[dim] = extend_r_local[dim] + factor*h_chaining[dim];
            //buckets_per_dim[dim] += 2*factor;
        }


        //extend the chaining mesh by one layer of chaining cells in each dimension
        rmin_m = extend_l_local-h_chaining;
        rmax_m = extend_r_local+h_chaining;
        buckets_per_dim+=2;

        //dmsg << "gammaz = " << gammaz << endl;
        //dmsg << "local domain width = " << domain_width_local << endl;
        //dmsg << "local extends with chaining: " << rmin_m << "\t" << rmax_m << endl;
        //dmsg << "h_chaining = " << h_chaining << endl;

        std::size_t Nbucket = buckets_per_dim[0]*buckets_per_dim[1]*buckets_per_dim[2];
        dmsg << "Nbuckets = " << Nbucket << endl;
        dmsg << "buckets = " << buckets_per_dim << endl;


        //std::vector<std::size_t> buckets(Nbucket);
        //std::vector<std::size_t> next(size);

        //std::fill(buckets.begin(), buckets.end(), END);
        //std::fill(next.begin(), next.end(), END);



        std::size_t *buckets = new size_t[Nbucket]; //index of first particle in this bucket
        std::size_t *next = new size_t[size]; //index of next particle in this bucket. END indicates last particle of bucket
        std::fill(buckets, buckets+Nbucket, END);
        std::fill(next, next+size, END);

        //in 3D we interact with 14 neighboring cells (including self cell interaction)
        unsigned neigh = 14;

        int offset[14][3] = {{ 1, 1, 1}, { 0, 1, 1}, {-1, 1, 1},
            { 1, 0, 1}, { 0, 0, 1}, {-1, 0, 1},
            { 1,-1, 1}, { 0,-1, 1}, {-1,-1, 1},
            { 1, 1, 0}, { 0, 1, 0}, {-1, 1, 0},
            { 1, 0, 0}, { 0, 0, 0}};

        //assign all particles to a bucket
        for(std::size_t i = 0;i<size;++i)
        {
            std::size_t bucket_id = get_bucket_id(i,pred);
            //dmsg << "we got bucket id = " << bucket_id << endl;
            if(bucket_id >= Nbucket) {
                std::cout << "Bucket with id: " << bucket_id << " is wrong" << std::endl;
                std::cout << "Rank: " << Ippl::myNode() << std::endl;
                std::cout << "Buckets: " << buckets_per_dim << std::endl;
                std::cout << "Particle coords: " << particles.R[i] << std::endl; 
                std::cout << "rmin_m: " << rmin_m << "rmax_m: " << rmax_m << std::endl;
            }
            next[i] = buckets[bucket_id];
            buckets[bucket_id] = i;
        }

        //dmsg << "Bucket id calculations finished " << endl;
        //std::cout << "rank: " << ippl::mynode() << "bucket id calculations finished " << std::endl;
        double part_count = 0;
        Vektor<double,3> shift(0,0,0);
        //loop over all buckets
        for (int bx=0; bx<buckets_per_dim[0]; ++bx) {
            for (int by=0; by<buckets_per_dim[1]; ++by) {
                for (int bz=0; bz<buckets_per_dim[2]; ++bz) {
                    unsigned bucket_id_self = bz*buckets_per_dim[1]*buckets_per_dim[0]+by*buckets_per_dim[0]+bx;
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

                            //i is index of particle considered in active cahining cell, j is index of neighbor particle considered
                            std::size_t i = buckets[bucket_id_self];
                            std::size_t j;

                            //loop over all particles in self cell
                            //self offset avoids double counting in self cell
                            int self_offset = 0;
                            part_count = 0;
                            while (i != END) {
                                part_count++;
                                j = buckets[bucket_id_neigh];
                                //increase offset by number of processed particles in self cell
                                for (int o=0;o<self_offset;o++){
                                    j = next[j];
                                }
                                //loop over all particles in neighbor cell
                                while(j != END) {
                                    if(pred(particles.R[i], particles.R[j])) {
                                        if (i!=j)
                                            op(i, j, particles, shift);
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
        
        //std::cout << "Rank: " << Ippl::myNode() << "PP calculations finished " << std::endl;

        //Ippl::Comm->barrier();
        delete[] buckets;
        delete[] next;
    }
private:

    //returns the bucket id of particle i
    template<class Pred>
    std::size_t get_bucket_id(std::size_t i, const Pred& /*pred*/)
    {
        //Inform dmsg("debug_msg:");

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
        
        std::size_t bucket_id = loc[2]*buckets_per_dim[1]*buckets_per_dim[0]+loc[1]*buckets_per_dim[0]+loc[0];


        //if(timestep_m == 8662) {
        //    std::cout << "Rank: " << Ippl::myNode() << "bucket id of particle " << i << "with coords " << particles.R[i] << " = [" << loc[0] << "," << loc[1] << "," << loc[2] << "] => bucket id = "  << bucket_id << std::endl;
        //}
        //dmsg << particles.R[i][0] << "," << particles.R[i][1] << "," << particles.R[i][2] << "," << bucket_id << endl;
        return bucket_id;
    }

    PBase &particles;
    double gammaz;
    long long timestep_m;
    Vektor<int,3> buckets_per_dim;
    Vektor<double,3> h_chaining;
    Vektor<double,3> rmin_m;
    Vektor<double,3> rmax_m;
    Vektor<double,3> hr_m;
};


#endif
