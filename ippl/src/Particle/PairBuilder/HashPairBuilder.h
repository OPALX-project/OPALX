#ifndef HASH_PAIR_BUILDER_H
#define HASH_PAIR_BUILDER_H

#include <algorithm>
#include <limits>
#include <cmath>
#include <set>

template<class PBase>
class HashPairBuilder
{
public:
	enum { Dim = PBase::Dim };  
	typedef typename PBase::Position_t      Position_t;
	
	HashPairBuilder(PBase &p) : particles(p) { }
	
	template<class Pred, class OP>
	void for_each(const Pred& pred, const OP &op)
	{
		const std::size_t END = std::numeric_limits<std::size_t>::max();
		int f[3];
		
		std::size_t size = particles.getLocalNum()+particles.getGhostNum();
		
		int edge = std::pow(size/8, 1./Dim);
		
		std::size_t Nbucket = 1;
		for(int d = 0;d<Dim;++d)
		{
		        f[d] = edge;
			Nbucket *= f[d];
			if(pred.getRange(d)==0)
			  return;
		}
		
		std::size_t *buckets = new size_t[Nbucket];
		std::size_t *next = new size_t[size];
		std::fill(buckets, buckets+Nbucket, END);
		std::fill(next, next+size, END);
		
		int neigh = 1;
		for(int d = 0;d<Dim;++d)
			neigh *= 3;
		neigh /= 2;
		
		int offset[14][3] = {{ 1, 1, 1}, { 0, 1, 1}, {-1, 1, 1},
							 { 1, 0, 1}, { 0, 0, 1}, {-1, 0, 1},
							 { 1,-1, 1}, { 0,-1, 1}, {-1,-1, 1},
							 { 1, 1, 0}, { 0, 1, 0}, {-1, 1, 0},
							 { 1, 0, 0}, { 0, 0, 0}};
		
		for(std::size_t i = 0;i<size;++i)
		{
			std::size_t pos = sum(i, pred, f, offset[13]);
			next[i] = buckets[pos];
			buckets[pos] = i;
		}	
		
		for(std::size_t i = 0;i<size;++i)
		{	
			std::size_t j = next[i];
			while(j != END)
			{
				if(pred(particles.R[i], particles.R[j]))
					op(i, j, particles);
				j = next[j];
			}
			
			for(int k = 0;k<neigh;++k)
			{
				std::size_t tmppos = sum(i, pred, f, offset[k]);
				
				j = buckets[tmppos];
				while(j != END)
				{
					if(pred(particles.R[i], particles.R[j]))
						op(i, j, particles);
					j = next[j];
				}
			}
			
		}
		delete[] buckets;
                delete[] next;
	}
private:

	template<class Pred>
	int sum(int i, const Pred& pred, int f[], int offset[])
	{
		int sum = 0;
		for(int d = 0;d<Dim;++d)
		{
			double scaled = particles.R[i][d]/pred.getRange(d);
			int pos = mod(int(floor(scaled+offset[d])), f[d]);
			for(int dd = 0;dd<d;++dd)
				pos*=f[dd];
			sum += pos;
		}
		return sum;
	}
	
	int mod(int x, int m)
	{
		if(x>=0)
			return x%m;
		else
			return (m - ((-x)%m))%m;
	}

	PBase &particles;
};


#endif
