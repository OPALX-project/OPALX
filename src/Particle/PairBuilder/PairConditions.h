#ifndef PAIR_CONDITIONS_H
#define PAIR_CONDITIONS_H

#include "Region/NDRegion.h"

#include <limits>

template<class T, unsigned Dim>
class TrueCondition
{
public:
	TrueCondition() { }
	template<class V>
	bool operator()(const V &a, const V &b) const
	{
		return true;
	}

	T getRange(unsigned dim) const { return std::numeric_limits<T>::max(); }
};

template<class T, unsigned Dim>
class RadiusCondition
{
public:
 RadiusCondition(T r) : sqradius(r*r), radius(r)
	{  }

	template<class V>
	bool operator()(const V &a, const V &b) const
	{
		T sqr = 0;
		for(int d = 0;d<Dim;++d)
		{
			sqr += (a[d]-b[d])*(a[d]-b[d]);
		}
		return sqr <= sqradius;
	}

	T getRange(unsigned d) const { return radius; }
private:
	T sqradius, radius;
};

template<class T, unsigned Dim>
class BoxCondition
{
public:
	BoxCondition(T b[Dim]) : box(b)
	{  }

	template<class V>
	bool operator()(const V &a, const V &b) const
	{
		for(int d = 0;d<Dim;++d)
		{
			T diff = a[d]-b[d];
			if(diff > box[d] || diff < -box[d])
				return false;
		}
		return true;
	}

	T getRange(unsigned d) const { return box[d]; }
private:
	T box[Dim];
};

#endif
