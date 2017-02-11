#ifndef OPAL_INDEXMAP_H
#define OPAL_INDEXMAP_H

#include <ostream>

#include "AbsBeamline/Component.h"
#include "Utilities/OpalException.h"

#include <map>
#include <set>
#include <utility>


class IndexMap
{
public:
    typedef std::pair<double, double> key_t;
    typedef std::set<std::shared_ptr<Component> > value_t;

    IndexMap();

    void add(key_t::first_type initialStep, key_t::second_type finalStep, const value_t &val);

    value_t query(key_t::first_type s, key_t::second_type ds);

    void tidyUp();

    void print(std::ostream&) const;
    void saveSDDS(double startS) const;
    size_t size() const;

    class OutOfBounds: public OpalException {
    public:
        OutOfBounds(const std::string &meth, const std::string &msg):
            OpalException(meth, msg) { }

        OutOfBounds(const OutOfBounds &rhs):
            OpalException(rhs) { }

        virtual ~OutOfBounds() { }

    private:
        OutOfBounds();
    };

private:
    class myCompare {
    public:
        bool operator()(const key_t x , const key_t y)
        {
            if (x.first < y.first) return true;

            if (x.first == y.first) {
                if (x.second < y.second) return true;
            }

            return false;
        }
    };

    typedef std::map<key_t, value_t, myCompare> map_t;
    map_t map_m;
    double totalPathLength_m;
    const double oneMinusEpsilon_m;
};

inline
void IndexMap::add(key_t::first_type initialS, key_t::second_type finalS, const value_t &val) {
    key_t key(initialS, finalS * oneMinusEpsilon_m);
    map_m.insert(std::pair<key_t, value_t>(key, val));
    totalPathLength_m = (*map_m.rbegin()).first.second;
}

inline
size_t IndexMap::size() const {
    return map_m.size();
}

inline
std::ostream& operator<< (std::ostream &out, const IndexMap &im)
{
    im.print(out);
    return out;
}

inline
Inform& operator<< (Inform &out, const IndexMap &im) {
    im.print(out.getStream());
    return out;
}

#endif