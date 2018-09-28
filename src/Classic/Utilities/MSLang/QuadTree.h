#ifndef MSLANG_QUADTREE_H
#define MSLANG_QUADTREE_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct QuadTree {
        int level_m;
        std::list<Base*> objects_m;
        BoundingBox bb_m;
        QuadTree *nodes_m;

        QuadTree():
            level_m(0),
            bb_m(),
            nodes_m(0)
        { }

        QuadTree(int l, const BoundingBox &b):
            level_m(l),
            bb_m(b),
            nodes_m(0)
        { }

        QuadTree(const QuadTree &right);

        ~QuadTree();

        void operator=(const QuadTree &right);

        void transferIfInside(std::list<Base*> &objs);
        void buildUp();

        void writeGnuplot(std::ofstream &out) const;

        bool isInside(const Vector_t &R) const;
    };
}

#endif