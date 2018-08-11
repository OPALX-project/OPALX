#include "Utilities/MSLang/QuadTree.h"

namespace mslang {
    QuadTree::QuadTree(const QuadTree &right):
        level_m(right.level_m),
        objects_m(right.objects_m.begin(),
                  right.objects_m.end()),
        bb_m(right.bb_m),
        nodes_m(0)
    {
        if (right.nodes_m != 0) {
            nodes_m = new QuadTree[4];
            for (unsigned int i = 0; i < 4u; ++ i) {
                nodes_m[i] = right.nodes_m[i];
            }
        }
    }

    QuadTree::~QuadTree() {
        for (Base *&obj: objects_m)
            obj = 0; // memory isn't handled by QuadTree class

        if (nodes_m != 0) {
            delete[] nodes_m;
        }
        nodes_m = 0;
    }

    void QuadTree::operator=(const QuadTree &right) {
        level_m = right.level_m;
        objects_m.insert(objects_m.end(),
                         right.objects_m.begin(),
                         right.objects_m.end());
        bb_m = right.bb_m;

        if (nodes_m != 0) delete[] nodes_m;
        nodes_m = 0;

        if (right.nodes_m != 0) {
            nodes_m = new QuadTree[4];
            for (unsigned int i = 0; i < 4u; ++ i) {
                nodes_m[i] = right.nodes_m[i];
            }
        }
    }

    void QuadTree::transferIfInside(std::list<Base*> &objs) {
        for (Base* &obj: objs) {
            if (bb_m.isInside(obj->bb_m)) {
                objects_m.push_back(obj);
                obj = 0;
            }
        }

        objs.remove_if([](const Base *obj) { return obj == 0; });
    }

    void QuadTree::buildUp() {
        QuadTree *next = new QuadTree[4];
        next[0] = QuadTree(level_m + 1,
                           BoundingBox(bb_m.center_m,
                                       Vector_t(bb_m.center_m[0] + 0.5 * bb_m.width_m,
                                                bb_m.center_m[1] + 0.5 * bb_m.height_m,
                                                0.0)));
        next[1] = QuadTree(level_m + 1,
                           BoundingBox(Vector_t(bb_m.center_m[0],
                                                bb_m.center_m[1] - 0.5 * bb_m.height_m,
                                                0.0),
                                       Vector_t(bb_m.center_m[0] + 0.5 * bb_m.width_m,
                                                bb_m.center_m[1],
                                                0.0)));
        next[2] = QuadTree(level_m + 1,
                           BoundingBox(Vector_t(bb_m.center_m[0] - 0.5 * bb_m.width_m,
                                                bb_m.center_m[1],
                                                0.0),
                                       Vector_t(bb_m.center_m[0],
                                                bb_m.center_m[1] + 0.5 * bb_m.height_m,
                                                0.0)));
        next[3] = QuadTree(level_m + 1,
                           BoundingBox(Vector_t(bb_m.center_m[0] - 0.5 * bb_m.width_m,
                                                bb_m.center_m[1] - 0.5 * bb_m.height_m,
                                                0.0),
                                       bb_m.center_m));

        for (unsigned int i = 0; i < 4u; ++ i) {
            next[i].transferIfInside(objects_m);
        }

        bool allEmpty = true;
        for (unsigned int i = 0; i < 4u; ++ i) {
            if (next[i].objects_m.size() != 0) {
                allEmpty = false;
                break;
            }
        }

        if (allEmpty) {
            for (unsigned int i = 0; i < 4u; ++ i) {
                objects_m.merge(next[i].objects_m);
            }

            delete[] next;
            return;
        }

        for (unsigned int i = 0; i < 4u; ++ i) {
            next[i].buildUp();
        }

        nodes_m = next;
    }

    void QuadTree::writeGnuplot(std::ofstream &out) const {
        out << "# level: " << level_m << ", size: " << objects_m.size() << std::endl;
        bb_m.writeGnuplot(out);
        out << "# num holes: " << objects_m.size() << std::endl;
        // for (const Base *obj: objects_m) {
        //     obj->writeGnuplot(out);
        // }
        out << std::endl;

        if (nodes_m != 0) {
            for (unsigned int i = 0; i < 4u; ++ i) {
                nodes_m[i].writeGnuplot(out);
            }
        }
    }

    bool QuadTree::isInside(const Vector_t &R) const {
        if (nodes_m != 0) {
            Vector_t X = R - bb_m.center_m;
            unsigned int idx = (X[1] >= 0.0 ? 0: 1);
            idx += (X[0] >= 0.0 ? 0: 2);

            if (nodes_m[idx].isInside(R)) {
                return true;
            }
        }

        for (Base* obj: objects_m) {
            if (obj->isInside(R)) {
                return true;
            }
        }

        return false;
    }
}