#ifndef MSLANG_ROTATION_H
#define MSLANG_ROTATION_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct Rotation: public Function {
        Function* func_m;
        double angle_m;

        virtual ~Rotation() {
            delete func_m;
        }

        virtual void print(int indentwidth);
        void applyRotation(std::vector<Base*> &bfuncs);
        virtual void apply(std::vector<Base*> &bfuncs);
        static bool parse_detail(iterator &it, const iterator &end, Function* &fun);
    };
}

#endif
