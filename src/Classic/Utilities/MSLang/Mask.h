#ifndef MSLANG_MASK_H
#define MSLANG_MASK_H

#include "Utilities/MSLang.h"
#include "Utilities/MSLang/Rectangle.h"

namespace mslang {
    struct Mask: public Function {
        static bool parse_detail(iterator &it, const iterator &end, Function* &fun);
        virtual void print(int ident);
        virtual void apply(std::vector<Base*> &bfuncs);

        std::vector<Rectangle> pixels_m;
    };
}

#endif