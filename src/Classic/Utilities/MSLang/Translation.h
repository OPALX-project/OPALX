#ifndef MSLANG_TRANSLATE_H
#define MSLANG_TRANSLATE_H

#include "Utilities/MSLang.h"

namespace mslang {
    struct Translation: public Function {
        Function* func_m;
        double shiftx_m;
        double shifty_m;

        virtual ~Translation() {
            delete func_m;
        }

        virtual void print(int indentwidth);
        void applyTranslation(std::vector<Base*> &bfuncs);
        virtual void apply(std::vector<Base*> &bfuncs);
        static bool parse_detail(iterator &it, const iterator &end, Function* &fun);
    };
}

#endif
