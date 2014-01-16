        // --- computeExternalFields();
        IpplTimings::startTimer(timeFieldEvaluation_m);
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {

            long ls = itsBunch->LastSection[i];
            itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], ls);
            if(ls != itsBunch->LastSection[i]) {
                if(!itsOpalBeamline_m.section_is_glued_to(itsBunch->LastSection[i], ls)) {
                    itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(ls), itsOpalBeamline_m.getSectionStart(ls));
                }
                itsBunch->LastSection[i] = ls;
            }

            Vector_t externalE, externalB;
            const unsigned long rtv = itsOpalBeamline_m.getFieldAt(i, itsBunch->R[i], ls, itsBunch->getT() + itsBunch->dt[i] / 2., externalE, externalB);
            if (rtv)
                ;
            itsBunch->Ef[i] += externalE;
            itsBunch->Bf[i] += externalB;
        }
        IpplTimings::stopTimer(timeFieldEvaluation_m);
        // --- computeExternalFields

