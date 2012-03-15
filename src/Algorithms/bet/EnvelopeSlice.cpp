/* slice.C
   slice class definition

   Project: Beam Envelope Tracker (BET)

   Revision history
   Date          Description                                     Programmer
   ------------  --------------------------------------------    --------------
   07-03-06      Created                                         Rene Bakker

   Last Revision:
   $Id: slice.C 175 2007-12-12 17:48:32Z bakker $
*/

#define SVN_DATE "$Date: 2007-12-12 18:48:32 +0100 (Wed, 12 Dec 2007) $"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <Physics/Physics.h>

#include "Algorithms/bet/error.h"
#include "Algorithms/bet/EnvelopeSlice.h"
/**
#include "utils.h"
*/

static double *vector(int nn)
{
    double *r;
    r = (double *) malloc(sizeof(double)*nn);
    if (!r) 
        writeError();
    return r;
}

EnvelopeSlice::EnvelopeSlice(FILE *f) 
{
    double g;

    /** 
      eList     = NULL;
      */
    p         = vector(SLNPAR);
    p_old     = vector(SLNPAR);
    p_scr     = NULL;
    p_local   = 1;

    if (p == NULL) {
        writeError();
    }

    if (f) {
        int nScan = read(f);
        if (nScan < (SLNPAR+1)) {
            writeError();
        }
    } else {
        g  = 1.0 + (1.0e6*Physics::q_e/(Physics::EMASS*Physics::c*Physics::c));
        p[SLI_beta] = sqrt(1.0/(1.0-(1.0/g*g)));
        p[SLI_z]    = 0.0;
        p[SLI_x]    = 1.0e-3;
        p[SLI_y]    = 1.0e-3;
        p[SLI_px]   = 0.0;
        p[SLI_py]   = 0.0;
        p[SLI_x0]   = 0.0;
        p[SLI_y0]   = 0.0;
        p[SLI_px0]  = 0.0;
        p[SLI_py0]  = 0.0;
        valid       = 1;
    }
    backup();
} /* Slice */

EnvelopeSlice::~EnvelopeSlice() 
{
    free(p_old);
    if (p_scr)   free(p_scr);
    if (p_local) free(p);
}


void EnvelopeSlice::setBuffer(double *b)
{
    int i;

    if(b) {
        memcpy(b,p,sizeof(double)*SLNPAR);
        if (p_local) free(p);
        p = b;
        p_local = 0;
    }else {
        writeError();
    }
}

void EnvelopeSlice::writeHeader(FILE *f) 
{
    fprintf(f,"%-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s\t %-20s\n",
            "z [m]","beta [-]","x [m]","px/beta/c [-]","y [m]","py/beta/c [-]",
            "x0 [m]","px0/beta/c [m/s]","y0 [m]","py0/beta/c [m/s]",
            "is_valid");
}


void EnvelopeSlice::write(FILE *f) 
{
    double bc = p[SLI_beta]*Physics::c;

    fprintf(f,
            "%20.14le \t %20.14le \t %20.14le \t %20.14le \t %20.14le \t %20.14le \t %20.14le \t %20.14le \t %20.14le \t %20.14le \t %2d\n",
            p[SLI_z],p[SLI_beta],
            p[SLI_x],p[SLI_px]/bc,
            p[SLI_y],p[SLI_py]/bc,
            p[SLI_x0],p[SLI_px0]/bc,
            p[SLI_y0],p[SLI_py0]/bc,
            valid);
}


int EnvelopeSlice::read(FILE *f) 
{
    char line[2048];
    int nScan;
    double bc;
    fgets(line,2048,f);
    nScan = sscanf(line,
            "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
            &p[SLI_z],&p[SLI_beta],
            &p[SLI_x],&p[SLI_px],
            &p[SLI_y],&p[SLI_py],
            &p[SLI_x0],&p[SLI_px0],
            &p[SLI_y0],&p[SLI_py0],
            &valid);

    bc = p[SLI_beta]*Physics::c;
    if(bc == 0.0) {
        if(nScan > 2) 
            nScan = 2;
        bc = Physics::c;
    }

    p[SLI_px]  *= bc;
    p[SLI_py]  *= bc;
    p[SLI_px0] *= bc;
    p[SLI_py0] *= bc;

    return nScan;
}

double EnvelopeSlice::gamma() 
{
    return sqrt(1.0/(1.0-pow(p[SLI_beta],2)));
}


/** 

  void Slice::findActive(double dt,Element *fe) {
  Element 
 *e = fe;
 double
 z    = p[SLI_z],
 zMax = z + c*dt;

 if (eList) {
 delete eList;
 eList = NULL;
 }

 while (e) {
 if (e->isActive(z,zMax)) {
 eList = new ElementList(e,eList);
 }
 e = e->getNext();
 }
 } 

*/

void EnvelopeSlice::backup() 
{
    was_valid = valid;
    memcpy(p_old,p,sizeof(double)*SLNPAR);
}

void EnvelopeSlice::restore() 
{
    valid = was_valid;
    memcpy(p,p_old,sizeof(double)*SLNPAR);
}

int EnvelopeSlice::mark(double z,double t) 
{
    if((p_old[SLI_z] <= z) && (p[SLI_z] >= z)) {
        int i;
        double z0 = p_old[SLI_z], dz = p[SLI_z] - z0;

        if(p_scr == NULL) 
            p_scr = vector(SLNPAR);
        if(dz == 0.0) {
            for(i=0; i<SLNPAR; i++) 
                p_scr[i] = 0.5*(p[i] + p_old[i]);
        } else {
            for(i=0; i<SLNPAR; i++) 
                p_scr[i] = p_old[i] + (p[i] - p_old[i])*(z - z0)/dz;
        }
        // replace z with time
        p_scr[SLI_z] = t - 2.0*(p[SLI_z]-z)/(Physics::c*(p[SLI_beta]+p_scr[SLI_beta]));
        return 1;
    } else if((p_old[SLI_z] > z) && (p_scr == NULL)) {
        return -1;
    } else {
        return (p_scr == NULL?0:1);
    }
}

void EnvelopeSlice::clear_mark() 
{
    if(p_scr) {
        free(p_scr);
        p_scr = NULL;
    }
}

int EnvelopeSlice::check() 
{
    int changed = 0;

    if(valid) {
        valid = (p[SLI_beta] > 0.0);
        if(!valid)
            changed = 1;
    }
    return changed;
}

int EnvelopeSlice::is_valid() 
{
    return valid;
}


