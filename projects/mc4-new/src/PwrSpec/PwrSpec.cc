#include "PwrSpec.hh"

template <class T, unsigned int Dim>
PwrSpec<T,Dim>::PwrSpec(ChargedParticles<T,Dim> *beam, SimData<T,Dim> simData):
    layout_m(&(beam->getFieldLayout())),
    simData_m(simData)
{
    doInit();
}

template <class T, unsigned int Dim>
void PwrSpec<T,Dim>::doInit() {

    gDomain_m = layout_m->getDomain();        // global domain
    lDomain_m = layout_m->getLocalNDIndex();  // local domain

    kmax_m = nint(sqrt(3*(simData_m.ng_comp*simData_m.ng_comp/4))) + 1 ;
    ng_m   = simData_m.ng_comp;

    spectra1D_m = (T *) malloc (kmax_m*sizeof(T));
    Nk_m = (int *) malloc (kmax_m*sizeof(int));

    for (int k=0; k<kmax_m; ++k) {
        spectra1D_m[k] = 0.0;
        Nk_m[k]=0;
    }
}

template <class T, unsigned int Dim>
void PwrSpec<T,Dim>::calc1D(ChargedParticles<T,Dim> *univ)
{
    unsigned int kk;

    for (int i=0;i<kmax_m;i++) {
        Nk_m[i]=0;
        spectra1D_m[i] = 0.0;
    }

    INFOMSG("Sum psp=real( ... " << sum(univ->rho_m) << " kmax= " << kmax_m << endl;);


    /*
       Here scatter_+^%^$#% takes place

*/
    NDIndex<Dim> loop;
    // Octant 1
    for (int i=gDomain_m[0].min(); i<(gDomain_m[0].max()+1)/2; ++i) {
        loop[0]=Index(i,i);
        for (int j=gDomain_m[1].min(); j<(gDomain_m[1].max()+1)/2; ++j) {
            loop[1]=Index(j,j);
            for (int k=gDomain_m[2].min(); k<(gDomain_m[2].max()+1)/2; ++k) {
                loop[2]=Index(k,k);
                if (lDomain_m.contains(loop)) {
                    kk=(int)nint(sqrt(i*i+j*j+k*k));
                    if (kk>kmax_m) {
                        INFOMSG("kk>kmax_m kk= " << kk << endl);
                        kk = kmax_m;
                    }
                    spectra1D_m[kk] += univ->rho_m.localElement(loop);
                    Nk_m[kk]++;
                }
            }
        }
    }
    // Octant 2
    for (int i=(gDomain_m[0].max()+1)/2; i<=gDomain_m[0].max(); ++i) {
        loop[0]=Index(i,i);
        for (int j=gDomain_m[1].min(); j<(gDomain_m[1].max()+1)/2; ++j) {
            loop[1]=Index(j,j);
            for (int k=gDomain_m[2].min(); k<(gDomain_m[2].max()+1)/2; ++k) {
                loop[2]=Index(k,k);
                if (lDomain_m.contains(loop)) {
                    kk=(int)nint(sqrt((i-ng_m)*(i-ng_m)+(j*j)+(k*k)));
                    if (kk>kmax_m) {
                        INFOMSG("kk>kmax_m kk= " << kk << endl);
                        kk = kmax_m;
                    }
                    spectra1D_m[kk] += univ->rho_m.localElement(loop);
                    Nk_m[kk]++;
                }
            }
        }
    }
    // Octant 3
    for (int i=gDomain_m[0].min(); i<(gDomain_m[0].max()+1)/2; ++i) {
        loop[0]=Index(i,i);
        for (int j=(gDomain_m[1].max()+1)/2; j<=gDomain_m[1].max(); ++j) {
            loop[1]=Index(j,j);
            for (int k=gDomain_m[2].min(); k<(gDomain_m[2].max()+1)/2; ++k) {
                loop[2]=Index(k,k);
                if (lDomain_m.contains(loop)) {
                    kk=(int)nint(sqrt((i*i)+(j-ng_m)*(j-ng_m)+(k*k)));
                    if (kk>kmax_m) {
                        INFOMSG("kk>kmax_m kk= " << kk << endl);
                        kk = kmax_m;
                    }
                    spectra1D_m[kk] += univ->rho_m.localElement(loop);
                    Nk_m[kk]++;
                }
            }
        }
    }
    //Octant 4
    for (int i=gDomain_m[0].min(); i<(gDomain_m[0].max()+1)/2; ++i) {
        loop[0]=Index(i,i);
        for (int j=gDomain_m[1].min(); j<(gDomain_m[1].max()+1)/2; ++j) {
            loop[1]=Index(j,j);
            for (int k=(gDomain_m[2].max()+1)/2; k<=gDomain_m[2].max(); ++k) {
                loop[2]=Index(k,k);
                if (lDomain_m.contains(loop)) {
                    kk=(int)nint(sqrt((i*i)+(j*j)+(k-ng_m)*(k-ng_m)));
                    if (kk>kmax_m) {
                        INFOMSG("kk>kmax_m kk= " << kk << endl);
                        kk = kmax_m;
                    }
                    spectra1D_m[kk] += univ->rho_m.localElement(loop);
                    Nk_m[kk]++;
                }
            }
        }
    }
    //Octant 5
    for (int i=gDomain_m[0].min(); i<(gDomain_m[0].max()+1)/2; ++i) {
        loop[0]=Index(i,i);
        for (int j=(gDomain_m[1].max()+1)/2; j<=gDomain_m[1].max(); ++j) {
            loop[1]=Index(j,j);
            for (int k=(gDomain_m[2].max()+1)/2; k<=gDomain_m[2].max(); ++k) {
                loop[2]=Index(k,k);
                if (lDomain_m.contains(loop)) {
                    kk=(int)nint(sqrt((i*i)+(j-ng_m)*(j-ng_m)+(k-ng_m)*(k-ng_m)));
                    if (kk>kmax_m) {
                        INFOMSG("kk>kmax_m kk= " << kk << endl);
                        kk = kmax_m;
                    }
                    spectra1D_m[kk]+=univ->rho_m.localElement(loop);
                    Nk_m[kk]++;
                }
            }
        }
    }
    //Octant 6
    for (int i=(gDomain_m[0].max()+1)/2; i<=gDomain_m[0].max(); ++i) {
        loop[0]=Index(i,i);
        for (int j=(gDomain_m[1].max()+1)/2; j<=gDomain_m[1].max(); ++j) {
            loop[1]=Index(j,j);
            for (int k=gDomain_m[2].min(); k<(gDomain_m[2].max()+1)/2; ++k) {
                loop[2]=Index(k,k);
                if (lDomain_m.contains(loop)) {
                    kk=(int)nint(sqrt((i-ng_m)*(i-ng_m)+(j-ng_m)*(j-ng_m)+(k*k)));
                    if (kk>kmax_m) {
                        INFOMSG("kk>kmax_m kk= " << kk << endl);
                        kk = kmax_m;
                    }
                    spectra1D_m[kk]+=univ->rho_m.localElement(loop);
                    Nk_m[kk]++;
                }
            }
        }
    }
    //Octant 7
    for (int i=(gDomain_m[0].max()+1)/2; i<=gDomain_m[0].max(); ++i) {
        loop[0]=Index(i,i);
        for (int j=gDomain_m[1].min(); j<(gDomain_m[1].max()+1)/2; ++j) {
            loop[1]=Index(j,j);
            for (int k=(gDomain_m[2].max()+1)/2; k<=gDomain_m[2].max(); ++k) {
                loop[2]=Index(k,k);
                if (lDomain_m.contains(loop)) {
                    kk=(int)nint(sqrt((i-ng_m)*(i-ng_m)+(j*j)+(k-ng_m)*(k-ng_m)));
                    if (kk>kmax_m) {
                        INFOMSG("kk>kmax_m kk= " << kk << endl);
                        kk = kmax_m;
                    }
                    spectra1D_m[kk]+=univ->rho_m.localElement(loop);
                    Nk_m[kk]++;
                }
            }
        }
    }
    //Octant 8
    for (int i=(gDomain_m[0].max()+1)/2; i<=gDomain_m[0].max(); ++i) {
        loop[0]=Index(i,i);
        for (int j=(gDomain_m[1].max()+1)/2; j<=gDomain_m[1].max(); ++j) {
            loop[1]=Index(j,j);
            for (int k=(gDomain_m[2].max()+1)/2; k<=gDomain_m[2].max(); ++k) {
                loop[2]=Index(k,k);
                if (lDomain_m.contains(loop)) {
                    kk=(int)nint(sqrt((i-ng_m)*(i-ng_m)+(j-ng_m)*(j-ng_m)+(k-ng_m)*(k-ng_m)));
                    if (kk>kmax_m) {
                        INFOMSG("kk>kmax_m kk= " << kk << endl);
                        kk = kmax_m;
                    }
                    spectra1D_m[kk]+=univ->rho_m.localElement(loop);
                    Nk_m[kk]++;
                }
            }
        }
    }
    INFOMSG("Loops done" << endl);
    reduce( &(Nk_m[0]), &(Nk_m[0]) + kmax_m , &(Nk_m[0]) ,OpAddAssign());
    reduce( &(spectra1D_m[0]), &(spectra1D_m[0]) + kmax_m, &(spectra1D_m[0]) ,OpAddAssign());
}


template <class T, unsigned int Dim>
void PwrSpec<T,Dim>::save1DSpectra(string fn)
{

    Inform* fdip = new Inform(NULL,fn.c_str(),Inform::OVERWRITE,0);
    Inform& fdi = *fdip;
    setInform(fdi);
    setFormat(9,1,0);

    T tpiL = 8.0*atan(1.0)/simData_m.rL;   // 2 pi / rL

    // Renormalize power spectrum to match mc2.
    int sumNk = 0;
    for (int i=0; i<kmax_m;i++) {
        sumNk+=Nk_m[i];
        spectra1D_m[i] /= 1.0*Nk_m[i];
        spectra1D_m[i] = spectra1D_m[i]*std::pow((T)(1.0*simData_m.np),(T)3.0);
    }

    fdi << "# " ;
    for (int i=0; i < kmax_m; i++)
        fdi << "Nk[" << i << "]= " << Nk_m[i] << " ";
    fdi << " sum Nk= " << sumNk << endl;

    for (int i=0; i<simData_m.np/2;i++) {
        spectra1D_m[i] = spectra1D_m[i+1];
        fdi  << (i+1)*tpiL << "   " << spectra1D_m[i]*std::pow((T)(simData_m.rL/simData_m.ng_comp),(T)3.0) << "  " << endl;
    }
    delete fdip;
}

