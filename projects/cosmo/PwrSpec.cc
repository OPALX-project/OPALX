

template <class T, unsigned int Dim>
PwrSpec<T,Dim>::PwrSpec(ChargedParticles<T,Dim> *beam):
  layout_m(&(beam->getFieldLayout())),
  mesh_m(&(beam->getMesh()))
{
  gDomain_m = layout_m->getDomain();        // global domain
  lDomain_m = layout_m->getLocalNDIndex();  // local domain

  fC_m.initialize(*mesh_m, *layout_m);
  pwrSpec_m.initialize(*mesh_m, *layout_m);

  // create the FFT object
  bool compressTemps = true;
  fft_m = new FFT_t(layout_m->getDomain(), compressTemps);
  fft_m->setDirectionName(+1, "forward");
  fft_m->setDirectionName(-1, "inverse");
  
}


template <class T, unsigned int Dim>
void PwrSpec<T,Dim>::cacl1D(ChargedParticles<T,Dim> *b) 
{
    
  /// b->rho_m == charge density scattered onto grid with CIC
  fC_m = b->rho_m;
  
  fft_m->transform( "forward", fC_m);
  
  pwrSpec_m = fC_m*conj(fC_m);


  /// now calculate the 1D pwr-spectrum
  NDIndex<3> elem;
  NDIndex<1> el1d;

  
  for (int i=lDomain_m[0].min(); i<=lDomain_m[0].max(); ++i) {
    elem[0]=Index(i,i);
    el1d[0]=Index(i);
    for (int j=lDomain_m[1].min(); j<=lDomain_m[1].max(); ++j) {
      elem[1]=Index(j,j);
      for (int k=lDomain_m[2].min(); k<=lDomain_m[2].max(); ++k) {
        elem[2]=Index(k,k);
	//	f1D[el1d] += pwrSpec_m.localElement(elem);
      }
    }
  }
}


template <class T, unsigned int Dim>
void PwrSpec<T,Dim>::save1DSpectra(string fn) 
{
  Inform o;
  NDIndex<1> el;
  
  for (int i=lDomain_m[0].min(); i<=lDomain_m[0].max(); ++i) {
    el[0]=Index(i);
    //    o << i << " \t" << f1D_m.localElement(el) << endl;
  }
}
