#include "Distribution.h"
#include "Ippl.h"

/*
  Issues: - scan option i.e. when the particles are already allcated


*/

Distribution::Distribution(DistrTypeT actDist, size_t Np, double dEBins, int nBins, Vector_t distCutOff, 
			   bool hasLaserProfile,
			   string laserProfileFn,
			   string laserImage,
			   double intensityCut,
			   double tpulsefwhm, double trise, double tfall, Vector_t m,
			   double ekin,
			   Vektor<double,7> corr,
			   Vector_t r,
			   Vector_t p,
			   double avrgt,
			   double avrgpt):
    pbin_m(NULL),
    doEmission_m(false),
    nBins_m(nBins),
    distrType_m(actDist),
    laserProfileFn_m(laserProfileFn),
    laserImage_m(laserImage),
    intensityCut_m(intensityCut), 
    distCutOff_m(distCutOff),         //Attributes::getReal(itsAttr[TRANSVCUTOFF]);
    sigx_m(r),    
    sigp_m(p),
    corr_m(corr)
{ 
    if (hasLaserProfile)
	laserProfile_m = new LaserProfile(laserProfileFn_m, laserImage_m, intensityCut_m);
    else
	laserProfile_m = NULL;
	
    rGen_m = new RANLIB_class(265314159, 4);

    switch(distrType_m) {

    case GUNGAUSSFLATTOPTH:
    
	tPulseLengthFWHM_m = tpulsefwhm; // Attributes::getReal(itsAttr[TPULSEFWHM]);
	tRise_m = trise;                 // Attributes::getReal(itsAttr[TRISE]);
	tFall_m = tfall;                 // Attributes::getReal(itsAttr[TFALL]);

	double tratio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));
	sigmaRise_m = tRise_m / tratio;
	sigmaFall_m = tFall_m / tratio;
	
	tEmission_m = tPulseLengthFWHM_m + (distCutOff_m(2) - sqrt(2.0 * log(2.0))) * (sigmaRise_m + sigmaFall_m);
	tBin_m = tEmission_m / nBins_m;
	
	/*
	  prepare quantities for thermal emittance calculation
	*/
	workf_m    = 0.0;      // eV
	siglaser_m = 0.0;      // m
	elaser_m   = 0.0;      // eV
	fe_m       = 0.0;      // Fermi energy eV
	ag_m       = 0.0;      // Acceleration gradient eV/m

	schottky_m = 0.0;      // eV

	ekin_m = ekin;                            // eV Attributes::getReal(itsAttr[EK]);
	ptot_m = eVtoBetaGamma(ekin_m, 0.511);    // beta gamma
	break;
   
    case GAUSS: 
        avrgpt_m = avrgpt; // eVtoBetaGamma(Attributes::getReal(itsAttr[PT]), beam.getM());
	avrgt_m  = avrgt;   // Attributes::getReal(itsAttr[T]);

	gauss_offset_m[0] = 0.0;
	gauss_offset_m[1] = 0.0;
	break;
    
    case BINOMIAL: 
	binc_m = m;
	for(int j = 0; j < 3; j++) {
	    double chi = asin(corr_m[j]);
	    emit_m[j]  = sigx_m[j] * sigx_m[j] * cos(chi);
	    beta_m[j]  = sigx_m[j] * sigx_m[j] / emit_m[j] * distCutOff_m[j];  // distCutOff was 4.0 before
	    gamma_m[j] = sigp_m[j] * sigp_m[j] / emit_m[j] * distCutOff_m[j];  // distCutOff was 4.0 before;
	    alpha_m[j] = -corr_m[j] * sqrt(beta_m[j] * abs(gamma_m[j]));
	    
	    M_m[j]       =  sqrt(emit_m[j] * beta_m[j]);
	    PM_m[j]      =  sqrt(emit_m[j] * gamma_m[j]);
	    COSCHI_m[j]  =  sqrt(1.0 / (1.0 + alpha_m[j] * alpha_m[j]));
	    SINCHI_m[j]  = -alpha_m[j] * COSCHI_m[j];
	    CHI_m[j]     =  atan2(SINCHI_m[j], COSCHI_m[j]);
	    AMI_m[j]     =  1.0 / binc_m[j];
	    L_m[j]       =  sqrt((binc_m[j] + 1.0) / 2.0) * M_m[j];
	    PL_m[j]      =  sqrt((binc_m[j] + 1.0) / 2.0) * PM_m[j];
	}
	break;
    }
    if (nBins_m > 0) {
	doEmission_m = true;
	pbin_m = new PartBins((int) nBins_m);
	pbin_m->setRebinEnergy(1.0);
	h_m = gsl_histogram_alloc(nBins_m);
	createTimeBins(Np);
	pbin_m->setHistogram(h_m);
    }
    else {
	pbin_m = NULL;
	doEmission_m = false;
    }
}


Distribution::~Distribution() { 
    if(pbin_m) {
        delete pbin_m;
        pbin_m = NULL;
    }

    if(laserProfile_m) {
        delete laserProfile_m;
        laserProfile_m = NULL;
    }
}

pair<Vector_t, Vector_t> Distribution::sample(int binNumber) {
    Inform msg("sample  ");
    Vector_t r(0.0);
    Vector_t p(0.0);

    switch(distrType_m) {
	
    case GUNGAUSSFLATTOPTH: 

	double phi;
	double theta;
    
	double x, y, px, py, pt ;
	double xy = 6;    
	
	if(laserProfile_m != NULL) {
	    laserProfile_m->GetXY(&x, &y);
	    x = 2 * x - 1.0;
	    y = 2 * y - 1.0;
	} else {
	    while(xy > 1) {
		x  = rGen_m->uniform(-1.0, 1.0);
		y  = rGen_m->uniform(-1.0, 1.0);
		xy = sqrt(x * x + y * y);
	    }
	}
	
	x   *=  sigx_m(0);
	y   *=  sigx_m(1);

	/*
	  Now calculate the thermal emittances
	*/
	p = calcThermalEmittance (ptot_m);

	const double dtr = tBin_m*rGen_m->uniform(0, 1.0);
	const double t = binNumber*tBin_m;

	r = Vector_t(x, y, t + dtr);
	
	return pair<Vector_t,Vector_t>(r,p);
	break;

    case GAUSS:
	return sampleGauss();
	break;
    
    case BINOMIAL: 
	return sampleBinom();
	break;
    default:
	ERRORMSG("Distribution unknown" << endl;);
    }
}

void Distribution::create(myPartBunch &beam, size_t Np) {

    Inform msg("Distribution::create ");

    switch(distrType_m) {
    case BINOMIAL:    
        msg << "About to create Binomial distribution " << endl;
	createBinom(beam, Np);	
        msg << ".... done creating " << beam.getTotalNum() << " particles" << endl;
	break;

    case GAUSS:
        msg << "About to create Gaussian distribution (4 \sigma)" << endl;

        for(int i = 0; i < Np; i++) {
            double x, y;      // generate independent Gaussians, then correlate and finaly scale
            x  = rGen_m->gauss(0.0, 1.0);
            y  = rGen_m->gauss(0.0, 1.0);
            double xx = x;
            double yy = y;
            double px0  = x * corr_m[0] + y * sqrt(1.0 - corr_m[0] * corr_m[0]);
            double x0   =  x * sigx_m(0);  
            px0 *= sigp_m(0);

            x  = rGen_m->gauss(0.0, 1.0);
            y  = rGen_m->gauss(0.0, 1.0);
            double py0  = x * corr_m[1] + y * sqrt(1.0 - corr_m[1] * corr_m[1]);
            double y0   =  x * sigx_m(1);  
            py0 *= sigp_m(1);

            double del0;
            double psi0;

	    x  = rGen_m->gauss(0.0, 1.0);
	    y  = rGen_m->gauss(0.0, 1.0);
	    double l32 = (corr_m[6] - corr_m[0] * corr_m[5]) / sqrt(1.0 - corr_m[0] * corr_m[0]);
	    double l33 = sqrt(1 - corr_m[5] * corr_m[5] - l32 * l32);
	    psi0 = xx * corr_m[5] + yy * l32 + x * l33;
	    double l42 = (corr_m[4] - corr_m[0] * corr_m[3]) / sqrt(1.0 - corr_m[0] * corr_m[0]);
	    double l43 = (corr_m[2] - corr_m[5] * corr_m[3] - l42 * l32) / l33;
	    double l44 = sqrt(1 - corr_m[3] * corr_m[3] - l42 * l42 - l43 * l43);
	    del0 = xx * corr_m[3] + yy * l42 + x * l43 + y * l44;
	    psi0 = avrgt_m + psi0 * sigx_m(2);
	    del0 = avrgpt_m + sigp_m(2) * del0;

	    if (doEmission_m) {
		vector<double> tmp;
		tmp.push_back(x0);
		tmp.push_back(y0);
		tmp.push_back(0.0);
		tmp.push_back(px0);
		tmp.push_back(py0);
		tmp.push_back(del0);
		tmp.push_back(0);
		pbin_m->fill(tmp);
	    }
	    else { 
		/*
		  beam.create(1);
		  beam.R[count] = Vector_t(x0, y0, psi0);
		  beam.P[count] = Vector_t(px0, py0, del0);
		  beam.Bin[count] = 0; // not initialized
		  beam.PType[count] = 0;
		*/
		beam.R.push_back(Vector_t(x0, y0, psi0));
		beam.P.push_back(Vector_t(px0, py0, del0));
		beam.Q.push_back(beam.getChargePerParticle());
		beam.PType.push_back(0);
		beam.Bin.push_back(0);
	    }
        }
	break;

    case GUNGAUSSFLATTOPTH:    
	pbin_m->sortArrayT();
	beam.setPBins(pbin_m);
	break;
    }
}

void Distribution::sampleGauss(myPartBunch &beam, size_t Np) {

    for(size_t i = beam.getTotalNum(); i < Np; i++) {

	pair<Vector_t, Vector_t> xp = sampleGauss();	
	const Vector_t x = xp.first;
	const Vector_t p = xp.first;
	
	/*
	  beam.create(1);
	  beam.R[count] = Vector_t(x0, y0, psi0);
	  beam.P[count] = Vector_t(px0, py0, del0);
	  beam.Bin[count] = 0; // not initialized
	  beam.Q[count] = beam.getChargePerParticle();
	  beam.PType[count] = 0;
	*/
	beam.R.push_back(x);
	beam.P.push_back(p);
	beam.Q.push_back(beam.getChargePerParticle());
	beam.PType.push_back(0);
	beam.Bin.push_back(0);
    }
}


pair<Vector_t, Vector_t> Distribution::sampleGauss() {

    double x, y;      // generate independent Gaussians, then correlate and finaly scale
    x  = rGen_m->gauss(0.0, 1.0);
    y  = rGen_m->gauss(0.0, 1.0);
    double xx = x;
    double yy = y;
    double px0  = x * corr_m[0] + y * sqrt(1.0 - corr_m[0] * corr_m[0]);
    double x0   = x * sigx_m(0) + gauss_offset_m[0];
    px0 *= sigp_m(0);
    
    x  = rGen_m->gauss(0.0, 1.0);
    y  = rGen_m->gauss(0.0, 1.0);
    double py0  = x * corr_m[1] + y * sqrt(1.0 - corr_m[1] * corr_m[1]);
    double y0   = x * sigx_m(1) + gauss_offset_m[1];
    py0 *= sigp_m(1);
    
    double del0;
    double psi0;
    
    x  = rGen_m->gauss(0.0, 1.0);
    y  = rGen_m->gauss(0.0, 1.0);
    const double l32 = (corr_m[6] - corr_m[0] * corr_m[5]) / sqrt(1.0 - corr_m[0] * corr_m[0]);
    const double l33 = sqrt(1 - corr_m[5] * corr_m[5] - l32 * l32);
    psi0 = xx * corr_m[5] + yy * l32 + x * l33;
    const double l42 = (corr_m[4] - corr_m[0] * corr_m[3]) / sqrt(1.0 - corr_m[0] * corr_m[0]);
    const double l43 = (corr_m[2] - corr_m[5] * corr_m[3] - l42 * l32) / l33;
    const double l44 = sqrt(1 - corr_m[3] * corr_m[3] - l42 * l42 - l43 * l43);
    del0 = xx * corr_m[3] + yy * l42 + x * l43 + y * l44;
    psi0 = avrgt_m + psi0 * sigx_m(2);
    del0 = avrgpt_m + sigp_m(2) * del0;
    
    return pair<Vector_t,Vector_t>(Vector_t(x0, y0, psi0),Vector_t(px0, py0, del0));
}



void Distribution::createTimeBins(const int Np) {

    gsl_histogram_set_ranges_uniform(h_m, 0, tEmission_m);
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    const double sq2pi = sqrt(2.0 * Physics::pi);
    const double tFlat = tPulseLengthFWHM_m - sqrt(2.0 * log(2.0)) * (sigmaRise_m + sigmaFall_m);
    const double totA = tFlat + 0.5 * sq2pi * (sigmaRise_m + sigmaFall_m);
    const int nrise = Np * 0.5 * gsl_sf_erf(distCutOff_m(2) / sqrt(2.0)) * sq2pi * sigmaRise_m / totA;
    const int nfall = Np * 0.5 * gsl_sf_erf(distCutOff_m(2) / sqrt(2.0)) * sq2pi * sigmaFall_m / totA;
    const int nflat = Np - nrise - nfall;
    
    // Rise: [0, c\sigma_R]
    for(int i = 0; i < nrise; i++) {
	double r1 = gsl_ran_gaussian_tail(r, 0, sigmaRise_m);
	while(r1 > distCutOff_m(2) * sigmaRise_m)
	    r1 = gsl_ran_gaussian_tail(r, 0, sigmaRise_m);
	gsl_histogram_increment(h_m, -r1 + distCutOff_m(2) * sigmaRise_m);
    }
    // Fall: [c\sigma_R + tFlat, c\sigma_R + tFlat + c\sigma_F]
    for(int i = 0; i < nfall; i++) {
	double r1 = gsl_ran_gaussian_tail(r, 0, sigmaFall_m);
	while(r1 > distCutOff_m(2) * sigmaFall_m)
	    r1 = gsl_ran_gaussian_tail(r, 0, sigmaFall_m);
	gsl_histogram_increment(h_m, r1 + distCutOff_m(2) * sigmaRise_m + tFlat);
    }
    // Flattop: [c\sigma_R, c\sigma_R + tFlat]
    for(int i = 0; i < nflat; i++) {
	const double r1 = gsl_ran_flat(r, 0, tFlat);
	gsl_histogram_increment(h_m, r1 + distCutOff_m(2) * sigmaRise_m);
    }
    gsl_rng_free(r);

    FILE *fp;
    fp = fopen("thisto.dat","w");
    gsl_histogram_fprintf(fp,h_m,"%g","%g");
    fclose(fp);
}

pair<Vector_t, Vector_t> Distribution::sampleBinom() {
    const double two_pi = 2.0 * 4.0 * atan(1.0);
    Inform msg("Distribution ");
    
    Vector_t x;
    Vector_t p;

    double betagamma_part, pos_part;
    
    if(doEmission_m) {
        const double mass = 1.0; // FixMe NO MASS HERE
        const double gamma_part = 1. + ekin_m / mass;
        
	betagamma_part = sqrt(ekin_m * ekin_m / (mass * mass) + 2.* ekin_m / mass);
	pbin_m->setGamma(gamma_part);
    } else {
        betagamma_part = sigp_m[2];
        pos_part = sigx_m[2];
    }

    double S1, S2;
    double A, AL, U, Ucp, V, Vcp;

    S1 = IpplRandom();
    S2 = IpplRandom();
    if (binc_m[0] <= 10000) {
	A = sqrt(1.0 - pow(S1, AMI_m[0]));
	AL = two_pi * S2;
	U = A * cos(AL);
	V = A * sin(AL);
	Ucp = U;
	Vcp = V;
	x[0] = L_m[0] * U;
	p[0] = PL_m[0] * (U * SINCHI_m[0] + V * COSCHI_m[0]);
    }
    else {
	A = sqrt(2.0) / 2.0 * sqrt(-log(S1));
	AL = two_pi * S2;
	U = A * cos(AL);
	V = A * sin(AL);
	Ucp = U;
	Vcp = V;
	x[0] = M_m[0] * U;
	p[0] = PM_m[0] * (U * SINCHI_m[0] + V * COSCHI_m[0]);
    }

    S1 = IpplRandom();
    S2 = IpplRandom();
    if (binc_m[1] <= 10000) {
	A = sqrt(1.0 - pow(S1, AMI_m[1]));
	AL = two_pi * S2;
	U = A * cos(AL);
	V = A * sin(AL);
	x[1] = L_m[1] * U;
	p[1] = PL_m[1] * (U * SINCHI_m[1] + V * COSCHI_m[1]);
    }
    else {
	A = sqrt(2.0) / 2.0 * sqrt(-log(S1));
	AL = two_pi * S2;
	U = A * cos(AL);
	V = A * sin(AL);
	x[1] = M_m[1] * U;
	p[1] = PM_m[1] * (U * SINCHI_m[1] + V * COSCHI_m[1]);
    }

    S1 = IpplRandom();
    S2 = IpplRandom();
    if (binc_m[2] <= 10000) {
	A = sqrt(1.0 - pow(S1, AMI_m[2]));
	AL = two_pi * S2;
	U = A * cos(AL);
	V = A * sin(AL);
	const double l32 = (corr_m[6] - corr_m[0] * corr_m[5]) / sqrt(1.0 - corr_m[0] * corr_m[0]);
	const double l33 = sqrt(1 - corr_m[5] * corr_m[5] - l32 * l32);
	x[2] = Ucp * corr_m[5] + Vcp * l32 + U * l33;
	
	const double l42 = (corr_m[4] - corr_m[0] * corr_m[3]) / sqrt(1.0 - corr_m[0] * corr_m[0]);
	const double l43 = (corr_m[2] - corr_m[5] * corr_m[3] - l42 * l32) / l33;
	const double l44 = sqrt(1 - corr_m[3] * corr_m[3] - l42 * l42 - l43 * l43);
	
	p[2] = Ucp * corr_m[3] + Vcp * l42 + U * l43 + V * l44;
	x[2] *= L_m[2];
	p[2] *= PL_m[2];
	
    } else {
	A = sqrt(2.0) / 2.0 * sqrt(-log(S1));
	AL = two_pi * S2;
	U = A * cos(AL);
	V = A * sin(AL);
	const double l32 = (corr_m[6] - corr_m[0] * corr_m[5]) / sqrt(1.0 - corr_m[0] * corr_m[0]);
	const double l33 = sqrt(1 - corr_m[5] * corr_m[5] - l32 * l32);
	x[2] = Ucp * corr_m[5] + Vcp * l32 + U * l33;
	const double l42 = (corr_m[4] - corr_m[0] * corr_m[3]) / sqrt(1.0 - corr_m[0] * corr_m[0]);
	const double l43 = (corr_m[2] - corr_m[5] * corr_m[3] - l42 * l32) / l33;
	const double l44 = sqrt(1 - corr_m[3] * corr_m[3] - l42 * l42 - l43 * l43);
	p[2] = Ucp * corr_m[3] + Vcp * l42 + U * l43 + V * l44;
	x[2]  *= M_m[2];
	p[2] *= PM_m[2];
	
    }
    // p[2] += betagamma_part;
    // x[2] += pos_part;
    
    return pair<Vector_t,Vector_t>(x,p);
}

void Distribution::createBinom(myPartBunch &beam, size_t particles) {

    const double two_pi = 2.0 * 4.0 * atan(1.0);

    Vector_t x;
    Vector_t p;

    double betagamma_part, pos_part;
    
    if(doEmission_m) {
        const double mass = beam.getM();
        const double gamma_part = 1. + ekin_m / mass;
        
	betagamma_part = sqrt(ekin_m * ekin_m / (mass * mass) + 2.* ekin_m / mass);
	pbin_m->setGamma(gamma_part);
    } else {
        betagamma_part = sigp_m[2];
        pos_part = sigx_m[2];
    }

    for(size_t n = 0; n < particles; ++n) {

	pair<Vector_t, Vector_t> xp = sampleBinom();	

	const Vector_t x = xp.first;
	const Vector_t p = xp.second;

	if(doEmission_m) {
	    vector<double> tmp;
	    tmp.push_back(x[0]);
	    tmp.push_back(x[1]);
	    tmp.push_back(x[2]);
	    tmp.push_back(p[0]);
	    tmp.push_back(p[1]);
	    tmp.push_back(p[2]);
	    tmp.push_back(0);
	    pbin_m->fill(tmp);
	} else {
	    /*
	      beam.create(1);
	      beam.R[count] = x;
	      beam.P[count] = p;
	      beam.Q[count] = beam.getChargePerParticle();
	      beam.PType[count] = 0;
	    */
	    beam.R.push_back(x);
	    beam.P.push_back(p);
	    beam.Q.push_back(beam.getChargePerParticle());
	    beam.PType.push_back(0);
	    beam.Bin.push_back(0);
	}
    }
}


Vector_t Distribution::calcThermalEmittance (double pTherm) {

    const double bega  = 0.0;
    const double phi   = 2.0 * acos(sqrt(rGen_m->uniform(0.0, 1.0)));
    const double theta = 2.0 * Physics::pi * rGen_m->uniform(0.0, 1.0);

    const Vector_t p (pTherm * sin(phi) * cos(theta),
		      pTherm * sin(phi) * sin(theta),
	      bega + (pTherm * abs(cos(phi))));
    return p;
}


Inform &Distribution::print(Inform &os) const {

    switch (distrType_m) {

    case GAUSS:
	os << "* DISTRIBUTION is GAUSS (3D) ---------------------------------------------------------" << endl;
	os << "* distribution cut-off  " << distCutOff_m << " [sigma] " << endl;
	
	os << "* sigmax=\t" << sigx_m(0) << " [m];\t"
	   << "sigmay=\t"   << sigx_m(1) << " [m];\t"
	   << "sigmat=\t"   << sigx_m(2) << " [m];" << endl;
 
	os << "* sigmapx=\t"     << sigp_m(0) << " [eV];\t"
	   << "sigmapy=\t"       << sigp_m(1) << " [eV];\t"
	   << "pt +- sigmapt=\t" << pt_m << "+-" <<  sigp_m(2) << " [eV]" << endl;
 
	os << "* corr x-px=\t" << corr_m[0];
	os << "\t corr y-py=\t" << corr_m[1];
	os << "\t corr t-pt=\t" << corr_m[2] << endl;

	os << "*  R61=\t" << corr_m[3];
	os << "\t R62=\t" << corr_m[4];
	os << "\t R51=\t" << corr_m[5];
	os << "\t R52=\t" << corr_m[6] << endl;

	break;

    case BINOMIAL:
	os << "* DISTRIBUTION is BINOMIAL ---------------------------------------------------------" << endl;
	os << "* Mx=My=Mt " << binc_m << endl;
	os << "* distribution cut-off  " << distCutOff_m << " [sigma] " << endl;
	
	os << "* sigmax=\t" << sigx_m(0) << " [m];\t"
	   << "sigmay=\t"   << sigx_m(1) << " [m];\t"
	   << "sigmat=\t"   << sigx_m(2) << " [m];" << endl;
 
	os << "* sigmapx=\t"     << sigp_m(0) << " [eV];\t"
	   << "sigmapy=\t"       << sigp_m(1) << " [eV];\t"
	   << "pt +- sigmapt=\t" << pt_m << "+-" <<  sigp_m(2) << " [eV]" << endl;
 
	os << "* corr x-px=\t" << corr_m[0];
	os << "\t corr y-py=\t" << corr_m[1];
	os << "\t corr t-pt=\t" << corr_m[2] << endl;

	os << "*  R61=\t" << corr_m[3];
	os << "\t R62=\t" << corr_m[4];
	os << "\t R51=\t" << corr_m[5];
	os << "\t R52=\t" << corr_m[6] << endl;

	break;

    case GUNGAUSSFLATTOPTH:    
	os << "* DISTRIBUTION is GUNGAUSSFLATTOPTH ------------------------------------------------" << endl;
	os << "* GUNGAUSS FLAT TOP &  THERMAL EMITTANCE a la ASTRA " << endl;
	os << "* Kinetic energy (thermal emittance) = " << ekin_m << " [eV]  " << endl;
	os << "* tBin = " << tBin_m << " [sec]  nBins = " << nBins_m << " tEmission =  " << tEmission_m << " [sec] " << endl;
	os << "* distribution cut-off  " << distCutOff_m << " [sigma] " << endl;

	os << "* sigmax=\t" << sigx_m(0) << " [m];\t"
	   << "sigmay=\t"   << sigx_m(1) << " [m];\t"
	   << "sigmat=\t"   << sigx_m(2) << " [m];" << endl;
	
	os << "* sigmapx=\t"     << sigp_m(0) << " [eV];\t"
	   << "sigmapy=\t"       << sigp_m(1) << " [eV];\t"
	   << "pt +- sigmapt=\t" << pt_m << "+-" <<  sigp_m(2) << " [eV]" << endl;
	
	os << "* corr x-px=\t" << corr_m[0];
	os << "\t corr y-py=\t" << corr_m[1];
	os << "\t corr t-pt=\t" << corr_m[2] << endl;
	
	os << "*  R61=\t" << corr_m[3];
	os << "\t R62=\t" << corr_m[4];
	os << "\t R51=\t" << corr_m[5];
	os << "\t R52=\t" << corr_m[6] << endl;

	if(laserProfile_m) {
	    os << "* Laser profile: " << laserProfileFn_m 
	       << " Image: " << laserImage_m 
	       << " Intensity cut: " << intensityCut_m << endl;
	}
	break;
     
    default:
	os << "* DISTRIBUTION is UNKNOWN ----------------------------------------------------------" << endl;
    }
    os << "* --------------------------------------------------------------------------------- " << endl;
}

//inline Inform &operator<<(Inform &os, const Distribution &d) {
//    return d.print(os);
//}





