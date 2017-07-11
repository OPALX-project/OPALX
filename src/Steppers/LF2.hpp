template <typename FieldFunction, typename ... Arguments>
bool LF2<FieldFunction, Arguments ...>::advance(PartBunch* bunch,
                                                const size_t& i,
                                                const double& t,
                                                const double dt,
                                                Arguments& ... args) const
{
    bool flagNoDeletion = true;
    
    // push for first LF2 half step
    push_m(bunch->R[i], bunch->P[i], 0.5 * dt * 1.0e-9);  // ns --> s
    
    //BEGIN REMOVE
    bunch->setT(bunch->getT() + dt);

    // Path length update
    double dotP = dot(bunch->P[0], bunch->P[0]);
    double gamma = sqrt(1.0 + dotP);
    double PathLength_m = bunch->getLPath() + dt * 1.0e-9 * sqrt(dotP) * Physics::c / gamma;
    bunch->setLPath(PathLength_m);
    //END REMOVE
    
    flagNoDeletion = kick_m(bunch, i, t, dt * 1.0e-9, args ...);
    
    // push for second LF2 half step
    push_m(bunch->R[i], bunch->P[i], 0.5 * dt * 1.0e-9);  // ns --> s
    
    //BEGIN REMOVE
    bunch->setT(bunch->getT() + dt * 1.0e-9);

    // Path length update
    dotP = dot(bunch->P[0], bunch->P[0]);
    gamma = sqrt(1.0 + dotP);
    PathLength_m += dt * 1.0e-9 * sqrt(dotP) * Physics::c / gamma;
    bunch->setLPath(PathLength_m);
    //END REMOVE
    
    return flagNoDeletion;
}


template <typename FieldFunction, typename ... Arguments>
void LF2<FieldFunction, Arguments ...>::push_m(Vector_t& R, const Vector_t& P,
                                               const double& h) const
{
    double const gamma = sqrt(1.0 + dot(P, P));
    double const c_gamma = Physics::c / gamma;
    Vector_t const v = P * c_gamma;
    R += h * v;
}


template <typename FieldFunction, typename ... Arguments>
bool LF2<FieldFunction, Arguments ...>::kick_m(PartBunch* bunch, const size_t& i,
                                               const double& t, const double& h,
                                               Arguments& ... args) const
{
    Vector_t externalE = Vector_t(0.0, 0.0, 0.0);
    Vector_t externalB = Vector_t(0.0, 0.0, 0.0);
    
    bool outOfBound = this->fieldfunc_m(t, i, externalE, externalB, args ...);
    
    if ( outOfBound )
        return false;
    
    
    double const q = bunch->Q[0] / Physics::q_e; // For now all particles have the same charge
    double const M = bunch->M[0] * 1.0e9; // For now all particles have the same rest energy
    double const h12Halfqc_M = 0.5 * h * q * Physics::c / M;
    double const h12Halfqcc_M = h12Halfqc_M * Physics::c;
    
    // Half step E
    bunch->P[i] += h12Halfqc_M * externalE;

    // Full step B
    double const gamma = sqrt(1.0 + dot(bunch->P[i], bunch->P[i]));
    Vector_t const r = h12Halfqcc_M * externalB / gamma;
    Vector_t const w = bunch->P[i] + cross(bunch->P[i], r);
    Vector_t const s = 2.0 / (1.0 + dot(r, r)) * r;
    bunch->P[i] += cross(w, s);

    // Half step E
    bunch->P[i] += h12Halfqc_M * externalE;
    
    return true;
}
