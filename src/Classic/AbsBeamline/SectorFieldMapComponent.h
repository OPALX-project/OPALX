
class SectorFieldMapComponent: public SBend {

public:
    /// Constructor with given name.
    SectorFieldMapComponent(const std::string &name);

    SectorFieldMapComponent();
    SectorFieldMapComponent(const SectorFieldMapComponent &right);
    ~SectorFieldMapComponent();

    /// Return field.
    //  The representation of the electro-magnetic field of the component
    //  (version for non-constant object).
    virtual EMField &getField() = 0;

    /// Return field.
    //  The representation of the electro-magnetic field of the component
    //  (version for constant object).
    virtual const EMField &getField() const = 0;

    virtual bool apply(const size_t &i, const double &t, double E[], double B[]) = 0;

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) = 0;

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) = 0;

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) = 0;

    virtual void finalise() = 0;

    virtual bool bends() const = 0;

    virtual void getDimensions(double &zBegin, double &zEnd) const = 0;

private:    
    
};


