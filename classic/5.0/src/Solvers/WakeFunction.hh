#ifndef WAKEFUNCTION_HH
#define WAKEFUNCTION_HH

class ElementBase;
class PartBunch;

class WakeFunction
{
public:
    WakeFunction(string name, ElementBase *elref);
    virtual void apply(PartBunch &bunch) = 0;
    virtual const string getType() const = 0;
    const string getName();
    void updateElement(ElementBase *newref);
protected:
    ElementBase *element_ref_m;

private:
    const string name_m;
};

inline WakeFunction::WakeFunction(string name, ElementBase *elref): 
    name_m(name),
    element_ref_m(elref)
{}

inline const string WakeFunction::getName()
{
    return name_m;
}

class LineDensity: public vector<double>
{
public:
    LineDensity(int size = 0, double defaultValue = 0.0) : vector<double>(size,defaultValue) {}
    void getFirstDerivative(vector<double> &firstDerivative, const double &hz);
};

inline void LineDensity::getFirstDerivative(vector<double> &firstDerivative, const double &hz)
{
    const int size = this->size();
    if (firstDerivative.size() != size)
        firstDerivative.resize(size,0.0);

    firstDerivative[0] = ((*this)[1] - (*this)[0])/hz;
    for (int i = 1; i < size - 1; ++i)
        firstDerivative[i] = ((*this)[i + 1] - (*this)[i - 1])/hz;
    firstDerivative[size - 1] = ((*this)[size - 1] - (*this)[size - 2])/hz;
}

inline void WakeFunction::updateElement(ElementBase *newref)
{
    element_ref_m = newref;
}

#endif // WAKEFUNCTION_HH
