#ifndef WAKEFUNCTION_HH
#define WAKEFUNCTION_HH

class ElementBase;
class PartBunch;

class WakeFunction {
public:
    WakeFunction(std::string name, ElementBase *elref);
    virtual ~WakeFunction(){ };
    virtual void apply(PartBunch &bunch) = 0;
    virtual const std::string getType() const = 0;
    const std::string getName();
    void updateElement(ElementBase *newref);
protected:
    ElementBase *element_ref_m;

private:
    const std::string name_m;
};

inline WakeFunction::WakeFunction(std::string name, ElementBase *elref):
    element_ref_m(elref),
    name_m(name)
{}

inline const std::string WakeFunction::getName() {
    return name_m;
}

class LineDensity: public std::vector<double> {
public:
    LineDensity(int size = 0, double defaultValue = 0.0) : std::vector<double>(size, defaultValue) {}
    void getFirstDerivative(std::vector<double> &firstDerivative, const double &hz);
};

inline void LineDensity::getFirstDerivative(std::vector<double> &firstDerivative, const double &hz) {
    const size_t size = this->size();
    if(firstDerivative.size() != size)
        firstDerivative.resize(size, 0.0);

    firstDerivative[0] = ((*this)[1] - (*this)[0]) / hz;
    for(unsigned int i = 1; i + 1 < size; ++i)
        firstDerivative[i] = ((*this)[i + 1] - (*this)[i - 1]) / hz;
    firstDerivative[size - 1] = ((*this)[size - 1] - (*this)[size - 2]) / hz;
}

inline void WakeFunction::updateElement(ElementBase *newref) {
    element_ref_m = newref;
}

#endif // WAKEFUNCTION_HH
