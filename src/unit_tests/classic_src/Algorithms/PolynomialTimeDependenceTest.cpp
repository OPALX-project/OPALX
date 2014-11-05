#include "gtest/gtest.h"
#include "Utilities/GeneralClassicException.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "Algorithms/PolynomialTimeDependence.h"

TEST(PolynomialTimeDependenceTest, PolynomialTimeDependenceTest) {
    // Check empty polynomial coefficients always returns 0.
    std::vector<double> test;
    PolynomialTimeDependence time_dependence_1(test);
    EXPECT_DOUBLE_EQ(time_dependence_1.getValue(0.1), 0.);

    // Check constant term produces constant
    test.push_back(1.);
    PolynomialTimeDependence time_dependence_2(test);
    EXPECT_DOUBLE_EQ(time_dependence_2.getValue(0.1), 1.);

    // Check cubic terms
    test.push_back(2.);
    test.push_back(3.);
    PolynomialTimeDependence time_dependence_3(test);
    EXPECT_DOUBLE_EQ(time_dependence_3.getValue(0.1), 1.23);

    // Check clone produces same result
    PolynomialTimeDependence* time_dependence_clone = time_dependence_3.clone();
    EXPECT_DOUBLE_EQ(time_dependence_clone->getValue(0.1), 1.23);
    delete time_dependence_clone;
}

TEST(PolynomialTimeDependenceTest, TDMapTest) {
    // throw on empty value
    EXPECT_THROW(AbstractTimeDependence::getTimeDependence("name"),
                 GeneralClassicException);
    std::vector<double> test;

    // set/get time dependence
    PolynomialTimeDependence time_dep(test);
    std::shared_ptr<PolynomialTimeDependence> td1(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td1", td1);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td1);
    std::shared_ptr<PolynomialTimeDependence> td2(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td2", td2);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td2"), td2);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td1);
    // set time dependence overwriting existing time dependence
    // should overwrite, without memory leak
    std::shared_ptr<PolynomialTimeDependence> td3(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td1", td3);
    EXPECT_EQ(AbstractTimeDependence::getTimeDependence("td1"), td3);
}

TEST(PolynomialTimeDependenceTest, TDMapNameLookupTest) {
    EXPECT_THROW(AbstractTimeDependence::getName(NULL),
                 GeneralClassicException);
    PolynomialTimeDependence time_dep(std::vector<double>(1, 1));
    std::shared_ptr<PolynomialTimeDependence> td1(time_dep.clone());
    std::shared_ptr<PolynomialTimeDependence> td2(time_dep.clone());
    std::shared_ptr<PolynomialTimeDependence> td3(time_dep.clone());
    AbstractTimeDependence::setTimeDependence("td1", td1);
    AbstractTimeDependence::setTimeDependence("td2", td2);
    AbstractTimeDependence::setTimeDependence("td3", td2);
    std::string name1 = AbstractTimeDependence::getName(td1);
    EXPECT_EQ(name1, "td1");
    std::string name2 = AbstractTimeDependence::getName(td2);
    EXPECT_TRUE(name2 == "td2" || name2 == "td3");
    EXPECT_THROW(AbstractTimeDependence::getName(td3),
                 GeneralClassicException);

}
