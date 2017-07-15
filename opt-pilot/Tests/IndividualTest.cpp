#include "Optimizer/EA/Individual.h"
#include "Optimizer/EA/Individual.h"
#include "gtest/gtest.h"

#include "boost/smart_ptr.hpp"

namespace {

    // The fixture for testing class Foo.
    class IndividualTest : public ::testing::Test {
    protected:

        IndividualTest() {
            // You can do set-up work for each test here.
            lower_bound_ = -10.0;
            upper_bound_ = 19.79;
        }

        virtual ~IndividualTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:

        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        boost::shared_ptr<Individual> createIndividual(size_t num_genes) {

            Individual::bounds_t bounds;
            for(size_t i=0; i < num_genes; i++)
                bounds.push_back(
                    std::pair<double, double>(lower_bound_, upper_bound_));

            boost::shared_ptr<Individual> ind(new Individual(bounds));
            return ind;
        }

        double lower_bound_;
        double upper_bound_;

    };

    TEST_F(IndividualTest, IndividualRespectsBounds) {

        size_t num_genes = 1;
        boost::shared_ptr<Individual> ind = createIndividual(num_genes);
        double gene = ind->genes[0];

        EXPECT_LE(lower_bound_, gene) << "gene should respect lower bound";
        EXPECT_GE(upper_bound_, gene) << "gene should respect upper bound";

        size_t my_size = ind->genes.size();
        EXPECT_EQ(static_cast<size_t>(num_genes), my_size)
            << "individual should only have " << num_genes << " gene(s)";

    }

    TEST_F(IndividualTest, IndividualHasCorrectNumberOfGenes) {

        size_t num_genes = 12;
        boost::shared_ptr<Individual> ind = createIndividual(num_genes);

        size_t my_size = ind->genes.size();
        EXPECT_EQ(static_cast<size_t>(num_genes), my_size)
            << "individual should only have " << num_genes << " gene(s)";

    }

    TEST_F(IndividualTest, IndividualRandomGene) {

        size_t num_genes = 1;
        boost::shared_ptr<Individual> ind = createIndividual(num_genes);
        double gene = ind->genes[0];
        double new_gene = ind->new_gene(0);

        EXPECT_NE(gene, new_gene) << "new gene should be different";
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
