#include "gtest/gtest.h"

#include "boost/lambda/lambda.hpp"
#include "boost/function.hpp"

#include <iostream>
#include <vector>
#include <algorithm>

#include "Problems.h"

using namespace std;
using namespace boost::lambda;

class LambdaTest : public ::testing::Test {
protected:
	vector<double> v;

	virtual void SetUp(){
		v.push_back( 1.0 );
		v.push_back( 2.0 );
		v.push_back( 3.0 );
	}
};




// Tests boost lambda functionality
TEST_F(LambdaTest, Alive) {
	for_each(v.begin(), v.end(), cout << _1 << '\n');
	std::for_each(v.begin(), v.end(), _1 = _1 * _1);
	for_each(v.begin(), v.end(), cout << _1 << '\n');
}

//Test that the Problems wrapping technique works
TEST_F(LambdaTest, WrapProblem) {
	for_each(v.begin(), v.end(), cout << _1 << '\n');

	//Problem * P = new Basin(1, 3);
	Problem::Interface * P = Problem::Factory("BASIN", 1, 3);
	cout << (P->Objectives[0]) (v) << endl;

}

class VerticalIterTest : public ::testing::Test {
public:
	template< typename T >
	class vertical_iterator : public T::const_iterator {
		int column;
	public:
		vertical_iterator( int c ) : column( c ) {}
		typename T::value_type::value_type operator*() const { return (this->T::const_iterator::operator*())[column]; }
		vertical_iterator& operator= (const typename T::const_iterator & rhs) { this->T::const_iterator::operator=(rhs); }
	};
	typedef std::vector< std::vector< double > > Data;
};

#include <algorithm>

// Tests vertical iterator concept
TEST_F(VerticalIterTest, Alive) {
	std::vector< std::vector< double > > data;
	std::vector< double > v(3);
	v[0] = 0.0; v[1] = 1.0; v[2] = 2.0;
	data.push_back(v);
	v[0]++ ; v[1]++ ; v[2]++ ;
	data.push_back(v);

	//This isn't supported since we havent implemented a
	//comparison operator (vi != data.end() fails)
	/*
	vertical_iterator< Data > vi(1);
	for(vi = data.begin(); vi != data.end(); vi++){
		printf("%lf\n", *vi);
	}
	*/

	vertical_iterator< Data > vi1(1);
	vertical_iterator< Data > vi2(1);
	vi1 = data.begin();
	vi2 = data.end();
	printf("%lf\n", *std::min_element(vi1, vi2));
}

//Apparently this will have to wait until compiler support
//becomes more widespread...
/*
#include <future>
#include <thread>
 
// Tests futures support

TEST_F(FuturesTest, Alive) {

    // future from a packaged_task
    std::packaged_task<int()> task([](){ return 7; }); // wrap the function
    std::future<int> f1 = task.get_future();  // get a future
    std::thread(std::move(task)).detach(); // launch on a thread
 
    // future from an async()
    std::future<int> f2 = std::async(std::launch::async, [](){ return 8; });
 
    // future from a promise
    std::promise<int> p;
    std::future<int> f3 = p.get_future();
    std::thread( [](std::promise<int>& p){ p.set_value(9); }, 
                 std::ref(p) ).detach();
 
    std::cout << "Waiting...";
    f1.wait();
    f2.wait();
    f3.wait();
    std::cout << "Done!\nResults are: "
              << f1.get() << ' ' << f2.get() << ' ' << f3.get() << '\n';
}
*/