#include "tests.h"

int test1()
{
    //Overlapping parralel segments test
    Segment3D s1({0,0,0},{1,1,0});
    Segment3D s2({-0.5,-0.5,1},{0.5,0.5,1});
    //The distance betweemn these segments should be equal to 1.0
    real d = segments_distance(s1,s2);
    const real eps = 1000 * std::numeric_limits<real>::epsilon();
    real d1 = segments_distance(s2,s1);
    assert(fabs(d-d1) < eps);
    real err = fabs(d - 1);
    bool test_res = err <= eps;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);   
    std::cout << std::endl << "TEST1: Overlapping parallel segments" << std::endl;
    std::cout << "Calculating distanse between segments " << s1 << " and " << s2 << std::endl;
    std::cout << "The calculated distanse between segments is equal to " << d << std::endl;
    std::cout << "The reference distanse between segments is equal to " << 1.0 << std::endl;
    std::cout << "The absolute difference is equal to " << err << std::endl;
    std::cout << "The test1 was passed: " << (test_res ? "YES" : "NO") << std::endl << std::endl;
    return !test_res;
}

int test2()
{
    //Non-overlapping parralel segments test
    Segment3D s1({0,0,0},{1,1,0});
    Segment3D s2({-5,-5,1},{-0.5,-0.5,1});
    //The distance betweemn these segments should be equal to 1.0
    real d = segments_distance(s1,s2);
    real ref_d = s2.end().norm();
    const real eps = 1000 * std::numeric_limits<real>::epsilon();
    real d1 = segments_distance(s2,s1);
    assert(fabs(d-d1) < eps);
    real err = fabs(d - ref_d);
    bool test_res = err <= eps;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);   
    std::cout << std::endl << "TEST2: Non-overlapping parallel segments" << std::endl;
    std::cout << "Calculating distanse between segments " << s1 << " and " << s2 << std::endl;
    std::cout << "The calculated distanse between segments is equal to " << d << std::endl;
    std::cout << "The reference distanse between segments is equal to " << ref_d << std::endl;
    std::cout << "The absolute difference is equal to " << err << std::endl;
    std::cout << "The test2 was passed: " << (test_res ? "YES" : "NO") << std::endl << std::endl;
    return !test_res;
}


int test3()
{
    //Skew segments
    Segment3D s1({0,0,0},{1,1,0});
    Segment3D s2({1,0,1},{0,1,1});
    //The distance betweemn these segments should be equal to 1.0
    real d = segments_distance(s1,s2);
    real ref_d = distance({0.5,0.5,0},{0.5,0.5,1});
    const real eps = 1000 * std::numeric_limits<real>::epsilon();
    real d1 = segments_distance(s2,s1);
    assert(fabs(d-d1) < eps);
    real err = fabs(d - ref_d);
    bool test_res = err <= eps;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);   
    std::cout << std::endl << "TEST3: Skew segments" << std::endl;
    std::cout << "Calculating distanse between segments " << s1 << " and " << s2 << std::endl;
    std::cout << "The calculated distanse between segments is equal to " << d << std::endl;
    std::cout << "The reference distanse between segments is equal to " << ref_d << std::endl;
    std::cout << "The absolute difference is equal to " << err << std::endl;
    std::cout << "The test3 was passed: " << (test_res ? "YES" : "NO") << std::endl << std::endl;
    return !test_res;
}

int test4()
{
    //Skew segments, but not intersected in projection
    Segment3D s1({0,0,0},{1,1,0});
    Segment3D s2({1,0,1},{0.7,0.3,1});
    //The distance betweemn these segments should be equal to 1.0
    real d = segments_distance(s1,s2);
    real ref_d = distance_to_line(s2.end(), s1);
    const real eps = 1000 * std::numeric_limits<real>::epsilon();
    real d1 = segments_distance(s2,s1);
    assert(fabs(d-d1) < eps);
    real err = fabs(d - ref_d);
    bool test_res = err <= eps;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);   
    std::cout << std::endl << "TEST4: Skew segments" << std::endl;
    std::cout << "Calculating distanse between segments " << s1 << " and " << s2 << std::endl;
    std::cout << "The calculated distanse between segments is equal to " << d << std::endl;
    std::cout << "The reference distanse between segments is equal to " << ref_d << std::endl;
    std::cout << "The absolute difference is equal to " << err << std::endl;
    std::cout << "The test4 was passed: " << (test_res ? "YES" : "NO") << std::endl << std::endl;
    return !test_res;
}

int test5()
{
    //Skew segments, but not intersected in projection
    Segment3D s1({0,0,0},{1,1,0});
    Segment3D s2({2,0,1},{1.7,0.3,1});
    //The distance betweemn these segments should be equal to 1.0
    real d = segments_distance(s1,s2);
    real ref_d = distance(s1.end(), s2.end());
    const real eps = 1000 * std::numeric_limits<real>::epsilon();
    real d1 = segments_distance(s2,s1);
    assert(fabs(d-d1) < eps);
    real err = fabs(d - ref_d);
    bool test_res = err <= eps;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);   
    std::cout << std::endl << "TEST5: Skew segments" << std::endl;
    std::cout << "Calculating distanse between segments " << s1 << " and " << s2 << std::endl;
    std::cout << "The calculated distanse between segments is equal to " << d << std::endl;
    std::cout << "The reference distanse between segments is equal to " << ref_d << std::endl;
    std::cout << "The absolute difference is equal to " << err << std::endl;
    std::cout << "The test5 was passed: " << (test_res ? "YES" : "NO") << std::endl << std::endl;
    return !test_res;
}

int test6()
{
    //Skew segments
    Segment3D s1({0,0,0},{1,1,0});
    Segment3D s2({0,1,1},{0.3,0.7,1});
    //The distance betweemn these segments should be equal to 1.0
    real d = segments_distance(s1,s2);
    real ref_d = distance_to_line(s2.end(),s1);
    const real eps = 1000 * std::numeric_limits<real>::epsilon();
    real d1 = segments_distance(s2,s1);
    assert(fabs(d-d1) < eps);
    real err = fabs(d - ref_d);
    bool test_res = err <= eps;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);   
    std::cout << std::endl << "TEST6: Skew segments" << std::endl;
    std::cout << "Calculating distanse between segments " << s1 << " and " << s2 << std::endl;
    std::cout << "The calculated distanse between segments is equal to " << d << std::endl;
    std::cout << "The reference distanse between segments is equal to " << ref_d << std::endl;
    std::cout << "The absolute difference is equal to " << err << std::endl;
    std::cout << "The test6 was passed: " << (test_res ? "YES" : "NO") << std::endl << std::endl;
    return !test_res;
}

int test7()
{
    //Skew segments
    Segment3D s1({0,0,0},{1,1,0});
    Segment3D s2({-2,1,1},{-1.7,0.7,1});
    //The distance betweemn these segments should be equal to 1.0
    real d = segments_distance(s1,s2);
    real ref_d = distance(s2.end(),s1.begin());
    const real eps = 1000 * std::numeric_limits<real>::epsilon();
    real d1 = segments_distance(s2,s1);
    assert(fabs(d-d1) < eps);
    real err = fabs(d - ref_d);
    bool test_res = err <= eps;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(3);   
    std::cout << std::endl << "TEST7: Skew segments" << std::endl;
    std::cout << "Calculating distanse between segments " << s1 << " and " << s2 << std::endl;
    std::cout << "The calculated distanse between segments is equal to " << d << std::endl;
    std::cout << "The reference distanse between segments is equal to " << ref_d << std::endl;
    std::cout << "The absolute difference is equal to " << err << std::endl;
    std::cout << "The test7 was passed: " << (test_res ? "YES" : "NO") << std::endl << std::endl;
    return !test_res;
}
