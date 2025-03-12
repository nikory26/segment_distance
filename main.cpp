#include <iostream>
#include <iomanip>
#include "geometry.h"
#include "tests.h"

int main()
{
    // Point3D p = {1,1,1};
    // Segment3D s1({1,1,0},{3,0,1});
    // Segment3D s2({10,20,30},{5,4,3});
    // std::cout << "The distance btw segments is " << segments_distance(s1,s2);
    int res = 0;
    res += test1();
    res += test2();
    res += test3();
    res += test4();
    res += test5();
    res += test6();
    res += test7();
    return res;
}