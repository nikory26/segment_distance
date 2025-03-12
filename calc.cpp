#include <iostream>
#include "geometry.h"

int main()
{
    real data[4 * SpaceDim];
    using namespace std;
    cout << "Input coordinates for segment 1 begin: ";
    cin >> data[0] >> data[1] >> data[2];
    cout << "Input coordinates for segment 1 end: ";
    cin >> data[3] >> data[4] >> data[5];
    cout << "Input coordinates for segment 2 begin: ";
    cin >> data[6] >> data[7] >> data[8];
    cout << "Input coordinates for segment 2 end: ";
    cin >> data[9] >> data[10] >> data[11];
    Segment3D s1(data,data+3);
    Segment3D s2(data+6,data+9);
    cout<<"Segment 1: "<<s1<<endl;
    cout<<"Segment 2: "<<s2<<endl;
    cout<<"Distance : "<<segments_distance(s1,s2)<<endl;
    return 0;
}