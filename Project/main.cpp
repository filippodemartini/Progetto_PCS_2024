#include <iostream>
#include <vector>
#include "GeometryLibrary.hpp"
#include "Utils.hpp"

using namespace std;
using namespace Eigen;
using namespace GeometryLibrary;

int main()
{
    GeometryDFN dfn;
    string filename = "./FR_data.txt";
    if(!ImportFractures(filename, dfn)){
        return 1;
    }

    cout << endl;
    for(unsigned int i = 0; i<dfn.Number_Fractures; i++){
    Vector3d barycentre = FindBarycentre(dfn.Fractures_Vertices[i]);

        cout << "Barycentre : (" << barycentre[0] << ";" << barycentre[1] << ";" << barycentre[2] << ")" << endl;
    }

    bool int01 = FirstSelectionTraces(dfn.Fractures_Vertices[0], dfn.Fractures_Vertices[1],2*numeric_limits<double>::epsilon());
    bool int02 = FirstSelectionTraces(dfn.Fractures_Vertices[0], dfn.Fractures_Vertices[2],2*numeric_limits<double>::epsilon());
    bool int12 = FirstSelectionTraces(dfn.Fractures_Vertices[1], dfn.Fractures_Vertices[2],2*numeric_limits<double>::epsilon());

    cout << int01 << "   " << int02 << "   " << int12 << endl;

  return 0;
}
