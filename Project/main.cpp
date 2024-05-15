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
    string filename = "./FR10_data.txt";
    if(!ImportFractures(filename, dfn)){
        return 1;
    }

    cout << endl;
    for(unsigned int i = 0; i<dfn.Number_Fractures; i++){
    Vector3d barycentre = FindBarycentre(dfn.Fractures_Vertices[i]);

        cout << "Barycentre : (" << barycentre[0] << ";" << barycentre[1] << ";" << barycentre[2] << ")" << endl;
    }

  return 0;
}
