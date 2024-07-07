#include <iostream>
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

    FindTraces(dfn);

    FindTracesType(dfn);

    TracesLength(dfn);

    string fileOutputTracce = "./Tracce_FR10.txt";
    OutputTraces(dfn, fileOutputTracce);

    string fileOutputFratture = "./Fratture_FR10.txt";
    OutputFractures(dfn,fileOutputFratture);

  return 0;
}
