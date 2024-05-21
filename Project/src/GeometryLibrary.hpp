#pragma once

#include <iostream>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace GeometryLibrary {

struct GeometryDFN
{
    unsigned int Number_Fractures = 0;
    vector<unsigned int> Fractures_Id = {};
    vector<unsigned int> Fractures_Number_Vertices = {};
    map<unsigned int, vector<Vector3d>> Fractures_Vertices = {};

    unsigned int Number_Traces = 0;
    vector<unsigned int> Traces_Id = {};
    vector<Vector2i> Traces_Generator_Id;
    map<unsigned int, array<Vector3d,2>> Traces_Coordinates;
    vector<double> traces_length = {};
    bool Tips;

    double tol = 2*numeric_limits<double>::epsilon();
};
}

