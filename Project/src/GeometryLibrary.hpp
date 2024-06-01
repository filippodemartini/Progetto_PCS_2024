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
   // vector<unsigned int> Vertices_id = {};

    unsigned int Number_Traces = 0;
    vector<unsigned int> Traces_Id = {};
    vector<Vector2i> Traces_Generator_Id;
   // map<unsigned int, vector<Vector2i>> Id_traces_generator;
    map<unsigned int, array<Vector3d,2>> Traces_Coordinates;
    vector<double> traces_length = {};
    map<unsigned int, array<bool,2>> Tips;

    //double tol = numeric_limits<double>::epsilon();
};
}

