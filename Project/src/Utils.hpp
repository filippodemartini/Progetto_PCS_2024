#pragma once

#include <iostream>
#include "GeometryLibrary.hpp"

using namespace std;

namespace GeometryLibrary{

//bool ImportDFN(const string &filepath, GeometryDFN& dfn);

bool ImportFractures(const string &filename, GeometryDFN& dfn);

//bool FindTraces();

Vector3d FindBarycentre(vector<Vector3d>& vertices);

// bool TIPS ??

}
