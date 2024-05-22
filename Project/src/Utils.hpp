#pragma once

#include <iostream>
#include "GeometryLibrary.hpp"

using namespace std;

namespace GeometryLibrary{

bool ImportDFN(const string &filepath, GeometryDFN& dfn);

bool ImportFractures(const string &filename, GeometryDFN& dfn);

Vector3d FindBarycentre(vector<Vector3d>& vertices);

double CircleRadius(vector<Vector3d>& vertices);

Vector3d NormalToPlane(vector<Vector3d>& fracture);

bool FirstSelectionTraces(vector<Vector3d>& fracture_generator1, vector<Vector3d>& fracture_generator2, double tol);

void FindTraces(GeometryDFN& dfn);

bool parallel_planes(vector<Vector3d>& fracture1, vector<Vector3d>& fracture2);

Vector3d point_intersection_lines(GeometryDFN& dfn);

bool checkSegmentIntersection(vector<Vector3d>& intersections, const Vector3d normal_plane, Vector3d plane_point, Vector3d point1, Vector3d point2, double tol);

// bool TIPS ??

}
