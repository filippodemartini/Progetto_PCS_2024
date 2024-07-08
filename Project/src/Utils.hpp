#pragma once

#include "GeometryLibrary.hpp"

using namespace std;

namespace GeometryLibrary{

bool ImportFractures(const string &filename, GeometryDFN& dfn);

inline Vector3d FindBarycentre(vector<Vector3d>& vertices);

double CircleRadius(vector<Vector3d>& vertices);

inline Vector3d NormalToPlane(vector<Vector3d>& fracture);

bool FirstSelectionTraces(vector<Vector3d>& fracture_generator1, vector<Vector3d>& fracture_generator2, double tol);

void FindTraces(GeometryDFN& dfn);

bool point_on_line(Vector3d& line_origin, Vector3d& line_end, Vector3d& point_coord);

bool check_barycentre_coord(const Vector3d& p0, const Vector3d& p1, const Vector3d& p2, const Vector3d& p3);

bool check_inside_fracture(const Vector3d& point, vector<Vector3d>& fracture_vertex);

bool parallel_planes(vector<Vector3d>& fracture1, vector<Vector3d>& fracture2);

void FindTracesType(GeometryDFN& DFN);

void TracesLength(GeometryDFN& DFN);

vector<unsigned int> reorganiseLength(const vector<unsigned int>& traceIds, const GeometryDFN& DFN);

bool OutputTraces(const GeometryDFN& DFN, const string& fileOutput);

bool OutputFractures(const GeometryDFN& DFN, const string& fileOutput);

}
