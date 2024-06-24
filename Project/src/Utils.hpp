#pragma once

#include <iostream>
#include "GeometryLibrary.hpp"

using namespace std;

namespace GeometryLibrary{

// bool ImportDFN(const string &filepath, GeometryDFN& dfn);

bool ImportFractures(const string &filename, GeometryDFN& dfn);

inline Vector3d FindBarycentre(vector<Vector3d>& vertices);

double CircleRadius(vector<Vector3d>& vertices);

inline Vector3d NormalToPlane(vector<Vector3d>& fracture);

bool FirstSelectionTraces(vector<Vector3d>& fracture_generator1, vector<Vector3d>& fracture_generator2, double tol);

void FindTraces(GeometryDFN& dfn);

void TracesType(GeometryDFN& dfn);

bool point_on_line(Vector3d& line_origin, Vector3d& line_end, Vector3d& point_coord);

bool check_barycentre_coord(const Vector3d& p0, const Vector3d& p1, const Vector3d& p2, const Vector3d& p3);

bool check_inside_fracture(const Vector3d& point, vector<Vector3d>& fracture_vertex);

bool parallel_planes(vector<Vector3d>& fracture1, vector<Vector3d>& fracture2);

void calcolaTipologiaTracce(GeometryDFN& DFN);

map<unsigned int,array<bool, 2>> riordinaTracce(const vector<double>& length, map<unsigned int,array<bool, 2>>& type, vector<unsigned int>& trace_id);

void calcolaLunghezzaTracce(GeometryDFN& DFN);

bool OutputTracce(const GeometryDFN& DFN, const string& fileOutput);

bool OutputFratture(const GeometryDFN& DFN, const string& fileOutput);

//void lengthTraces(GeometryDFN& DFN);

//bool isFratturaIntersecante(unsigned int k, unsigned int i, const GeometryDFN& DFN);

//std::array<bool, 2> getTipologia(unsigned int i, unsigned int k, GeometryDFN& DFN);

//inline MatrixXd fracture_vertices_line(unsigned int id_vertex1, unsigned int id_vertex2, const vector<Vector3d>& coordinates);

//inline Vector2d alpha_beta_intersection(MatrixXd fr_v_line, MatrixXd intersection);

//inline Vector3d point_intersection_lines(GeometryDFN& dfn);

bool checkSegmentIntersection(vector<Vector3d>& intersections, const Vector3d normal_plane, Vector3d plane_point, Vector3d point1, Vector3d point2, double tol);

// bool TIPS ??

}
