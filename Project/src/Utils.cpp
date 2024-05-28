#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace GeometryLibrary{

bool ImportFractures(const string &filename, GeometryDFN& dfn)
{
    ifstream file;
    file.open(filename);

    if(file.fail())
        return false;

    string line;
     while(!file.eof())
     {
         getline(file,line);

         if(line[0] != '#')
             break;
    }

    istringstream convertNumFract;
    convertNumFract.str(line);
    convertNumFract >> dfn.Number_Fractures;

    dfn.Fractures_Id.reserve(dfn.Number_Fractures);
    dfn.Fractures_Number_Vertices.reserve(dfn.Number_Fractures);

    for(unsigned int i = 0; i<dfn.Number_Fractures; i++)
    {
        getline(file,line);
        getline(file,line);
        replace(line.begin(), line.end(), ';', ' ');
        istringstream convert(line);

        unsigned int Id;
        unsigned int Number_Vertices;
        convert >> Id >> Number_Vertices;

        dfn.Fractures_Id.push_back(Id);
        dfn.Fractures_Number_Vertices.push_back(Number_Vertices);

        getline(file,line);
        vector<Vector3d> Matrix_Vertices;
        Matrix_Vertices.resize(Number_Vertices);

        getline(file,line);
        vector<double> x_coordinates;
        istringstream convertX(line);
        string dataX;

        while(getline(convertX, dataX, ';')){
            x_coordinates.push_back(stod(dataX));
        }

        for(unsigned int j=0; j<Number_Vertices; j++){
            Matrix_Vertices[j][0] = x_coordinates[j];
        }

        getline(file,line);
        vector<double> y_coordinates;
        istringstream convertY(line);
        string dataY;

        while(getline(convertY, dataY, ';')){
            y_coordinates.push_back(stod(dataY));
        }

        for(unsigned int j=0; j<Number_Vertices; j++){
            Matrix_Vertices[j][1] = y_coordinates[j];
        }

        getline(file,line);
        vector<double> z_coordinates;
        istringstream convertZ(line);
        string dataZ;

        while(getline(convertZ, dataZ, ';')){
            z_coordinates.push_back(stod(dataZ));
        }

        for(unsigned int j=0; j<Number_Vertices; j++){
            Matrix_Vertices[j][2] = z_coordinates[j];
        }

        dfn.Fractures_Vertices[Id] = Matrix_Vertices;

    }
    file.close();

    for(unsigned int j = 0; j<dfn.Number_Fractures; j++){
        cout << "ID frattura: " << dfn.Fractures_Id[j] << endl;
        cout << "Numero vertici frattura: " << dfn.Fractures_Number_Vertices[j] << endl;

        cout << scientific << setprecision(16) << "Coordinate vertici: " << endl;
        for(unsigned int k = 0; k < dfn.Fractures_Number_Vertices[j]; k++){
            Vector3d vertice = dfn.Fractures_Vertices[dfn.Fractures_Id[j]][k];
            cout << "vertice " << k << " : " << vertice.transpose() << endl;
        }
    }

    return true;
}

inline Vector3d FindBarycentre(vector<Vector3d>& fracture)
{
    Vector3d sum_vertices = Vector3d::Zero();
    Vector3d barycentre = Vector3d::Zero();
    int num_vert = fracture.size();
    for(const auto& vertices : fracture){
        sum_vertices += vertices;
    }
    barycentre = sum_vertices/num_vert;

    return barycentre;
}

double CircleRadius(vector<Vector3d>& fracture)
{
    double radius = 0.0;
    for(const auto& vertices : fracture){
        Vector3d vector_distance = vertices - FindBarycentre(fracture);
        double distance = vector_distance.norm();
        if(distance>radius)
            radius = distance;
    }
    return radius;
}

inline Vector3d NormalToPlane(vector<Vector3d>& fracture)
{
    Vector3d p0 = fracture[0];
    Vector3d p1 = fracture[1];
    Vector3d p2 = fracture[2];

    Vector3d v = p1 - p0;
    Vector3d u = p2 - p0;

    Vector3d n = u.cross(v).normalized();
    return n;
}

bool FirstSelectionTraces(vector<Vector3d>& fracture_generator1, vector<Vector3d>& fracture_generator2, double tol)
{
    Vector3d barycentre1 = FindBarycentre(fracture_generator1);
    Vector3d barycentre2 = FindBarycentre(fracture_generator2);

    double barycentres_distance = (barycentre1 - barycentre2).squaredNorm();

    double radius1 = CircleRadius(fracture_generator1);
    double radius2 = CircleRadius(fracture_generator2);

    //unsigned int circle_intersections = 0;

    if ((radius1 + radius2 + tol) > barycentres_distance)
    {
        return true;
    }
    return false;
}


// CAMBIARE VOID
void FindTraces(GeometryDFN& dfn)
{
    unsigned int index_trace = 0;
    double tol = numeric_limits<double>::epsilon();
    // CONTROLLARE COSA CAMBIA SENZA IL -1; il -1 serve per l'id
    for(unsigned int i = 0; i < dfn.Number_Fractures - 1; i++){
        for(unsigned int j = i+1; j<dfn.Number_Fractures; j++){
            if(FirstSelectionTraces(dfn.Fractures_Vertices[i], dfn.Fractures_Vertices[j], tol)){
                Vector3d Normal_i = NormalToPlane(dfn.Fractures_Vertices[i]);
                Vector3d Normal_j = NormalToPlane(dfn.Fractures_Vertices[j]);
                Vector3d T = Normal_i.cross(Normal_j);
                Matrix3d Planes_Matrix;
                Planes_Matrix.row(0) = Normal_i;
                Planes_Matrix.row(1) = Normal_j;
                Planes_Matrix.row(2) = T;
                double d1 = Normal_i.dot(FindBarycentre(dfn.Fractures_Vertices[i]));
                double d2 = Normal_j.dot(FindBarycentre(dfn.Fractures_Vertices[j]));
                double d3 = 0;
                Vector3d b = {d1,d2,d3};

                double det_matrix = Planes_Matrix.determinant();
                if(det_matrix > tol){
                    Vector3d intersection = Planes_Matrix.fullPivLu().solve(b); // trovato vettore p1,p2,p3 del foglio geom comp 2
                }


            }
        }
    }
}

bool parallel_planes(vector<Vector3d>& fracture1, vector<Vector3d>& fracture2)
{
    bool result = false;
    double tol = numeric_limits<double>::epsilon();
    Vector3d normal1 = NormalToPlane(fracture1);
    Vector3d normal2 = NormalToPlane(fracture2);

    Vector3d direction = normal1.cross(normal2);

    if(direction.norm() < tol)
    {
        result = true;
    }
    return result;
}

inline MatrixXd fracture_vertices_line(unsigned int id_vertex1, unsigned int id_vertex2, const vector<Vector3d>& coordinates)
{
    double x1 = coordinates[id_vertex1][0];
    double y1 = coordinates[id_vertex1][1];
    double z1 = coordinates[id_vertex1][2];
    double x2 = coordinates[id_vertex2][0];
    double y2 = coordinates[id_vertex2][1];
    double z2 = coordinates[id_vertex2][2];

    Vector3d direction = {x2-x1, y2-y1, z2-z1};
    Vector3d line_origin = {x1, y1, z1};
    MatrixXd direction_and_line_origin;
    direction_and_line_origin.resize(2,3);
    direction_and_line_origin.row(0) = direction;
    direction_and_line_origin.row(1) = line_origin;

    return direction_and_line_origin;
}

inline Vector2d alpha_beta_intersection(MatrixXd fr_v_line, MatrixXd intersection)
{
    Vector3d direction1 = fr_v_line.row(0).transpose();
    Vector3d line_origin1 = fr_v_line.row(1).transpose();

    Vector3d direction2 = intersection.row(0).transpose();
    Vector3d line_origin2 = intersection.row(1).transpose();

    MatrixXd A = MatrixXd::Zero(3,2);
    A.col(0) = direction1;
    A.col(1) = -direction2;  // per tornare al valore positivo di beta dal -beta che abbiamo nella sottrazione tra equazioni delle due rette

    double b0 = line_origin2[0] - line_origin1[0];
    double b1 = line_origin2[1] - line_origin1[1];
    double b2 = line_origin2[2] - line_origin1[2];
    Vector3d b = {b0,b1,b2};

    Vector2d alpha_beta = A.householderQr().solve(b);

    return alpha_beta;
}

Matrix3d fracture_plane(const vector<Vector3d>& coordinates, vector<unsigned int> id_vertex)
{
    Matrix3d FP;
    Vector3d p0 = coordinates[id_vertex[0]];
    Vector3d p1 = coordinates[id_vertex[1]];
    Vector3d p2 = coordinates[id_vertex[2]];

    FP.row(0) = p0;
    FP.row(1) = p2-p0;
    FP.row(2) = p1-p0;

    return FP;
}

bool point_on_line(Vector3d& line_origin, Vector3d& line_end, Vector3d& point)
{ // se il punto è più distante dall'inizio del segmento rispetto alla fine del segmento stesso allora il punto non appartiene al segmento
    bool on_line = false;

    double tol = numeric_limits<double>::epsilon();
    double line_length = (line_origin - line_end).norm();
    double distance_point_origin = (line_origin - point).norm();
    double distance_point_end = (line_end - point).norm();
    double difference = distance_point_origin + distance_point_end - line_length;

    if(abs(difference)<tol)
    {
        on_line = true;
    }

    return on_line;
}

//Vector3d point_intersection_lines(GeometryDFN& dfn)

}
