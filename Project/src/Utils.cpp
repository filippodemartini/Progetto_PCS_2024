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

Vector3d FindBarycentre(vector<Vector3d>& fracture)
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

Vector3d NormalToPlane(vector<Vector3d>& fracture)
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

    if ((radius1 + radius2 + tol) > barycentres_distance)
    {
        return true;
    }
    return false;
}

// void FindTraces(GeometryDFN& dfn)
// {
//     const unsigned int n = 1000;

//     dfn.Traces_Id.reserve(n);
//     dfn.Traces_Generator_Id.reserve(n);


// }

}
