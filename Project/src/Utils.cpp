#include "Utils.hpp"
#include <iostream>
#include <array>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "MergeSort.hpp"

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
        return false;
    }
    return true;
}                  // CONTROLLA I RETURN DI FALSE E TRUE PER PRIMA SCREMATURA

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


// CAMBIARE VOID
void FindTraces(GeometryDFN& dfn)
{
    unsigned int index_trace = 0;
    double tol = numeric_limits<double>::epsilon();

    // CONTROLLARE COSA CAMBIA SENZA IL -1; il -1 serve per l'id
    for(unsigned int i = 0; i < dfn.Number_Fractures - 1; i++){
        for(unsigned int j = i+1; j<dfn.Number_Fractures; j++){

            if(parallel_planes(dfn.Fractures_Vertices[i], dfn.Fractures_Vertices[j])){
                cout << "Le fratture: " << dfn.Fractures_Id[i] << " e " << dfn.Fractures_Id[j] << " sono parallele" << endl;
                continue;           // CONTROLLA LA FUNZIONALITA DI CONTINUE
            }

            if(!FirstSelectionTraces(dfn.Fractures_Vertices[i], dfn.Fractures_Vertices[j], tol)){   // CONTROLLA IL !
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
                    Vector3d p_intersection = Planes_Matrix.fullPivLu().solve(b); // trovato vettore p1,p2,p3 del foglio geom comp 2
                    unsigned int num_traces_points = 0;
                    array<Vector3d,2> traces_points;
                    for(unsigned int k = 0; k < dfn.Fractures_Number_Vertices[i]; k++)
                    {                  // SE NON FUNZIONA SICURAMENTE QUA QUALCOSA NON VA BENE; CAPIRE SE QUESTO FOR DEVE GIRARE SU K O SU J

                        Vector3d v1 = dfn.Fractures_Vertices[i][k];
                        Vector3d v2 = dfn.Fractures_Vertices[i][(k + 1) % dfn.Fractures_Number_Vertices[i]]; // faccio sta roba per non andare fuori scope

                        MatrixXd A = MatrixXd::Zero(3,2);

                        A.col(0) = v2-v1;   // è una retta così come t
                        A.col(1) = T;  // per tornare al valore positivo di beta dal -beta che abbiamo nella sottrazione tra equazioni delle due rette

                        if((v2-v1).cross(T).squaredNorm() < tol)
                        {
                            continue;
                        }


                        double b0 = p_intersection[0] - v1[0];
                        double b1 = p_intersection[1] - v1[1];
                        double b2 = p_intersection[2] - v1[2];
                        Vector3d b = {b0,b1,b2};

                        Vector2d alpha_beta = A.householderQr().solve(b);

                        if(alpha_beta[0]>-tol && alpha_beta[0]<1+tol && num_traces_points < 2)
                        {//controllo che non prenda due volte lo stesso punto
                            Vector3d Intersection = v1+(alpha_beta[0]*(v2-v1));
                            if (num_traces_points == 0 && check_inside_fracture(Intersection, dfn.Fractures_Vertices[j]))
                                {
                                    traces_points[num_traces_points] = Intersection;
                                    num_traces_points += 1;
                                }
                            else if (num_traces_points == 1 && check_inside_fracture(Intersection, dfn.Fractures_Vertices[j]) && traces_points[0] !=Intersection)
                            { traces_points[num_traces_points]=Intersection;
                                num_traces_points+=1;
                                    dfn.Number_Traces+=1;
                                dfn.Traces_Id.push_back(index_trace);
                                Vector2i fracture = {dfn.Fractures_Id[i], dfn.Fractures_Id[j]};
                                dfn.Traces_Generator_Id.push_back(fracture);
                                dfn.Traces_Coordinates[index_trace] = traces_points;
                                index_trace += 1;
                            }

                        }

                    } // qui controllo sulla seconda frattura, almeno da trovare anche il secondo punto di intersezione nel caso rosso-blu in FR3
                    for(unsigned int k = 0; k < dfn.Fractures_Number_Vertices[j]; k++)
                    {                  // SE NON FUNZIONA SICURAMENTE QUA QUALCOSA NON VA BENE; CAPIRE SE QUESTO FOR DEVE GIRARE SU K O SU J

                        Vector3d v1 = dfn.Fractures_Vertices[j][k];
                        Vector3d v2 = dfn.Fractures_Vertices[j][(k + 1) % dfn.Fractures_Number_Vertices[i]];

                        MatrixXd A = MatrixXd::Zero(3,2);

                        A.col(0) = v2-v1;
                        A.col(1) = T;  // per tornare al valore positivo di beta dal -beta che abbiamo nella sottrazione tra equazioni delle due rette

                        if((v2-v1).cross(T).squaredNorm() < tol)
                        {
                            continue;
                        }


                        double b0 = p_intersection[0] - v1[0];
                        double b1 = p_intersection[1] - v1[1];
                        double b2 = p_intersection[2] - v1[2];
                        Vector3d b = {b0,b1,b2};

                        Vector2d alpha_beta = A.householderQr().solve(b);

                        if(alpha_beta[0]>-tol && alpha_beta[0]<1+tol && num_traces_points < 2)
                        {
                            Vector3d Intersection = v1+(alpha_beta[0]*(v2-v1));
                            if (num_traces_points == 0 && check_inside_fracture(Intersection, dfn.Fractures_Vertices[i]))
                            {
                                traces_points[num_traces_points] = Intersection;
                                num_traces_points += 1;
                            }
                            else if (num_traces_points == 1 && check_inside_fracture(Intersection, dfn.Fractures_Vertices[i]) && traces_points[0] !=Intersection)
                            { traces_points[num_traces_points]=Intersection;
                                num_traces_points+=1;
                                dfn.Number_Traces+=1;
                                dfn.Traces_Id.push_back(index_trace);
                                Vector2i fracture = {dfn.Fractures_Id[i], dfn.Fractures_Id[j]};
                                dfn.Traces_Generator_Id.push_back(fracture);
                                dfn.Traces_Coordinates[index_trace] = traces_points;
                                index_trace += 1;
                            }

                        }

                    }
                }
            }
        }
    }
}


// inline MatrixXd fracture_vertices_line(unsigned int id_vertex1, unsigned int id_vertex2, const vector<Vector3d>& coordinates)
// {
//     double x1 = coordinates[id_vertex1][0];
//     double y1 = coordinates[id_vertex1][1];
//     double z1 = coordinates[id_vertex1][2];
//     double x2 = coordinates[id_vertex2][0];
//     double y2 = coordinates[id_vertex2][1];
//     double z2 = coordinates[id_vertex2][2];

//     Vector3d direction = {x2-x1, y2-y1, z2-z1};
//     Vector3d line_origin = {x1, y1, z1};
//     MatrixXd direction_and_line_origin;
//     direction_and_line_origin.resize(2,3);
//     direction_and_line_origin.row(0) = direction;
//     direction_and_line_origin.row(1) = line_origin;

//     return direction_and_line_origin;
// }

// inline Vector2d alpha_beta_intersection(MatrixXd fr_v_line, MatrixXd intersection)
// {
//     double tol = numeric_limits<double>::epsilon();
//     Vector3d direction1 = fr_v_line.row(0).transpose();
//     Vector3d line_origin1 = fr_v_line.row(1).transpose();

//     Vector3d direction2 = intersection.row(0).transpose();
//     Vector3d line_origin2 = intersection.row(1).transpose();

//     MatrixXd A = MatrixXd::Zero(3,2);
//     A.col(0) = direction1;
//     A.col(1) = -direction2;  // per tornare al valore positivo di beta dal -beta che abbiamo nella sottrazione tra equazioni delle due rette

//     double b0 = line_origin2[0] - line_origin1[0];
//     double b1 = line_origin2[1] - line_origin1[1];
//     double b2 = line_origin2[2] - line_origin1[2];
//     Vector3d b = {b0,b1,b2};

//     if(direction1.cross(direction2).squaredNorm() < tol)
//     {
//         cout << "Le due rette direzionali sono parallele" << endl;
//     }

//     else
//     {
//        Vector2d alpha_beta = A.householderQr().solve(b);
//        return alpha_beta;
//     }
// }  //  CONTROLLA IL RETURN

// Matrix3d fracture_plane(const vector<Vector3d>& coordinates, vector<unsigned int> id_vertex)
// {
//     Matrix3d FP;
//     Vector3d p0 = coordinates[id_vertex[0]];
//     Vector3d p1 = coordinates[id_vertex[1]];
//     Vector3d p2 = coordinates[id_vertex[2]];

//     FP.row(0) = p0;
//     FP.row(1) = p2-p0;
//     FP.row(2) = p1-p0;

//     return FP;
// }


bool point_on_line(Vector3d& p1, Vector3d& p2, Vector3d& p3)
{ // se il punto è più distante dall'inizio del segmento rispetto alla fine del segmento stesso allora il punto non appartiene al segmento
    //bool on_line = false;

    // QUELLO CHE AVEVAMO PRIMA
    // double tol = numeric_limits<double>::epsilon();
    // double line_length = (line_origin - line_end).norm();
    // double distance_point_origin = (line_origin - point).norm();
    // double distance_point_end = (line_end - point).norm();

    Vector3d prodotto_Vettoriale = (p2-p1).cross(p3-p1);
    // Controllo se il punto p3 è sulla stessa retta dei punti p1 e p2 che rappresentano due vertici consecutivi della Frattura
    // Ovvero controllo che i tre punti siano collineari
    if (fabs(prodotto_Vettoriale[0]) < 1e-12 && fabs(prodotto_Vettoriale[1]) < 1e-12 && fabs(prodotto_Vettoriale[2]) < 1e-12)
    {
        return ((p3[0] >= min(p1[0], p2[0]) && p3[0] <= max(p1[0], p2[0]) &&
                 p3[1] >= min(p1[1], p2[1]) && p3[1] <= max(p1[1], p2[1]) &&
                 p3[2] >= min(p1[2], p2[2]) && p3[2] <= max(p1[2], p2[2])));
    }
    return false;
    //double difference = distance_point_origin + distance_point_end - line_length;

    //if(abs(difference)<tol)
    //{
    //  on_line = true;
    //}

    //return on_line;
    // oppure potremmo togliere questo if e mettere solamente
    //return abs((distance_point_origin +distance_point_end)-line_length)<tol;
}

// QUELLA CHE AVEVAMO PRIMA
// bool point_on_line(Vector3d& line_origin, Vector3d& line_end, Vector3d& point)
// { // se il punto è più distante dall'inizio del segmento rispetto alla fine del segmento stesso allora il punto non appartiene al segmento
//     bool on_line = false;

//     double tol = numeric_limits<double>::epsilon();
//     double line_length = (line_origin - line_end).norm();
//     double distance_point_origin = (line_origin - point).norm();
//     double distance_point_end = (line_end - point).norm();
//     double difference = distance_point_origin + distance_point_end - line_length;

//     if(abs(difference)<tol)
//     {
//         on_line = true;
//     }

//     return on_line;
// }


bool check_barycentre_coord(const Vector3d& p0, const Vector3d& p1, const Vector3d& p2, const Vector3d& p3)
{
    Vector3d vector0 = p1-p3;
    Vector3d vector1 = p1-p2;
    Vector3d vector2 = p1-p0;

    double dotProductMatrix[2][2]={
        {vector0.dot(vector0), vector0.dot(vector1)},
        {vector0.dot(vector1), vector1.dot(vector1)}
    };

    // double dot00 = vector0.dot(vector0);
    // double dot01 = vector0.dot(vector1);
    // double dot02 = vector0.dot(vector2);
    // double dot11 = vector1.dot(vector1);
    // double dot12 = vector1.dot(vector2);
    double det = dotProductMatrix[0][0]*dotProductMatrix[1][1]-dotProductMatrix[0][1]*dotProductMatrix[1][0];
    if (det==0){
        return false;
    }

    double barycentre_coord1 = (1 / det) * (dotProductMatrix[1][1]*vector0.dot(vector2)-dotProductMatrix[0][1]*vector1.dot(vector2));
    double barycentre_coord2 = (1 / det) * (dotProductMatrix[0][0]*vector1.dot(vector2)-dotProductMatrix[0][1]*vector0.dot(vector2));

    return((barycentre_coord1 >= 0) && (barycentre_coord2 >= 0) && (barycentre_coord1+barycentre_coord2 <= 1));  // questi sono gli alfa0 e alfa1 di geom comp
}

bool check_inside_fracture(const Vector3d& point, vector<Vector3d>& fracture_vertex)
{
    for(unsigned int i = 0; i < fracture_vertex.size()-1; i++)
    {
        if(check_barycentre_coord(point, fracture_vertex[0], fracture_vertex[i], fracture_vertex[i+1]))
        {
            return true;
        }
    }
    return false;
}

// QUELLA CHE AVEVAMO PRIMA
// bool check_barycentre_coord(const Vector3d& p0, const Vector3d& p1, const Vector3d& p2, const Vector3d& p3)
// {
//     Vector3d vector0 = p1-p3;
//     Vector3d vector1 = p1-p2;
//     Vector3d vector2 = p1-p0;

//     double dot00 = vector0.dot(vector0);
//     double dot01 = vector0.dot(vector1);
//     double dot02 = vector0.dot(vector2);
//     double dot11 = vector1.dot(vector1);
//     double dot12 = vector1.dot(vector2);

//     double barycentre_coord1 = (1 / (dot00*dot11 - dot01*dot01)) * (dot11*dot02 - dot01*dot12);
//     double barycentre_coord2 = (1 / (dot00*dot11 - dot01*dot01)) * (dot00*dot12 - dot01*dot02);

//     if ((barycentre_coord1 >= 0) && (barycentre_coord2 >= 0) && (barycentre_coord1+barycentre_coord2 <= 1))  // questi sono gli alfa0 e alfa1 di geom comp
//     {
//         return true;
//     }
//     return false;
// }

// bool check_inside_fracture(const Vector3d& point, vector<Vector3d>& fracture_vertex)
// {
//     for(unsigned int i = 0; i < fracture_vertex.size()-1; i++)
//     {
//         if(check_barycentre_coord(point, fracture_vertex[0], fracture_vertex[i], fracture_vertex[i+1]))
//         {
//             return true;
//         }
//     }
//     return false;
// }

//riguardare questa funzione che è sbagliata
// cicla su ogni traccia e vedere da id traccia i due valori booleani delle fratture; se le due fratture generatrici passanti o no, peggio per riodinare

// void TracesType(GeometryDFN& dfn)
// {
//     //dfn.Traces_Id.reserve(dfn.Number_Traces);
//     for(unsigned int i = 0; i < dfn.Number_Traces; i++)
//     {
//         unsigned int number=0;
//         for(unsigned int j = 0; j < dfn.Number_Traces; j++)
//         {
//             unsigned int fracture_id = i;
//             if(fracture_id==dfn.Traces_Generator_Id[j][0] || fracture_id==dfn.Traces_Generator_Id[j][1])
//             {
//                 Vector2i fracture_type;
//                 for(unsigned int k = 0; k<2; k++)  // itero sui due punti estremi della traccia in questione
//                 {
//                     Vector3d point = dfn.Traces_Coordinates[j][k];
//                     bool result = true;
//                     for(unsigned int l = 0; l<dfn.Fractures_Number_Vertices[i]; l++)
//                     {
//                         Vector3d point1 = dfn.Fractures_Vertices[i][l];
//                         Vector3d point2 = dfn.Fractures_Vertices[i][(l+1)%dfn.Fractures_Number_Vertices[i]];
//                         if(point_on_line(point1,point2,point))
//                         {
//                             result = false;
//                         }
//                     }
//                     fracture_type[k] = result;
//                 }
//                 if(fracture_type[0]==false && fracture_type[1]==false)
//                 {
//                     dfn;
//                 }

//             }
//         }
//     }
// }


// void calcolaTipologiaTracce(GeometryDFN& DFN) {
//     for (auto it = DFN.Traces_Tips.begin(); it != DFN.Traces_Tips.end(); ++it) {
//         unsigned int id_traccia = it->first;
//         array<bool, 2> tipologiaTracce = it->second;

//         for (unsigned int j = 0; j < 2; j++) { // ciclo sui due estremi della traccia
//             Vector3d p_traccia = DFN.Traces_Coordinates[it->first][j];



//                 if (DFN.Traces_Id == DFN.Traces_Generator_Id[it->first][0] || DFN.Traces_Id == DFN.Traces_Generator_Id[it->first][1]) { // Verifica frattura
//                     for (unsigned int l = 0; l < DFN.Fractures_Vertices[i].size(); l++) {
//                         Vector3d p1 = DFN.Fractures_Vertices[i][l];
//                         Vector3d p2;
//                         if (l + 1 < DFN.Fractures_Vertices[i].size()) {
//                             p2 = DFN.Fractures_Vertices[i][l + 1];
//                         } else {
//                             p2 = DFN.Fractures_Vertices[i][0];
//                         }
//                         if (point_on_line(p1, p2, p_traccia)) {
//                             tipologiaTracce[j] = false; // È passante
//                             break;
//                         }
//                         else{
//                             tipologiaTracce[j] = true;
//                         }
//                     }
//                 }
//                 if (!tipologiaTracce[j]) break; // Se già determinato come passante, esci dal ciclo
//             }
//         }

//         // Assegna i risultati a Traces_Tips
//         DFN.Traces_Tips[it->first] = tipologiaTracce;
//     }
// }


/*void calcolaTipologiaTracce(GeometryDFN& DFN) {
    for (unsigned int k = 0; k < DFN.Number_Traces; k++) { // ciclo sulle tracce
        array<bool, 2> tipologiaTracce = {true, true}; // Assume non passante all'inizio

        for (unsigned int j = 0; j < 2; j++) { // ciclo sui due estremi della traccia
            Vector3d p_traccia = DFN.Traces_Coordinates[k][j];

            for (unsigned int i = 0; i < DFN.Number_Fractures; i++) { // ciclo sulle fratture
                int idFrattura = i;
                if (idFrattura == DFN.Traces_Generator_Id[k][0] || idFrattura == DFN.Traces_Generator_Id[k][1]) { // Verifica frattura
                    for (unsigned int l = 0; l < DFN.Fractures_Vertices[i].size(); l++) {
                        Vector3d p1 = DFN.Fractures_Vertices[i][l];
                        Vector3d p2;
                        if (l + 1 < DFN.Fractures_Vertices[i].size()) {
                            p2 = DFN.Fractures_Vertices[i][l + 1];
                        } else {
                            p2 = DFN.Fractures_Vertices[i][0];
                        }
                        if (point_on_line(p1, p2, p_traccia)) {
                            tipologiaTracce[j] = false; // È passante
                            break;
                        }
                        else{
                            tipologiaTracce[j] = true;
                        }
                    }
                }
                if (!tipologiaTracce[j]) break; // Se già determinato come passante, esci dal ciclo
            }
        }

        // Assegna i risultati a Traces_Tips
        DFN.Traces_Tips[k] = tipologiaTracce;
    }
}*/

// PROVA GIUSTAAAAAAAAAAAAAAAAAAAAA
void calcolaTipologiaTracce(GeometryDFN& DFN) {
    for (const auto& fracture_entry : DFN.Fractures_Vertices) { // ciclo sulle fratture
        unsigned int i = fracture_entry.first;
        const vector<Vector3d>& vertices = fracture_entry.second;
        for (unsigned int k = 0; k < DFN.Number_Traces; k++) { // ciclo sulle tracce
            int idFrattura = i;
            if (idFrattura == DFN.Traces_Generator_Id[k][0] || idFrattura == DFN.Traces_Generator_Id[k][1]) {
                array<bool, 2> tipologia = {true, true};
                for (unsigned int j = 0; j < 2; j++) {
                    Vector3d p_traccia = DFN.Traces_Coordinates[k][j];
                    for (unsigned int l = 0; l < vertices.size(); l++) {
                        Vector3d p1 = vertices[l];
                        Vector3d p2 = (l + 1 < vertices.size()) ? vertices[l + 1] : vertices[0];
                        if (point_on_line(p1, p2, p_traccia)) {
                            tipologia[j] = false; // La traccia è passante
                            break;
                        }
                    }
                }
                DFN.Traces_Tips[k][0] = !tipologia[0] && !tipologia[1]; // La traccia è passante se entrambe sono false
                DFN.Traces_Tips[k][1] = tipologia[0] || tipologia[1];  // La traccia non è passante se almeno una è true

            }
        }
    }
}



// DECENTE
// void calcolaTipologiaTracce(GeometryDFN& DFN) {
//     //DFN.Traces_Tips.resize(DFN.Number_Fractures); // Assicura che tipoTraccia abbia la dimensione corretta
//     for (unsigned int k = 0; k < DFN.Number_Traces; k++) { // ciclo sulle tracce
//         //array<bool, 2> presenzaTracce = {};
//         unsigned int tracceCount = 0; // Contatore per tracce per la frattura k
//         for (unsigned int i = 0; i < DFN.Number_Fractures; i++) { // ciclo sulle fratture
//             int idFrattura = i;
//             if (idFrattura == DFN.Traces_Generator_Id[k][0] || idFrattura == DFN.Traces_Generator_Id[k][1]) {
//                 array<bool, 2> tipologia = {};
//                 for (unsigned int j = 0; j < 2; j++) {
//                     Vector3d p_traccia = DFN.Traces_Coordinates[k][j];
//                     for (unsigned int l = 0; l < DFN.Fractures_Vertices[i].size(); l++) {
//                         Vector3d p1 = DFN.Fractures_Vertices[i][l];
//                         Vector3d p2;
//                         if (l + 1 < DFN.Fractures_Vertices[i].size()) {
//                             p2 = DFN.Fractures_Vertices[i][l + 1];
//                         } else {
//                             p2 = DFN.Fractures_Vertices[i][0];
//                         }
//                         if (point_on_line(p1, p2, p_traccia)) {
//                             tipologia[j] = false;
//                             break;
//                         }
//                     }
//                 }
//                 if (tipologia[0] == false && tipologia[1] == false) {
//                     DFN.Traces_Tips[k][tracceCount] = {false};
//                     //presenzaTracce[0] = true;
//                 } else {
//                     DFN.Traces_Tips[k][tracceCount] = {true};
//                     //presenzaTracce[1] = true;
//                 }
//                 tracceCount++;
//             }
//         }

//     }
// }

// QUELLA CHE AVEVAMO PRIMA
// void calcolaTipologiaTracce(GeometryDFN& DFN) {
//     //DFN.Traces_Tips.resize(DFN.Number_Fractures); // Assicura che tipoTraccia abbia la dimensione corretta
//     for (unsigned int k = 0; k < DFN.Number_Fractures; k++) { // ciclo sulle Fratture
//         array<bool, 2> presenzaTracce = {false, false};
//         unsigned int tracceCount = 0; // Contatore per tracce per la frattura k
//         for (unsigned int i = 0; i < DFN.Number_Traces; i++) { // ciclo sulle tracce
//             int idFrattura = k;
//             if (idFrattura == DFN.Traces_Generator_Id[i][0] || idFrattura == DFN.Traces_Generator_Id[i][1]) {
//                 array<bool, 2> tipologia = {true, true};
//                 for (unsigned int j = 0; j < 2; j++) {
//                     Vector3d p_traccia = DFN.Traces_Coordinates[i][j];
//                     for (unsigned int l = 0; l < DFN.Fractures_Vertices[k].size(); l++) {
//                         Vector3d p1 = DFN.Fractures_Vertices[k][l];
//                         Vector3d p2;
//                         if (l + 1 < DFN.Fractures_Vertices[k].size()) {
//                             p2 = DFN.Fractures_Vertices[k][l + 1];
//                         } else {
//                             p2 = DFN.Fractures_Vertices[k][0];
//                         }
//                         if (point_on_line(p1, p2, p_traccia)) {
//                             tipologia[j] = false;
//                             break;
//                         }
//                     }
//                 }
//                 if (tipologia[0] == false && tipologia[1] == false) {
//                     DFN.Traces_Tips[k][tracceCount] = {false};
//                     presenzaTracce[0] = true;
//                 } else {
//                     DFN.Traces_Tips[k][tracceCount] = {true};
//                     presenzaTracce[1] = true;
//                 }
//                 tracceCount++;
//             }
//         }
//         if (presenzaTracce[0] || presenzaTracce[1]) {
//             DFN.Traces_Tips[k] = presenzaTracce;
//         }
//     }
// }



// void calcolaTipologiaTracce(GeometryDFN& DFN) {
//     for (unsigned int k = 0; k < DFN.Number_Fractures; k++) { // ciclo sulle Fratture
//         array<bool, 2> presenzaTracce = {false, false};
//         for (unsigned int i = 0; i < DFN.Number_Traces; i++) { // ciclo sulle tracce
//             if (isFratturaIntersecante(k, i, DFN)) {
//                 array<bool, 2> tipologia = getTipologia(i, k, DFN);
//                 if (tipologia[0] == false && tipologia[1] == false) {
//                     DFN.Traces_Tips[k].push_back({i, false});
//                     presenzaTracce[0] = true;
//                 } else {
//                     DFN.Traces_Tips[k].emplace_back(i, true);
//                     presenzaTracce[1] = true;
//                 }
//             }
//         }
//         if (presenzaTracce[0] || presenzaTracce[1]) {
//             DFN.Id_FrattureConTraccia[k] = presenzaTracce;
//         }
//     }
// }

// bool isFratturaIntersecante(unsigned int k, unsigned int i, GeometryDFN& DFN) {
//     int idFrattura = k;
//     return idFrattura == DFN.Traces_Generator_Id[i][0] || idFrattura == DFN.Traces_Generator_Id[i][1];
// }

// array<bool, 2> getTipologia(unsigned int i, unsigned int k, GeometryDFN& DFN) {
//     array<bool, 2> tipologia;
//     for (unsigned int j = 0; j < 2; j++) {
//         Vector3d p_traccia = DFN.Traces_Coordinates[i][j];
//         bool risposta = true;
//         // prendere due vertici consecutivi
//         for (unsigned int l = 0; l < DFN.Fractures_Number_Vertices[k]; l++) {
//             Vector3d p1 = DFN.Fractures_Vertices[k][l];
//             Vector3d p2;
//             if (l + 1 < DFN.Fractures_Number_Vertices[k]) {
//                 p2 = DFN.Fractures_Vertices[k][l + 1];
//             } else {
//                 p2 = DFN.Fractures_Vertices[k][0];
//             }
//             if (point_on_line(p1, p2, p_traccia)) {
//                 risposta = false;
//                 break;
//             }
//         }
//         tipologia[j] = risposta;
//     }
//     return tipologia;
// }



//vector<Vector2i> sorting(const vector<double>& length, const vector<Vector2i>& type)
map<unsigned int,array<bool, 2>> riordinaTracce(const vector<double>& length, map<unsigned int,array<bool, 2>>& type, vector<unsigned int>& trace_id)
{
    vector<pair<int,double>> pairLengthID;
    for (unsigned int i=0; i<length.size();i++)
    {
        pairLengthID.push_back({i,length[i]});
    }
    SortLibrary::MergeSort(pairLengthID);

    vector<array<bool, 2>> passante;
    vector<array<bool, 2>> non_passante;
    for (unsigned int i=0; i<pairLengthID.size();i++)
    {
       unsigned int index = pairLengthID[i].first;  // indice della traccia
       unsigned int t_id = trace_id[index];


       if (type.at(t_id)[1]==false)
       {
           passante.push_back(type.at(t_id));
       }
        else if (type.at(t_id)[1]==true)
        {
           non_passante.push_back(type.at(t_id));
        }

    }
    map<unsigned int,array<bool,2>> sorted;
    unsigned int passanteIdx = 0, nonPassanteIdx = 0;

    // Riempie il map con le tracce passanti
    for (unsigned int i = 0; i < passante.size(); i++)
    {
        unsigned int index = pairLengthID[i].first; // indice della traccia
        unsigned int t_id = trace_id[index];
        sorted[t_id] = passante[passanteIdx++];
    }

    // Riempie il map con le tracce non passanti
    for (unsigned int i = 0; i < non_passante.size(); i++)
    {
        unsigned int index = pairLengthID[passante.size() + i].first; // indice della traccia
        unsigned int t_id = trace_id[index];
        sorted[t_id] = non_passante[nonPassanteIdx++];
    }

    return sorted;
    // sorted.insert(sorted.end(),passante.begin(),passante.end());
    // sorted.insert(sorted.end(),non_passante.begin(),non_passante.end());
}

void calcolaLunghezzaTracce(GeometryDFN& DFN)
{
    // DFN.traces_length.resize(DFN.Number_Traces);
    // unsigned int index = 0;

    for (auto& entry : DFN.Traces_Coordinates) {
        double lunghezza = (entry.second[1] - entry.second[0]).norm();
        DFN.traces_length.push_back(lunghezza);
        DFN.Traces_Id.push_back(entry.first);
    }

    DFN.Traces_Tips = GeometryLibrary::riordinaTracce(DFN.traces_length, DFN.Traces_Tips, DFN.Traces_Id);

    //vector<unsigned int> trace_id;
    // for (const auto& item : DFN.Traces_Tips) {

    // }
    // for(auto& entry : DFN.Traces_Tips)
    // {
    //     entry.second = riordinaTracce(DFN.traces_length, DFN.Traces_Tips, DFN.Traces_Id);
    // }
    //entry.second = riordinaTracce(DFN.traces_length, DFN.Traces_Tips, trace_id);

}

bool OutputTracce(const GeometryDFN& DFN, const string& fileOutput)
{
    ofstream file;
    file.open(fileOutput);
    if (file.fail())
    {
        cerr << "errore nell'aprire il file di output: " << fileOutput << endl;
        return false;
    }
    file << "# Number of Traces" << endl;
    file.precision(16);
    file << scientific;
    file << DFN.Number_Traces << endl;
    file << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;

    for (unsigned int i = 0; i < DFN.Number_Traces; i++)
    {
        unsigned int traceId = DFN.Traces_Id[i];
        Vector2i fractureIds = DFN.Traces_Generator_Id[i];
        array<Vector3d, 2> coordinates = DFN.Traces_Coordinates.at(traceId);

        file << traceId << "; " << fractureIds[0] << "; " << fractureIds[1] << "; "
             << coordinates[0][0] << "; " << coordinates[0][1] << "; " << coordinates[0][2] << "; "
             << coordinates[1][0] << "; " << coordinates[1][1] << "; " << coordinates[1][2] << endl;
    }
    file.close();
    return true;
}


bool OutputFratture(const GeometryDFN& DFN, const string& fileOutput)
{
    ofstream file;
    file.open(fileOutput);
    if (file.fail())
    {
        cerr << "errore nell'aprire il file di output: " << fileOutput << endl;
        return false;
    }
    file.precision(16);
    file << scientific;

    for (unsigned int i = 0; i < DFN.Number_Fractures; i++)
    {
        unsigned int fractureId = DFN.Fractures_Id[i];
        vector<unsigned int> traceIds;

        for (unsigned int j = 0; j < DFN.Number_Traces; j++)
        {
            Vector2i fractureIds = DFN.Traces_Generator_Id[j];
            if (fractureIds[0] == fractureId || fractureIds[1] == fractureId)
            {
                traceIds.push_back(DFN.Traces_Id[j]);
            }
        }

        if (!traceIds.empty())
        {
            file << "# FractureId; NumTraces" << endl;
            file << fractureId << "; " << traceIds.size() << endl;
            file << "# TraceId; Tips; Length" << endl;

            for (const auto& traceId : traceIds)
            {
                array<bool, 2> tips = DFN.Traces_Tips.at(traceId);
                double length = DFN.traces_length[traceId];



                // Determinare se è passante o non passante
                //int tipType = (tips[0] && tips[1]) ? 1 : 0;
                int tipType;
                if (tips[0] && tips[1]) {
                    tipType = 1;
                } else {
                    tipType = 0;
                }

                file << traceId << "; " << tipType << "; " << length << endl;
            }
        }
    }
    file.close();
    return true;
}


// bool OutputFratture(const GeometryDFN& DFN, const string& fileOutput)
// {
//     ofstream file;
//     file.open(fileOutput);
//     if (file.fail())
//     {
//         cerr << "errore nell'aprire il file di output: " << fileOutput << endl;
//         return false;
//     }
//     file.precision(16);
//     file << scientific;

//     for (unsigned int i = 0; i < DFN.Number_Fractures; i++)
//     {
//         unsigned int fractureId = DFN.Fractures_Id[i];
//         vector<unsigned int> traceIds;

//         for (unsigned int j = 0; j < DFN.Number_Traces; j++)
//         {
//             Vector2i fractureIds = DFN.Traces_Generator_Id[j];
//             if (fractureIds[0] == fractureId || fractureIds[1] == fractureId)
//             {
//                 traceIds.push_back(DFN.Traces_Id[j]);
//             }
//         }

//         if (!traceIds.empty())
//         {
//             file << "# FractureId; NumTraces" << endl;
//             file << fractureId << "; " << traceIds.size() << endl;
//             file << "# TraceId; Tips; Length" << endl;

//             for (const auto& traceId : traceIds)
//             {
//                 array<bool, 2> tips = DFN.Traces_Tips.at(traceId);
//                 double length = DFN.traces_length[traceId];
//                 file << traceId << "; " << tips[0] << ", " << tips[1] << "; " << length << endl;
//             }
//         }
//     }
//     file.close();
//     return true;
// }



}


// metti in ordine prima in base alla lunghezza e poi guardi i booleani

