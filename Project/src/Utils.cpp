#include "Utils.hpp"
#include <iostream>
#include <array>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "Eigen/Eigen"
#include "Eigen/Dense"

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

    if ((radius1 + radius2 + tol) > barycentres_distance)
    {
        return false;
    }
    return true;
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


// COSI FUNZIONAVA TUTTO
// void FindTraces(GeometryDFN& dfn)
// {
//     unsigned int index_trace = 0;
//     double tol = numeric_limits<double>::epsilon();

//     // CONTROLLARE COSA CAMBIA SENZA IL -1; il -1 serve per l'id
//     for(unsigned int i = 0; i < dfn.Number_Fractures - 1; i++){
//         for(unsigned int j = i+1; j<dfn.Number_Fractures; j++){

//             if(parallel_planes(dfn.Fractures_Vertices[i], dfn.Fractures_Vertices[j])){
//                 cout << "Le fratture: " << dfn.Fractures_Id[i] << " e " << dfn.Fractures_Id[j] << " sono parallele" << endl;
//                 continue;           // CONTROLLA LA FUNZIONALITA DI CONTINUE
//             }

//             if(!FirstSelectionTraces(dfn.Fractures_Vertices[i], dfn.Fractures_Vertices[j], tol)){   // CONTROLLA IL !
//                 Vector3d Normal_i = NormalToPlane(dfn.Fractures_Vertices[i]);
//                 Vector3d Normal_j = NormalToPlane(dfn.Fractures_Vertices[j]);
//                 Vector3d T = Normal_i.cross(Normal_j);
//                 Matrix3d Planes_Matrix;
//                 Planes_Matrix.row(0) = Normal_i;
//                 Planes_Matrix.row(1) = Normal_j;
//                 Planes_Matrix.row(2) = T;
//                 double d1 = Normal_i.dot(FindBarycentre(dfn.Fractures_Vertices[i]));
//                 double d2 = Normal_j.dot(FindBarycentre(dfn.Fractures_Vertices[j]));
//                 double d3 = 0;
//                 Vector3d b = {d1,d2,d3};

//                 double det_matrix = Planes_Matrix.determinant();
//                 if(det_matrix > tol){
//                     Vector3d p_intersection = Planes_Matrix.fullPivLu().solve(b); // trovato vettore p1,p2,p3 del foglio geom comp 2
//                     unsigned int num_traces_points = 0;
//                     array<Vector3d,2> traces_points;
//                     for(unsigned int k = 0; k < dfn.Fractures_Number_Vertices[i]; k++)
//                     {                  // SE NON FUNZIONA SICURAMENTE QUA QUALCOSA NON VA BENE; CAPIRE SE QUESTO FOR DEVE GIRARE SU K O SU J

//                         Vector3d v1 = dfn.Fractures_Vertices[i][k];
//                         Vector3d v2 = dfn.Fractures_Vertices[i][(k + 1) % dfn.Fractures_Number_Vertices[i]]; // faccio sta roba per non andare fuori scope

//                         MatrixXd A = MatrixXd::Zero(3,2);

//                         A.col(0) = v2-v1;   // è una retta così come t
//                         A.col(1) = T;  // per tornare al valore positivo di beta dal -beta che abbiamo nella sottrazione tra equazioni delle due rette

//                         if((v2-v1).cross(T).squaredNorm() < tol)
//                         {
//                             continue;
//                         }


//                         double b0 = p_intersection[0] - v1[0];
//                         double b1 = p_intersection[1] - v1[1];
//                         double b2 = p_intersection[2] - v1[2];
//                         Vector3d b = {b0,b1,b2};

//                         Vector2d alpha_beta = A.householderQr().solve(b);

//                         if(alpha_beta[0]>-tol && alpha_beta[0]<1+tol && num_traces_points < 2)
//                         {//controllo che non prenda due volte lo stesso punto
//                             Vector3d Intersection = v1+(alpha_beta[0]*(v2-v1));
//                             if (num_traces_points == 0 && check_inside_fracture(Intersection, dfn.Fractures_Vertices[j]))
//                                 {
//                                     traces_points[num_traces_points] = Intersection;
//                                     num_traces_points += 1;
//                                 }
//                             else if (num_traces_points == 1 && check_inside_fracture(Intersection, dfn.Fractures_Vertices[j]) && traces_points[0] !=Intersection)
//                             { traces_points[num_traces_points]=Intersection;
//                                 num_traces_points+=1;
//                                     dfn.Number_Traces+=1;
//                                 dfn.Traces_Id.push_back(index_trace);
//                                 Vector2i fracture = {dfn.Fractures_Id[i], dfn.Fractures_Id[j]};
//                                 dfn.Traces_Generator_Id.push_back(fracture);
//                                 dfn.Traces_Coordinates[index_trace] = traces_points;
//                                 index_trace += 1;
//                             }

//                         }

//                     } // qui controllo sulla seconda frattura, almeno da trovare anche il secondo punto di intersezione nel caso rosso-blu in FR3
//                     for(unsigned int k = 0; k < dfn.Fractures_Number_Vertices[j]; k++)
//                     {                  // SE NON FUNZIONA SICURAMENTE QUA QUALCOSA NON VA BENE; CAPIRE SE QUESTO FOR DEVE GIRARE SU K O SU J

//                         Vector3d v1 = dfn.Fractures_Vertices[j][k];
//                         Vector3d v2 = dfn.Fractures_Vertices[j][(k + 1) % dfn.Fractures_Number_Vertices[i]];

//                         MatrixXd A = MatrixXd::Zero(3,2);

//                         A.col(0) = v2-v1;
//                         A.col(1) = T;  // per tornare al valore positivo di beta dal -beta che abbiamo nella sottrazione tra equazioni delle due rette

//                         if((v2-v1).cross(T).squaredNorm() < tol)
//                         {
//                             continue;
//                         }


//                         double b0 = p_intersection[0] - v1[0];
//                         double b1 = p_intersection[1] - v1[1];
//                         double b2 = p_intersection[2] - v1[2];
//                         Vector3d b = {b0,b1,b2};

//                         Vector2d alpha_beta = A.householderQr().solve(b);

//                         if(alpha_beta[0]>-tol && alpha_beta[0]<1+tol && num_traces_points < 2)
//                         {
//                             Vector3d Intersection = v1+(alpha_beta[0]*(v2-v1));
//                             if (num_traces_points == 0 && check_inside_fracture(Intersection, dfn.Fractures_Vertices[i]))
//                             {
//                                 traces_points[num_traces_points] = Intersection;
//                                 num_traces_points += 1;
//                             }
//                             else if (num_traces_points == 1 && check_inside_fracture(Intersection, dfn.Fractures_Vertices[i]) && traces_points[0] !=Intersection)
//                             { traces_points[num_traces_points]=Intersection;
//                                 num_traces_points+=1;
//                                 dfn.Number_Traces+=1;
//                                 dfn.Traces_Id.push_back(index_trace);
//                                 Vector2i fracture = {dfn.Fractures_Id[i], dfn.Fractures_Id[j]};
//                                 dfn.Traces_Generator_Id.push_back(fracture);
//                                 dfn.Traces_Coordinates[index_trace] = traces_points;
//                                 index_trace += 1;
//                             }

//                         }

//                     }
//                 }
//             }
//         }
//     }
// }


void FindTraces(GeometryDFN& dfn)
{
    unsigned int index_trace = 0;
    double tol = numeric_limits<double>::epsilon();

    for(unsigned int i = 0; i < dfn.Number_Fractures - 1; ++i){
        for(unsigned int j = i + 1; j < dfn.Number_Fractures; ++j){

            if(parallel_planes(dfn.Fractures_Vertices[i], dfn.Fractures_Vertices[j])){
                cout << "Fractures " << dfn.Fractures_Id[i] << " and " << dfn.Fractures_Id[j] << " are parallel" << endl;
                continue;
            }

            if(!FirstSelectionTraces(dfn.Fractures_Vertices[i], dfn.Fractures_Vertices[j], tol)){
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
                Vector3d b(d1, d2, d3);

                double det_matrix = Planes_Matrix.determinant();
                if(det_matrix > tol){
                    Vector3d p_intersection = Planes_Matrix.fullPivLu().solve(b);
                    array<Vector3d, 2> traces_points;
                    unsigned int num_traces_points = 0;

                    auto find_traces = [&](const vector<Vector3d>& vertices, unsigned int fracture_id, unsigned int other_fracture_id) {
                        for(unsigned int k = 0; k < vertices.size(); ++k) {
                            Vector3d v1 = vertices[k];
                            Vector3d v2 = vertices[(k + 1) % vertices.size()];

                            MatrixXd A(3, 2);
                            A.col(0) = v2 - v1;
                            A.col(1) = T;

                            if((v2 - v1).cross(T).squaredNorm() < tol) {
                                continue;
                            }

                            Vector3d b(p_intersection[0] - v1[0], p_intersection[1] - v1[1], p_intersection[2] - v1[2]);
                            Vector2d alpha_beta = A.householderQr().solve(b);

                            if(alpha_beta[0] > -tol && alpha_beta[0] < 1 + tol && num_traces_points < 2) {
                                Vector3d Intersection = v1 + (alpha_beta[0] * (v2 - v1));
                                if(num_traces_points == 0 && check_inside_fracture(Intersection, dfn.Fractures_Vertices[other_fracture_id])) {
                                    traces_points[num_traces_points++] = Intersection;
                                } else if(num_traces_points == 1 && check_inside_fracture(Intersection, dfn.Fractures_Vertices[other_fracture_id]) && traces_points[0] != Intersection) {
                                    traces_points[num_traces_points++] = Intersection;
                                    dfn.Number_Traces++;
                                    dfn.Traces_Id.push_back(index_trace);
                                    dfn.Traces_Generator_Id.push_back({dfn.Fractures_Id[i], dfn.Fractures_Id[j]});
                                    dfn.Traces_Coordinates[index_trace++] = traces_points;
                                }
                            }
                        }
                    };

                    find_traces(dfn.Fractures_Vertices[i], i, j);
                    find_traces(dfn.Fractures_Vertices[j], j, i);
                }
            }
        }
    }
}

bool point_on_line(Vector3d& line_origin, Vector3d& line_end, Vector3d& point)
{
    bool on_line = false;

    double tol = 1e-12;
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


bool check_barycentre_coord(const Vector3d& p0, const Vector3d& p1, const Vector3d& p2, const Vector3d& p3)
{
    Vector3d vector0 = p1-p3;
    Vector3d vector1 = p1-p2;
    Vector3d vector2 = p1-p0;

    double dotProductMatrix[2][2]={
        {vector0.dot(vector0), vector0.dot(vector1)},
        {vector0.dot(vector1), vector1.dot(vector1)}
    };

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


// FUNZIONA QUASI PERFETTAMENTE
void FindTracesType(GeometryDFN& DFN) {
    for (const auto& fracture_entry : DFN.Fractures_Vertices) { // ciclo sulle fratture
        unsigned int i = fracture_entry.first;
        const vector<Vector3d>& vertices = fracture_entry.second;
        for (unsigned int k = 0; k < DFN.Number_Traces; k++) { // ciclo sulle tracce
            int idFrattura = i;
            if (idFrattura == DFN.Traces_Generator_Id[k][0]) {
                array<bool, 2> type = {true, true};
                //unsigned int contatore1 = 0;
                for (unsigned int j = 0; j < 2; j++) {
                    Vector3d p_traccia = DFN.Traces_Coordinates[k][j];
                    for (unsigned int l = 0; l < vertices.size(); l++) {
                        Vector3d p1 = vertices[l];
                        Vector3d p2 = (l + 1 < vertices.size()) ? vertices[l + 1] : vertices[0];
                        if (point_on_line(p1, p2, p_traccia)) {
                            type[j] = false; // La traccia è passante
                            //contatore1 += 1;
                            break;
                        }
                    }
                }
                if (type[0]==false && type[1]==false){
                    DFN.Traces_Tips[k][0] = false; // La traccia è passante se entrambe sono false
                }
                else{
                    DFN.Traces_Tips[k][0] = true;
                }
            }
            if ( idFrattura == DFN.Traces_Generator_Id[k][1]) {
                array<bool, 2> type = {true, true};
                for (unsigned int j = 0; j < 2; j++) {
                    Vector3d p_traccia = DFN.Traces_Coordinates[k][j];
                    for (unsigned int l = 0; l < vertices[idFrattura].size(); l++) {
                        Vector3d p1 = vertices[l];
                        Vector3d p2 = (l + 1 < vertices.size()) ? vertices[l + 1] : vertices[0];
                        if (point_on_line(p1, p2, p_traccia)) {
                            type[j] = false; // La traccia è passante
                            break;
                        }
                    }
                }
                if (type[0]==false && type[1]==false){
                    DFN.Traces_Tips[k][1] = false; // La traccia è passante se entrambe sono false
                }
                else{
                    DFN.Traces_Tips[k][1] = true;
                }
            }

        }
    }
}


vector<unsigned int> reorganiseLength(const vector<unsigned int>& traceIds, const GeometryDFN& DFN)
{
    vector<pair<double, unsigned int>> pairLengthID;
    for (const auto& traceId : traceIds)
    {
        auto it = find(DFN.Traces_Id.begin(), DFN.Traces_Id.end(), traceId);
        if (it != DFN.Traces_Id.end())
        {
            unsigned int index = distance(DFN.Traces_Id.begin(), it);
            double length = DFN.traces_length[index];
            pairLengthID.push_back({length, traceId});
        }
    }

    sort(pairLengthID.begin(), pairLengthID.end(), [](const auto& left, const auto& right) {
        return left.first > right.first;
    });

    vector<unsigned int> sortedTraceIDs;
    for (const auto& pair : pairLengthID)
    {
        sortedTraceIDs.push_back(pair.second);
    }

    return sortedTraceIDs;
}


void TracesLength(GeometryDFN& DFN)
{

    for (auto& entry : DFN.Traces_Coordinates) {
        double length = (entry.second[1] - entry.second[0]).norm();
        DFN.traces_length.push_back(length);
        DFN.Traces_Id.push_back(entry.first);
    }
}


bool OutputTraces(const GeometryDFN& DFN, const string& fileOutput)
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


bool OutputFractures(const GeometryDFN& DFN, const string& fileOutput)
{
    ofstream file;
    file.open(fileOutput);
    if (file.fail())
    {
        cerr << "Errore nell'aprire il file di output: " << fileOutput << endl;
        return false;
    }
    file.precision(16);
    file << scientific;

    for (unsigned int i = 0; i < DFN.Number_Fractures; i++)
    {
        int fractureId = DFN.Fractures_Id[i];
        vector<unsigned int> passanti;
        vector<unsigned int> nonPassanti;

        for (unsigned int j = 0; j < DFN.Number_Traces; j++)
        {
            Vector2i fractureIds = DFN.Traces_Generator_Id[j];
            if (fractureIds[0] == fractureId || fractureIds[1] == fractureId)
            {
                unsigned int traceId = DFN.Traces_Id[j];
                array<bool, 2> tips = DFN.Traces_Tips.at(traceId);

                if (tips[0] && tips[1]) {
                    nonPassanti.push_back(traceId);
                } else {
                    passanti.push_back(traceId);
                }
            }
        }

        // Ordinare per lunghezza
        auto sortedPassanti = reorganiseLength(passanti, DFN);
        auto sortedNonPassanti = reorganiseLength(nonPassanti, DFN);


        if (!sortedPassanti.empty() || !sortedNonPassanti.empty())
        {
            file << "# FractureId; NumTraces" << endl;
            file << fractureId << "; " << sortedPassanti.size() + sortedNonPassanti.size() << endl;
            file << "# TraceId; Tips; Length" << endl;



            for (const auto& traceId : sortedPassanti)
            {
                double length = DFN.traces_length[traceId];
                int tipType = 0;

                file << traceId << "; " << tipType << "; " << length << endl;
            }

            for (const auto& traceId : sortedNonPassanti)
            {
                double length = DFN.traces_length[traceId];
                int tipType = 1;

                file << traceId << "; " << tipType << "; " << length << endl;
            }
        }
    }
    file.close();
    return true;
}

}
