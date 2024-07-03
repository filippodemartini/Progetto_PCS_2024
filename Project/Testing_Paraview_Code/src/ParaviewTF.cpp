#include "ParaviewTF.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

namespace DFNLibrary{

bool stampaDatiSulFileFrattureParaview(const string &percorsoFileFrattureParaview, GeometryLibrary::GeometryDFN &Fratt){

    ofstream fileFrattureParaview;
    fileFrattureParaview.open(percorsoFileFrattureParaview);

    if(fileFrattureParaview.fail()){

        cerr << "errore: impossibile aprire il file " << percorsoFileFrattureParaview << endl;
        return false;

    }else{

        // Scrivo l'intestazione
        fileFrattureParaview << "# vtk DataFile Version 3.0" << endl;
        fileFrattureParaview << "VTK file for polygon visualization" << endl;
        fileFrattureParaview << "ASCII" << endl;
        fileFrattureParaview << "DATASET POLYDATA" << endl;

        unsigned int puntiTotali = 0;
        unsigned int poligoniTotali = 0;

        // Calcolo il numero totale di punti e poligoni
        for(unsigned int idFrattura = 0; idFrattura < Fratt.Fractures_Vertices.size(); idFrattura++){
            MatrixXd matrCoordinateFrattura = Fratt.Fractures_Vertices[idFrattura];
            puntiTotali += matrCoordinateFrattura.cols();
            poligoniTotali += 1;  // Supponendo che ogni frattura rappresenti un poligono
        }

        // Scrivo i punti
        fileFrattureParaview << "POINTS " << puntiTotali << " double" << endl;
        for(unsigned int idFrattura = 0; idFrattura < Fratt.Fractures_Vertices.size(); idFrattura++){
            MatrixXd matrCoordinateFrattura = Fratt.Fractures_Vertices[idFrattura];
            for(unsigned int i = 0; i < matrCoordinateFrattura.cols(); i++){
                fileFrattureParaview << fixed << scientific << setprecision(16) << matrCoordinateFrattura(0, i) << " "
                                     << matrCoordinateFrattura(1, i) << " " << matrCoordinateFrattura(2, i) << endl;
            }
        }

        // Scrivo i poligoni
        fileFrattureParaview << "POLYGONS " << poligoniTotali << " " << poligoniTotali * 5 << endl;
        unsigned int indicePunto = 0;
        for(unsigned int idFrattura = 0; idFrattura < Fratt.Fractures_Vertices.size(); idFrattura++){
            MatrixXd matrCoordinateFrattura = Fratt.Fractures_Vertices[idFrattura];
            fileFrattureParaview << "4 "
                                 << indicePunto << " "
                                 << indicePunto + 1 << " "
                                 << indicePunto + 2 << " "
                                 << indicePunto + 3 << endl;
            indicePunto += matrCoordinateFrattura.cols();
        }

        // Aggiungo dati scalari ai punti
        fileFrattureParaview << "POINT_DATA " << puntiTotali << endl;
        fileFrattureParaview << "SCALARS FractureId int 1" << endl;
        fileFrattureParaview << "LOOKUP_TABLE default" << endl;
        for(unsigned int idFrattura = 0; idFrattura < Fratt.Fractures_Vertices.size(); idFrattura++){
            MatrixXd matrCoordinateFrattura = Fratt.Fractures_Vertices[idFrattura];
            for(unsigned int i = 0; i < matrCoordinateFrattura.cols(); i++){
                fileFrattureParaview << idFrattura << endl;
            }
        }

        // Aggiungo dati scalari ai poligoni per il colore
        fileFrattureParaview << "CELL_DATA " << poligoniTotali << endl;
        fileFrattureParaview << "SCALARS FractureColor int 1" << endl;
        fileFrattureParaview << "LOOKUP_TABLE default" << endl;
        for(unsigned int idFrattura = 0; idFrattura < Fratt.Fractures_Vertices.size(); idFrattura++){
            fileFrattureParaview << idFrattura << endl;
        }

        fileFrattureParaview.close();
        return true;
    }
}

bool stampaDatiSulFileTracceParaview(const string &percorsoFileTracceParaview, GeometryLibrary::GeometryDFN &Trac){

    ofstream fileTracceParaview;
    fileTracceParaview.open(percorsoFileTracceParaview);

    if(fileTracceParaview.fail()){
        cerr << "errore: impossibile aprire il file " << percorsoFileTracceParaview << endl;
        return false;
    }else{

        // Scrivo l'intestazione
        fileTracceParaview << "# vtk DataFile Version 3.0" << endl;
        fileTracceParaview << "VTK file for point and line visualization" << endl;
        fileTracceParaview << "ASCII" << endl;
        fileTracceParaview << "DATASET POLYDATA" << endl;

        // Scrivo i punti
        fileTracceParaview << "POINTS " << 2 * Trac.coordinateIntersezioniTracce.size() << " double" << endl;
        for(const auto& intersezioni : Trac.coordinateIntersezioniTracce){
            const MatrixXd& coords = intersezioni.second;
            fileTracceParaview << fixed << scientific << setprecision(16) << coords(0, 0) << " " << coords(0, 1) << " " << coords(0, 2) << endl;
            fileTracceParaview << fixed << scientific << setprecision(16) << coords(1, 0) << " " << coords(1, 1) << " " << coords(1, 2) << endl;
        }

        // Scrivo i segmenti
        fileTracceParaview << "LINES " << Trac.coordinateIntersezioniTracce.size() << " " << 3 * Trac.coordinateIntersezioniTracce.size() << endl;
        for(unsigned int i = 0; i < Trac.coordinateIntersezioniTracce.size(); i++){
            fileTracceParaview << "2 " << 2 * i << " " << 2 * i + 1 << endl;
        }

        // Aggiungo dati scalari ai punti
        fileTracceParaview << "POINT_DATA " << 2 * Trac.coordinateIntersezioniTracce.size() << endl;

        fileTracceParaview << "SCALARS TraceId int 1" << endl;
        fileTracceParaview << "LOOKUP_TABLE default" << endl;
        for(unsigned int i = 0; i < Trac.coordinateIntersezioniTracce.size(); i++){
            fileTracceParaview << i << endl;
            fileTracceParaview << i << endl;
        }

        fileTracceParaview << "SCALARS FractureId1 int 1" << endl;
        fileTracceParaview << "LOOKUP_TABLE default" << endl;
        for(const auto& intersezioni : Trac.coordinateIntersezioniTracce){
            fileTracceParaview << intersezioni.first.first << endl;
            fileTracceParaview << intersezioni.first.first << endl;
        }

        fileTracceParaview << "SCALARS FractureId2 int 1" << endl;
        fileTracceParaview << "LOOKUP_TABLE default" << endl;
        for(const auto& intersezioni : Trac.coordinateIntersezioniTracce){
            fileTracceParaview << intersezioni.first.second << endl;
            fileTracceParaview << intersezioni.first.second << endl;
        }

        // Aggiungo vettori ai punti (direzioni dei segmenti)
        fileTracceParaview << "VECTORS SegmentDirection double" << endl;
        for(const auto& intersezioni : Trac.coordinateIntersezioniTracce){
            const MatrixXd& coords = intersezioni.second;
            double dx = coords(1, 0) - coords(0, 0);
            double dy = coords(1, 1) - coords(0, 1);
            double dz = coords(1, 2) - coords(0, 2);
            fileTracceParaview << fixed << scientific << setprecision(16) << dx << " " << dy << " " << dz << endl;
            fileTracceParaview << fixed << scientific << setprecision(16) << dx << " " << dy << " " << dz << endl;
        }

        // Aggiungo dati scalari ai segmenti
        fileTracceParaview << "CELL_DATA " << Trac.coordinateIntersezioniTracce.size() << endl;

        fileTracceParaview << "SCALARS SegmentLength double 1" << endl;
        fileTracceParaview << "LOOKUP_TABLE default" << endl;
        for(const auto& intersezioni : Trac.coordinateIntersezioniTracce){
            const MatrixXd& coords = intersezioni.second;
            Vector3d punto1(coords(0, 0), coords(0, 1), coords(0, 2));
            Vector3d punto2(coords(1, 0), coords(1, 1), coords(1, 2));
            double quadratoDistanza = (punto2 - punto1).squaredNorm();
            fileTracceParaview << fixed << scientific << setprecision(16) << quadratoDistanza << endl;
        }

        fileTracceParaview << "SCALARS SegmentId int 1" << endl;
        fileTracceParaview << "LOOKUP_TABLE default" << endl;
        for(unsigned int i = 0; i < Trac.coordinateIntersezioniTracce.size(); i++){
            fileTracceParaview << i << endl;
        }

        fileTracceParaview.close();
        return true;
    }
}
}
