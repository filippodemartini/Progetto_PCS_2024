#include <iostream>
#include <vector>
#include "GeometryLibrary.hpp"
#include "Utils.hpp"

using namespace std;
using namespace Eigen;
using namespace GeometryLibrary;

int main()
{
    GeometryDFN dfn;
    string filename = "./FR3_data.txt";
    if(!ImportFractures(filename, dfn)){
        return 1;
    }

    cout << endl;

    FindTraces(dfn);
    for(unsigned int i = 0; i<dfn.Number_Traces; i++)
    {
        cout << " ID: " << dfn.Traces_Id[i] << " Id intersecanti: " << dfn.Traces_Generator_Id[i][0] << ";" << dfn.Traces_Generator_Id[i][1]
             << " Coordinate: " << dfn.Traces_Coordinates[i][0] << ";" << dfn.Traces_Coordinates[i][1] << endl;
    }

    calcolaTipologiaTracce(dfn);
    map<unsigned int, vector<string>> fratture_tracce;

    for (unsigned int k = 0; k < dfn.Number_Traces; k++) {
        for (unsigned int i = 0; i < dfn.Number_Fractures; i++) {
            if (i == dfn.Traces_Generator_Id[k][0] || i == dfn.Traces_Generator_Id[k][1]) {
                string tipo_traccia;
                if (dfn.Traces_Tips[k][0]) {
                    tipo_traccia = "non passante";
                } else {
                    tipo_traccia = "passante";
                }
                fratture_tracce[i].push_back("traccia tipo " + tipo_traccia);
            }
        }
    }

    for (const auto& entry : fratture_tracce) {
        unsigned int frattura_id = entry.first;
        const vector<string>& tracce = entry.second;
        cout << "Frattura " << frattura_id << ": ";
        for (const string& traccia : tracce) {
            cout << traccia << " ";
        }
        cout << endl;
    }

    calcolaLunghezzaTracce(dfn);

    cout << "Lunghezze delle tracce:" << endl;
    for (const auto& lunghezza : dfn.traces_length) {
        cout << lunghezza << endl;
    }

    // calcolaLunghezzaTracce(dfn);
    // auto sorted_traces=dfn.Traces_Tips;
    // cout << "Tracce passanti: " << endl;
    // for (const auto& trace: sorted_traces){
    //     if (!trace.second[1]){
    //         unsigned int id =trace.first;
    //         double length=dfn.traces_length[id];
    //         cout <<"Traccia ID: " << id << ", Lunghezza: " << length << endl;
    //     }
    // }
    // cout << "Tracce non passanti:" << endl;
    // for (const auto& trace : sorted_traces) {
    //     if (trace.second[1]) {
    //         unsigned int id = trace.first;
    //         double length = dfn.traces_length[id];
    //         cout << "Traccia ID: " << id << ", Lunghezza: " << length << endl;
    //     }
    // }



    // PER I FILE DI OUTPUT
    // string fileOutputTracce = "./Tracce_FR10.txt";
    // OutputTracce(dfn, fileOutputTracce);

    // string fileOutputFratture = "./Fratture_FR10.txt";
    // OutputFratture(dfn,fileOutputFratture);



    // calcolaLunghezzaTracce(dfn);

    // cout << "Risultati di Traces_Tips dopo il calcolo delle lunghezze e il riordino:" << endl;
    // for (const auto& entry : dfn.Traces_Tips)
    // {
    //     cout << "ID Traccia: " << entry.first << ", Passante: " << entry.second[0] << ", Non Passante: " << entry.second[1] << endl;
    // }

    // for(unsigned int i = 0; i<dfn.Number_Fractures; i++){
    // Vector3d barycentre = FindBarycentre(dfn.Fractures_Vertices[i]);

    //     cout << "Barycentre : (" << barycentre[0] << ";" << barycentre[1] << ";" << barycentre[2] << ")" << endl;
    // }

    // bool int01 = FirstSelectionTraces(dfn.Fractures_Vertices[0], dfn.Fractures_Vertices[1],numeric_limits<double>::epsilon());
    // bool int02 = FirstSelectionTraces(dfn.Fractures_Vertices[0], dfn.Fractures_Vertices[2],numeric_limits<double>::epsilon());
    // bool int12 = FirstSelectionTraces(dfn.Fractures_Vertices[1], dfn.Fractures_Vertices[2],numeric_limits<double>::epsilon());

    // cout << int01 << "   " << int02 << "   " << int12 << endl;

  return 0;
}
