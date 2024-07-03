#pragma once

#include "GeometryLibrary.hpp"
#include "Utils.hpp"
#include <gtest/gtest.h>
#include <vector>
#include "Eigen/Eigen"
#include "cmath"
#include "algorithm"


using namespace std;
using namespace GeometryLibrary;

// test funzione che calcola il centroide
TEST(FindBarycentreTest, quadrato) {
    vector<Vector3d> quadrato = {
        Vector3d(1, 1, 0),
        Vector3d(1, 2, 0),
        Vector3d(2, 2, 0),
        Vector3d(2, 1, 0)
    };
    Vector3d risultato(1.5, 1.5, 0);
    ASSERT_EQ(FindBarycentre(quadrato), risultato);
}


TEST(FindBarycentreTest, rettangolo) {
    vector<Vector3d> rettangolo = {
        Vector3d(1, 1, 1),
        Vector3d(1, 3, 2),
        Vector3d(4, 3, 3),
        Vector3d(4, 1, 2)
    };
    Vector3d risultato2(2.5, 2, 2);
    ASSERT_EQ(FindBarycentre(rettangolo), risultato2);
}

TEST(FindBarycentreTest, pentagono) {
    vector<Vector3d> pentagono = {
        Vector3d(1, 2, 3),
        Vector3d(2, 4, 6),
        Vector3d(5, 7, 8),
        Vector3d(6, 2, 5),
        Vector3d(3, 1, 4)
    };
    Vector3d risultato3(3.4, 3.2, 5.2);
    ASSERT_EQ(FindBarycentre(pentagono), risultato3);
}

// test per raggio della circonferenza che inscrive una frattura
TEST(CircleRadiusTest, rettangolo) {
    vector<Vector3d> rettangolo = {
       Vector3d(1, 1, 1),
       Vector3d(1, 3, 2),
       Vector3d(4, 3, 3),
       Vector3d(4, 1, 2)
    };
    double raggio(2.06155);
    ASSERT_TRUE(std::abs(CircleRadius(rettangolo)-raggio)<1.e-5*abs(raggio));
}

// test per calcolo della normale ad un piano
TEST(NormalToPlaneTest, triangoloPiano) {
    vector<Vector3d> triangolo = {
        Vector3d(0, 0, 0),
        Vector3d(1, 0, 0),
        Vector3d(0, 1, 0)
    };
    Vector3d normale(0, 0, -1);
    ASSERT_EQ(NormalToPlane(triangolo), normale);
}

TEST(FirstSelectionTracesTest, PoligoniVicini){
    vector<Vector3d> poligono1 = {
      Vector3d(1, 0, 0),
      Vector3d(0, 1, 0),
      Vector3d(0, 0, 0),
      Vector3d(1, 1, 0)
    };
    vector<Vector3d> poligono2 = {
      Vector3d(0.7, -0.2, 0.3),
      Vector3d(0, 0, 0),
      Vector3d(-0.1, 0.5, 0.7),
      Vector3d(0, 0.2, 0.9)
    };
    const double tolleranza1 = 0.5;
    ASSERT_FALSE(FirstSelectionTraces(poligono1, poligono2, tolleranza1));
}

TEST(FirstSelectionTracesTest, PoligoniLontani){
    vector<Vector3d> poligono3 = {
      Vector3d(1, 0, 0),
      Vector3d(0, 1, 0),
      Vector3d(0, 0, 0),
      Vector3d(1, 1, 0)
    };
    vector<Vector3d> poligono4 = {
        Vector3d(3, 2, 4),
        Vector3d(4.8, 2, 4),
        Vector3d(3.2, 3, 4),
        Vector3d(3, 3, 4),
        Vector3d(3, 2.5, 3.5),
    };
    const double tolleranza1 = 0.5;
    ASSERT_TRUE(FirstSelectionTraces(poligono3, poligono4, tolleranza1));
}

TEST(pointonlineTest, TestTrue){
    Vector3d p1(1, 1, 1);
    Vector3d p2(0, 0, 0);
    Vector3d p3(0.8, 0.8, 0.8);
    ASSERT_TRUE(point_on_line(p1, p2, p3));
}

TEST(pointonlineTest, TestFalse){
    Vector3d p1(0, 0, 0);
    Vector3d p2(1, 0, 0);
    Vector3d p3(0, 1, 1);
    ASSERT_FALSE(point_on_line(p1, p2, p3));
}

TEST(checkbarycentrecoordTest, TriangoloTrue){
    Vector3d p1(-0.5, 1, 3.5);
    Vector3d p2(1, 2, 3);
    Vector3d p3(-2, 0, 4);
    Vector3d p4(0, -1, 2);
    ASSERT_TRUE(check_barycentre_coord(p1, p2, p3, p4));
}

TEST(checkbarycentrecoordTest, TriangoloFalse){
    Vector3d p1(3, 3, 3);
    Vector3d p2(1, 2, 3);
    Vector3d p3(-2, 0, 4);
    Vector3d p4(0, -1, 2);
    ASSERT_FALSE(check_barycentre_coord(p1, p2, p3, p4));
}

TEST(checkinsidefractureTest, PoligonoInternoTrue){
    vector<Vector3d> poligono = {
        Vector3d(0, 0, 0),
        Vector3d(1, 0, 0),
        Vector3d(1, 1, 0),
        Vector3d(0, 1, 0)
    };
    Vector3d punto(0.2, 0.2, 0);
    ASSERT_TRUE(check_inside_fracture(punto, poligono));
}

TEST(PuntoInternoPoligonoTest, PoligonoInternoFalse){
    vector<Vector3d> poligono = {
        Vector3d(0, 0, 0),
        Vector3d(1, 0, 0),
        Vector3d(1, 1, 0),
        Vector3d(0, 1, 0)
    };
    Vector3d punto(2, 1.5, 0.3);
    ASSERT_FALSE(check_inside_fracture(punto, poligono));
}

TEST(PianiParalleliTest, pianiparalleli){
    vector<Vector3d> poligono1 = {
        Vector3d(0, 0, 0),
        Vector3d(2, 0, 0),
        Vector3d(2, 1, 0),
        Vector3d(0, 1, 0)
    };
    vector<Vector3d> poligono2 = {
        Vector3d(4, 0, 0),
        Vector3d(6, 0, 0),
        Vector3d(6, 1, 0),
        Vector3d(4, 1, 0)
    };
    ASSERT_TRUE(parallel_planes(poligono1,poligono2));
}

TEST(PianiParalleliTest, pianinonparalleli){
    vector<Vector3d> poligono1 = {
        Vector3d(0, 0, 0),
        Vector3d(2, 0, 0),
        Vector3d(2, 1, 0),
        Vector3d(0, 1, 0)
    };
    vector<Vector3d> poligono2 = {
        Vector3d(0, 1, 2),
        Vector3d(3, 0, 2),
        Vector3d(2, 1, 3),
        Vector3d(0, 0, 3)
    };
    ASSERT_FALSE(parallel_planes(poligono1,poligono2));
}

TEST(CalcolaLunghezzaTracceTest, TestLunghezza) {
    GeometryDFN DFN;
    DFN.Number_Traces = 2;
    //DFN.Traces_Tips[2] = {Vector2i{0, false}, Vector2i{1, true}};
    DFN.Traces_Coordinates[0] = {Vector3d{0, 0, 0}, Vector3d{3, 4, 0}};
    DFN.Traces_Coordinates[1] = {Vector3d{1, 1, 1}, Vector3d{4, 5, 5}};

    calcolaLunghezzaTracce(DFN);

    ASSERT_EQ(DFN.traces_length.size(), 2);
    EXPECT_DOUBLE_EQ(DFN.traces_length[0], 5.0);
    EXPECT_DOUBLE_EQ(DFN.traces_length[1], sqrt(41.0));
}

TEST(CalcolaLunghezzaTracceTest, TestLunghezzaZero) {
    GeometryDFN DFN;
    DFN.Number_Traces = 2;
    DFN.Traces_Coordinates[0] = {Vector3d{1, 1, 1}, Vector3d{1, 1, 1}};
    DFN.Traces_Coordinates[1] = {Vector3d{1, 1, 1}, Vector3d{4, 5, 5}};

    calcolaLunghezzaTracce(DFN);

    ASSERT_EQ(DFN.traces_length.size(), 2);
    EXPECT_DOUBLE_EQ(DFN.traces_length[0], 0.0);  // La lunghezza dovrebbe essere zero
    EXPECT_DOUBLE_EQ(DFN.traces_length[1], sqrt(41.0));
}

TEST(RiordinaLunghezzaTracceTest, TestBase) {
    GeometryDFN DFN;
    DFN.Number_Traces = 3;
    DFN.Traces_Coordinates[0] = {Vector3d{0, 0, 0}, Vector3d{3, 0, 0}};  // Lunghezza = 3.0
    DFN.Traces_Coordinates[1] = {Vector3d{0, 0, 0}, Vector3d{0, 4, 0}};  // Lunghezza = 4.0
    DFN.Traces_Coordinates[2] = {Vector3d{0, 0, 0}, Vector3d{0, 0, 5}};  // Lunghezza = 5.0

    calcolaLunghezzaTracce(DFN);

    const vector<unsigned int> traceIds = {0, 1, 2};
    vector<unsigned int> sortedTraceIds = riordinaLunghezzaTracce(traceIds, DFN);

    ASSERT_EQ(sortedTraceIds.size(), 3);
    EXPECT_EQ(sortedTraceIds[0], 2);
    EXPECT_EQ(sortedTraceIds[1], 1);
    EXPECT_EQ(sortedTraceIds[2], 0);
}

TEST(TrovaTracceTest, Test) {
    GeometryDFN DFN;
    DFN.Number_Fractures = 2;
    DFN.Fractures_Vertices[0] = {Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(1, 1, 0), Vector3d(0, 1, 0)};
    DFN.Fractures_Vertices[1] = {Vector3d(0.5, 0.25, -1), Vector3d(0.5, 0.25, 1), Vector3d(1.5, 0.75, 1), Vector3d(1.5, 0.75, -1)};
    DFN.Fractures_Number_Vertices = {4, 4};
    DFN.Number_Traces = 0;
    DFN.Fractures_Id = {0, 1};

    FindTraces(DFN);

    array<Vector3d, 2> coordinateTraccia = {Vector3d(1.0, 0.5, 0), Vector3d(0.5, 0.25, 0)};
    EXPECT_EQ(DFN.Number_Traces, 1);
    EXPECT_EQ(DFN.Traces_Id.size(), 1);
    EXPECT_EQ(DFN.Traces_Id[0], 0);
    EXPECT_EQ(DFN.Traces_Generator_Id.size(), 1);
    EXPECT_EQ(DFN.Traces_Generator_Id[0], Vector2i(0, 1));
    EXPECT_EQ(DFN.Traces_Coordinates[0], coordinateTraccia);
}

TEST(CalcolaTipologiaTracceTest, TestNonPassante) {
    GeometryDFN DFN;
    DFN.Number_Fractures = 2;
    DFN.Number_Traces = 1;
    DFN.Fractures_Number_Vertices = {4, 4};

    DFN.Traces_Coordinates[0] = {Vector3d(0, 0.5, 0), Vector3d(0.3161837, 0.5, 0)};
    DFN.Fractures_Vertices[0] = {Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(1, 1, 0), Vector3d(0, 1, 0)};
    DFN.Fractures_Vertices[1] = {Vector3d(-0.237778, 0.5, -0.34444),
                                 Vector3d(0.3161837, 0.5, -0.34444),
                                 Vector3d(0.3161837, 0.5, 0.4528389),
                                 Vector3d(-0.237778, 0.5, 0.4528389)};

    DFN.Traces_Generator_Id.push_back(Vector2i(0, 1));

    calcolaTipologiaTracce(DFN);

    vector<int> Id_FrattureConTraccia_attese = {0, 1};

    ASSERT_EQ(DFN.Traces_Tips.size(), DFN.Number_Traces);
    for (const auto& tip : DFN.Traces_Tips) {
        EXPECT_EQ(tip.second[0], true);
        EXPECT_EQ(tip.second[1], true);
    }
}

TEST(CalcolaTipologiaTracceTest, TestPassante) {
    GeometryDFN DFN;
    DFN.Number_Fractures = 2;
    DFN.Number_Traces = 1;
    DFN.Fractures_Number_Vertices = {4, 4};

    DFN.Traces_Coordinates[0] = {Vector3d(0.8, 0, 0), Vector3d(0.8, 1, 0)};
    DFN.Fractures_Vertices[0] = {Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(1, 1, 0), Vector3d(0, 1, 0)};
    DFN.Fractures_Vertices[1] = {Vector3d(0.8, 0, -1), Vector3d(0.8, 0, 0.3), Vector3d(0.8, 1, 0.3), Vector3d(0.8, 1, -1)};

    DFN.Traces_Generator_Id.push_back(Vector2i(0, 1));

    calcolaTipologiaTracce(DFN);

    vector<int> Id_FrattureConTraccia_attese = {0, 1};

    // Verifica
    ASSERT_EQ(DFN.Traces_Tips.size(), DFN.Number_Traces);
    for (const auto& tip : DFN.Traces_Tips) {
        EXPECT_EQ(tip.second[0], false);
        EXPECT_EQ(tip.second[1], false);
    }
}
