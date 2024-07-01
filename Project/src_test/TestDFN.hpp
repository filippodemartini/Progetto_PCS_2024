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

