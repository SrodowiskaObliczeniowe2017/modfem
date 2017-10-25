#include "Enums.h"
#include "elem_tables.h"
#include "Point3D.h"
#include<stdio.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <assert.h>
#include <string.h>
using namespace FemViewer;

const int iTetraEdgeConnect[6][2] = {
	{0,1}, {1,2}, {0,2}, {0,3}, {1,3}, {2,3}
};

const int iTetraVertcesOfEdges[4][3] = { // celowo zmieniona podstawa
	{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}
};

int iTetraEdgeTbl[16]=
  {
    0x000, /*0000 0000 inactive*/
    0x00d, /*0000 1101*/
    0x013, /*0001 0011*/
    0x01e, /*0001 1110*/
    0x026, /*0010 0110*/
    0x02b, /*0010 1011*/
    0x035, /*0011 0101*/
    0x038, /*0011 1000*/
    0x038, /*0011 1000*/
    0x035, /*0011 0101*/
    0x02b, /*0010 1011*/
    0x026, /*0010 0110*/
    0x01e, /*0001 1110*/
    0x013, /*0001 0011*/
    0x00d, /*0000 1101*/
    0x000, /*0000 0000 active*/
  };


/*Tetrahedron Triangles Connection Table*/

/*
This table lists all configurations for triangulation of edge intersection points for tetrahedon.
For each of the case there is a specific triangulation in the form
of 0-5 egde tiples terminated by -1. There is maximum 2 trinagles per tetrahedon in 4
topologies. There are no ambiguous cases. 
*/

int iTetrahedronTrianglesCT[16][7] =
  {
    {-1, -1, -1, -1, -1, -1, -1}, /* 0 */
    { 0,  3,  2, -1, -1, -1, -1}, /* 1 */
    { 0,  1,  4, -1, -1, -1, -1}, /* 2 */
    { 1,  4,  3,  1,  3,  2, -1}, /* 3 */

    { 1,  2,  5, -1, -1, -1, -1}, /* 4 */
    { 0,  3,  5,  0,  5,  1, -1}, /* 5 */
    { 0,  2,  5,  0,  5,  4, -1}, /* 6 */
    { 5,  4,  3, -1, -1, -1, -1}, /* 7 */

    { 3,  4,  5, -1, -1, -1, -1}, /* 8 */
    { 4,  5,  0,  5,  2,  0, -1}, /* 9 */
    { 1,  5,  0,  5,  3,  0, -1}, /* 10 */
    { 5,  2,  1, -1, -1, -1, -1}, /* 11 */

    { 2,  3,  1,  3,  4,  1, -1}, /* 12 */
    { 4,  1,  0, -1, -1, -1, -1}, /* 13 */
    { 2,  3,  0, -1, -1, -1, -1}, /* 14 */
    {-1, -1, -1, -1, -1, -1, -1}, /* 15 */
  };

int iSideTrainglesCT[16][4][6] = 
{
	{ {-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1} }, /*0*/
	{ { 0, 6, 4,-1,-1,-1}, { 0, 4, 7,-1,-1,-1}, { 0, 7, 6,-1,-1,-1}, {-1,-1,-1,-1,-1,-1} }, /*1*/
	{ { 1, 4, 5,-1,-1,-1}, { 8, 4, 1,-1,-1,-1}, {-1,-1,-1,-1,-1,-1}, { 1, 5, 8,-1,-1,-1} }, /*2*/
	{ { 6, 5, 1, 6, 1, 0}, { 8, 7, 0, 8, 0, 1}, { 7, 6, 0,-1,-1,-1}, { 1, 5, 8,-1,-1,-1} }, /*3*/
	{ { 2, 5, 6,-1,-1,-1}, {-1,-1,-1,-1,-1,-1}, { 2, 6, 9,-1,-1,-1}, { 9, 5, 2,-1,-1,-1} }, /*4*/
	{ { 0, 2, 5, 0, 5, 4}, { 0, 4, 7,-1,-1,-1}, { 7, 9, 2, 7, 2, 0}, { 2, 9, 5,-1,-1,-1} }, /*5*/
	{ { 6, 2, 1, 6, 1, 4}, { 8, 4, 1,-1,-1,-1}, { 2, 6, 9,-1,-1,-1}, { 8, 1, 2, 8, 2, 9} }, /*6*/
	{ { 0, 2, 1,-1,-1,-1}, { 7, 0, 1, 7, 1, 8}, { 0, 7, 9, 0, 9, 2}, { 8, 1, 2, 8, 2, 9} }, /*7*/
	{ {-1,-1,-1,-1,-1,-1}, { 3, 7, 8,-1,-1,-1}, { 9, 7, 3,-1,-1,-1}, { 3, 8, 9,-1,-1,-1} }, /*8*/
	{ { 6, 4, 0,-1,-1,-1}, { 4, 8, 3, 4, 3, 0}, { 9, 6, 0, 9, 0, 3}, { 3, 8, 9,-1,-1,-1} }, /*9*/
	{ { 1, 4, 5,-1,-1,-1}, { 1, 3, 7, 1, 7, 4}, { 9, 7, 3,-1,-1,-1}, { 5, 9, 3, 5, 3, 1} }, /*10*/
	{ { 6, 5, 1, 6, 1, 0}, { 0, 1, 3,-1,-1,-1}, { 0, 3, 9, 0, 9, 6}, { 5, 9, 3, 5, 3, 1} }, /*11*/
	{ { 2, 5, 6,-1,-1,-1}, { 3, 7, 8,-1,-1,-1}, { 6, 7, 3, 6, 3, 2}, { 2, 3, 8, 2, 8, 5} }, /*12*/
	{ { 5, 4, 0, 5, 0, 2}, { 4, 8, 3, 4, 3, 0}, { 0, 3, 2,-1,-1,-1}, { 2, 3, 8, 2, 8, 5} }, /*13*/
	{ { 4, 6, 2, 4, 2, 1}, { 1, 3, 7, 1, 7, 4}, { 6, 7, 3, 6, 3, 2}, { 1, 2, 3,-1,-1,-1} }, /*14*/
	{ { 0, 2, 1,-1,-1,-1}, { 0, 1, 3,-1,-1,-1}, { 0, 3, 2,-1,-1,-1}, { 1, 2, 3,-1,-1,-1} }  /*15*/
};


/*Squares Triangles Connection Table*/



/* Prism */
/* 9-bit key */
const static int iPrizmEdgeTbl[64] = {
	/*1*/0x0,  	0x045,	0x083,	0x0c6,	0x106,	0x143,	0x185,	0x1c0,
	/*2*/0x068,	0x02d,	0x0eb,	0x03e,	0x16e,	0x11d,	0x1ed,	0x1a8,
	/*3*/0x098,	0x0dd,	0x01b,	0x073,	0x19e,	0x1db,	0x12b,	0x158,
	/*4*/0x0f0,	0x0ae,	0x05e,	0x036,	0x1f6,	0x1b3,	0x175,	0x130,
	/*5*/0x130,	0x175,	0x1b3,	0x1f6,	0x036,	0x05e,	0x0ae,	0x0f0,
	/*6*/0x158,	0x12b,	0x1db,	0x19e,	0x073,	0x01b,	0x0dd,	0x098,
	/*7*/0x1a8,	0x1ed,	0x11d,	0x16e,	0x03e,	0x0eb,	0x02d,	0x068,
	/*8*/0x1c0, 0x185,  0x143,	0x106,	0x0c6,	0x083,	0x045,	0x0,
};

#define X -1
/* The first pozition in each row dinote tne number of triangles */
const static int iPrizmtriTable[64][10] = {
	/*1*/
	{ 0, X, X, X, X, X, X, X, X, X },
	{ 1, 6, 2, 0, X, X, X, X, X, X },
	{ 1, 0, 1, 7, X, X, X, X, X, X },
	{ 2, 2, 1, 7, 7, 6, 2, X, X, X },
	{ 1, 1, 2, 8, X, X, X, X, X, X },
	{ 2, 1, 0, 6, 6, 8, 1, X, X, X },
	{ 2, 0, 2, 8, 8, 7, 0, X, X, X },
	{ 1, 8, 7, 6, X, X, X, X, X, X },
	/*2*/
	{ 1, 3, 5, 6, X, X, X, X, X, X },
	{ 2, 0, 3, 5, 5, 2, 0, X, X, X },
	{ 2, 0, 1, 7, 3, 5, 6, X, X, X },
	{ 3, 1, 3, 4, 1, 4, 5, 5, 2, 1 },
	{ 2, 1, 2, 8, 3, 5, 6, X, X, X },
	{ 3, 0, 3, 4, 0, 4, 8, 0, 8, 2 },
	{ 3, 0, 2, 8, 8, 7, 0, 3, 5, 6 },
	{ 2, 3, 5, 8, 8, 7, 3, X, X, X },
	/*3*/
	{ 1, 7, 4, 3, X, X, X, X, X, X },
	{ 2, 6, 2, 0, 7, 4, 3, X, X, X },
	{ 2, 0, 1, 4, 4, 3, 0, X, X, X },
	{ 3, 0, 1, 4, 4, 6, 0, 6, 4, 5 },
	{ 2, 1, 2, 8, 7, 4, 3, X, X, X },
	{ 3, 7, 4, 3, 0, 6, 8, 8, 1, 0 },
	{ 3, 0, 1, 3, 3, 1, 8, 8, 5, 3 },
	{ 2, 8, 4, 3, 3, 6, 8, X, X, X },
	/*4*/
	{ 2, 4, 5, 6, 6, 7, 4, X, X, X },
	{ 3, 5, 2, 3, 3, 2, 1, 1, 7, 3 },
	{ 3, 6, 2, 3, 3, 2, 1, 1, 4, 3 },
	{ 2, 1, 4, 5, 5, 2, 1, X, X, X },
	{ 3, 1, 2, 8, 4, 5, 6, 6, 7, 4 },
	{ 2, 7, 1, 0, 4, 5, 8, X, X, X },
	{ 2, 0, 2, 6, 4, 5, 8, X, X, X },
	{ 1, 4, 5, 8, X, X, X, X, X, X },
	/*5*/
	{ 1, 8, 5, 4, X, X, X, X, X, X },
	{ 2, 8, 5, 4, 6, 2, 0, X, X, X },
	{ 2, 8, 5, 4, 0, 1, 7, X, X, X },
	{ 3, 4, 7, 6, 6, 5, 4, 8, 2, 1 },
	{ 2, 1, 2, 5, 5, 4, 1, X, X, X },
	{ 3, 3, 4, 1, 1, 2, 3, 2, 6, 3 },
	{ 3, 1, 2, 7, 3, 7, 2, 5, 3, 2 },
	{ 2, 4, 7, 6, 6, 5, 4, X, X, X },
	/*6*/
	{ 2, 8, 6, 3, 3, 4, 8, X, X, X },
	{ 3, 0, 3, 1, 3, 8, 1, 8, 3, 5 },
	{ 3, 0, 1, 8, 8, 6, 0, 3, 4, 7 },
	{ 2, 3, 4, 7, 8, 2, 1, X, X, X },
	{ 3, 0, 4, 1, 4, 0, 6, 6, 5, 4 },
	{ 2, 0, 3, 4, 4, 1, 0, X, X, X },
	{ 2, 3, 4, 7, 0, 2, 6, X, X, X },
	{ 1, 3, 4, 7, X, X, X, X, X, X },
	/*7*/
	{ 2, 3, 7, 8, 8, 5, 3, X, X, X },
	{ 3, 6, 5, 3, 0, 7, 8, 8, 2, 0 },
	{ 3, 0, 2, 3, 3, 2, 8, 8, 4, 3 },
	{ 2, 6, 5, 3, 8, 2, 1, X, X, X },
	{ 3, 1, 2, 5, 5, 4, 1, 4, 3, 1 },
	{ 2, 6, 5, 3, 7, 1, 0, X, X, X },
	{ 2, 0, 2, 5, 5, 3, 0, X, X, X },
	{ 1, 0, 2, 6, X, X, X, X, X, X },
	/*8*/
	{ 1, 6, 7, 8, X, X, X, X, X, X },
	{ 2, 0, 7, 8, 8, 2, 0, X, X, X },
	{ 2, 1, 8, 6, 6, 0, 1, X, X, X },
	{ 1, 8, 2, 1, X, X, X, X, X, X },
	{ 2, 2, 6, 7, 7, 1, 2, X, X, X },
	{ 1, 7, 1, 0, X, X, X, X, X, X },
	{ 1, 0, 2, 6, X, X, X, X, X, X },
	{ 0, X, X, X, X, X, X, X, X, X }
};

#undef X


/* maps ede ID to coresponding vertixes */
const int iPrismEdges[9][2] = {
	/*0*/{ 0, 1 },
	/*1*/{ 1, 2 },
	/*2*/{ 0, 2 },
	/*3*/{ 3, 4 },
	/*4*/{ 4, 5 },
	/*5*/{ 3, 5 },
	/*6*/{ 0, 3 },
	/*7*/{ 1, 4 },
	/*8*/{ 2, 5 }
};

/* maps prism vertices into edges */
const int iPrismVerticesOfEdges[6][3] = {
	{0,2,6}, {0,1,7}, {1,2,8}, {3,5,6}, {3,4,7}, {4,5,8}
};

/* maps prism edges into faces */
const int iPrismEdgesOfFaces[9][2] = {
	{0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,4}, {2,4}, {2,3}, {3,4}
};

int iQuadEdgesIds[4][2] = {
	{0,1},{1,2},{2,3},{3,0}
};

int iSLinesCT[16][5]=
  {
    {-1, -1, -1, -1, -1},
    {0, 3, -1, -1, -1},
    {0, 1, -1, -1, -1},
    {3, 1, -1, -1, -1},
    {2, 1, -1, -1, -1},
    {0, 1,  3,  2, -1}, /* two cases */
    {2, 0, -1, -1, -1},
    {2, 3, -1, -1, -1},
    {2, 3, -1, -1, -1},
    {2, 0, -1, -1, -1},
    {1, 0,  3,  2, -1}, /* two cases */
    {2, 1, -1, -1, -1},
    {3, 1, -1, -1, -1},
    {0, 1, -1, -1, -1},
    {3, 0, -1, -1, -1},
    {-1, -1, -1, -1, -1}
  };

/* Table of flags of edges for square */
int iQuadFlags[16] = {
	0x00,
	0x09,
	0x03,
	0x0a,
	0x06,
	0x0f,
	0x05,
	0x0c,
	0x0c,
	0x05,
	0x0f,
	0x06,
	0x0a,
	0x03,
	0x09,
	0x00
};


int iQuadEdges[16][7] = {
	{0,-1,-1,-1,-1,-1,-1}, /*0*/
	{3, 0, 4, 7,-1,-1,-1}, /*1*/
	{3, 1, 5, 4,-1,-1,-1}, /*2*/
	{4, 0, 1, 5, 7,-1,-1}, /*3*/
	{3, 2, 6, 5,-1,-1,-1}, /*4*/
	{6, 0, 4, 7, 2, 6, 5}, /*5*/
	{4, 1, 2, 6, 4,-1,-1}, /*6*/
	{5, 0, 1, 2, 6, 7,-1}, /*7*/
	{3, 3, 7, 6,-1,-1,-1}, /*8*/
	{4, 0, 4, 6, 3,-1,-1}, /*9*/
	{6, 1, 5, 4, 3, 7, 6}, /*10*/
	{5, 0, 1, 5, 6, 3,-1}, /*11*/
	{4, 2, 3, 7, 5,-1,-1}, /*12*/
	{5, 0, 4, 5, 2, 3,-1}, /*13*/
	{5, 1, 2, 3, 7, 4,-1}, /*14*/
	{4, 0, 1, 2, 3,-1,-1}, /*15*/
};

const int iPrismFaces[5][4] = {{0,2,1,-1},{3,4,5,-1},{0,1,4,3},{1,2,5,4},{0,3,5,2}};
const int iPrizmMidPts[5][4] = {{6,7,8,-1},{9,10,11,-1},{6,13,9,12},{7,14,10,13},{8,12,11,14}};


/*
int  DoSliceTerahedron(const FemViewer::CutPlane* plane,const int faces[5],double tetra[16], std::vector<FemViewer::GraphElement2<double> >& grEls,FILE* fp)
{
	int key_id(0),ed_flg,v0,v1,m;
	double coor[60],u;  // moze byc wolne?

	::memcpy(coor,tetra,sizeof(double) * 16);

	for (int i=0;i<4;i++)
	{
		register int m = (i << 2);
		if (plane->CheckLocation(coor[m],coor[m+1],coor[m+2]) < 1)
			key_id |= (1<<i);
	}

	ed_flg = iTetraEdgeTbl[key_id];
	if (ed_flg == 0) return 0;

	for (int i=0;i<6;++i)
	{
		if (ed_flg & (1 << i))
		{
			register int m = (i << 2) + 16;
			v0 = (iTetraEdgeConnect[i][0] << 2);
			v1 = (iTetraEdgeConnect[i][1] << 2);
		    //plane->IntersectWithLine(&coor[v1],&coor[v0],&coor[m],u);
			coor[m+3]  = coor[v0+3];
			coor[m+3] += u*(coor[v1+3] - coor[v0+3]);
		}
	}
	//std::cout<<"w elemtables>>>>>>>\n";
	int cnt(0);
	//double out[3];
	FemViewer::Point3D<double> pt;
	FemViewer::GraphElement2<double> triangle;
	for (int i(0);i<2;++i)
	{
		v0 = 3*i;
		if (iTetrahedronTrianglesCT[key_id][v0]<0) break;
		
		triangle.Clear();
		register int m; 
		if (key_id >> 3)
		{
			v1  = iTetrahedronTrianglesCT[key_id][v0+2];
			m = (v1 << 2) + 16; 
			pt.x = coor[m];
			pt.y = coor[m+1];
			pt.z = coor[m+2]; 
			triangle.Add(pt,coor[m+3]);

			v1  = iTetrahedronTrianglesCT[key_id][v0];
			m = (v1 << 2) + 16; 
			pt.x = coor[m];
			pt.y = coor[m+1];
			pt.z = coor[m+2]; 
			triangle.Add(pt,coor[m+3]);

			v1  = iTetrahedronTrianglesCT[key_id][v0+1];
			m = (v1 << 2) + 16; 
			pt.x = coor[m];
			pt.y = coor[m+1];
			pt.z = coor[m+2]; 
			triangle.Add(pt,coor[m+3]);
		}	
		else for (int j(0);j<3;++j)
		{
			v1  = iTetrahedronTrianglesCT[key_id][v0+j];
			m = (v1 << 2) + 16; 
			pt.x = coor[m];
			pt.y = coor[m+1];
			pt.z = coor[m+2]; 
			triangle.Add(pt,coor[m+3]);
		}

		if(!triangle.Check()) {
			grEls.push_back(triangle);
			++cnt;
		} else {
			std::cout<<"cos jest zle w tetra slice\n";
		}
		
	}

	if (faces == NULL) return cnt;

	int v3;
	//FemViewer::Point3D<double> pt1,pt2,pt3,pt4;
	// Do rest of faces
	for (int i(0);i<faces[0];++i)
	{
		//std::cout<<"w II\n";
		v0 = faces[i+1];
		assert(v0 < 4);
		for (int j(0);j<2;++j)
		{
			v1 = 3*j;
			if (iSideTrainglesCT[key_id][v0][v1] < 0) break;
			
			triangle.Clear();
			register int m;
			for (int k(0);k<3;++k)
			{
				v3 = iSideTrainglesCT[key_id][v0][v1+k];
				assert(0 <= v3 && v3 <= 9);
				m = (v3 << 2);
				pt.x = coor[m];
				pt.y = coor[m+1];
				pt.z = coor[m+2]; 
				assert(0.0 <= coor[m] && coor[m] <= 1.0);
				triangle.Add(pt,coor[m+3]);
			}
			if(!triangle.Check()) {
				grEls.push_back(triangle);
				++cnt;
			} else {
				std::cout<<"cos jest zle tetra sides\n";
			}
			
		}// end for
	}

	return cnt;
}


int  DoSlicePrizm(const FemViewer::CutPlane* plane,const int faces[5],double prizm[24], std::vector<FemViewer::GraphElement2<double> >& grEls,FILE* fp)
{
	int key_id(0),ed_flg,v0,v1;
	double coor[60],u;  // moze byc wolne?

	// to trzeba zmienic na zewnetrzna tablice o wiekszej pjemnosci
	::memcpy(coor,prizm,sizeof(double) * 24);

	for (int i(0);i<6;i++)
	{
		register int m = (i << 2);
		if (plane->CheckLocation(coor[m],coor[m+1],coor[m+2]) < 1)
			key_id |= (1<<i);
	}
		
	ed_flg = iPrizmEdgeTbl[key_id];
	if (ed_flg == 0) return 0;

	for (int i(0);i<9;++i)
	{
		register int m = (i << 2) + 24;
		if (ed_flg & (1 << i))
		{
			v0 = iPrismEdges[i][0] << 2;
			v1 = iPrismEdges[i][1] << 2;

			//??plane->IntersectWithLine(&coor[v1],&coor[v0],&coor[m],u);
			m += 3;
			coor[m]  = coor[v0+3];
			coor[m] += u*(coor[v1+3] - coor[v0+3]);
		}
	}

	int cnt(0);
	FemViewer::Point3D<double> pt;
	FemViewer::GraphElement2<double> triangle;
	for (int i(0);i<iPrizmtriTable[key_id][0];++i)
	{
		v0 = 3*i + 1;
		triangle.Clear();

		 // key_id / 8
		for (int j(0);j<3;++j)
		{
			v1  = iPrizmtriTable[key_id][v0+j];
			register int m = (v1 << 2) + 24; 
			pt.x = coor[m];
			pt.y = coor[m+1];
			pt.z = coor[m+2]; 
			triangle.Add(pt,coor[m+3]);
		}
		
		if(!triangle.Check()) {
			grEls.push_back(triangle);
			++cnt;
		}
	}

	if (faces == NULL) return cnt;

	FemViewer::Point3D<double> pt1,pt2,pt3,pt4;
		
	int result[4];
	double p[32];
	for(unsigned int j = 0; j<5 ;++j)
	{
		if (faces[j] != 1) continue;

		if (j < 2)
		{
			register int pa = iPrismFaces[j][0] << 2;
			register int pb = iPrismFaces[j][1] << 2;
			register int pc = iPrismFaces[j][2] << 2;
			
			result[0] = plane->CheckLocation(&coor[pa]);
			result[1] = plane->CheckLocation(&coor[pb]);
			result[2] = plane->CheckLocation(&coor[pc]);
			
			triangle.Clear();

			if (result[0] <= 0 && result[1] <= 0 && result[2] <= 0) {
				pt1.x = static_cast<float>(coor[pa++]);
				pt1.y = static_cast<float>(coor[pa++]);
				pt1.z = static_cast<float>(coor[pa++]);

				pt2.x = static_cast<float>(coor[pb++]);
				pt2.y = static_cast<float>(coor[pb++]);
				pt2.z = static_cast<float>(coor[pb++]);

				pt3.x = static_cast<float>(coor[pc++]);
				pt3.y = static_cast<float>(coor[pc++]);
				pt3.z = static_cast<float>(coor[pc++]);
				
				triangle.Add(pt1,coor[pa]);
				triangle.Add(pt2,coor[pb]);
				triangle.Add(pt3,coor[pc]);
				if(!triangle.Check()) {
					grEls.push_back(triangle);
					++cnt;
				}

				continue;
			}

			if (result[0]*result[1]>=0) {
				std::swap(result[2],result[0]);
				std::swap(pc,pa);
				std::swap(result[1],result[2]);
				std::swap(pb,pc);
			} else if (result[0]*result[2]>=0) {
				std::swap(result[1],result[0]);
				std::swap(pb,pa);
				std::swap(result[1],result[2]);
				std::swap(pb,pc);
			}
				
			// powtorka
			//??plane->IntersectWithLine(&coor[pb],&coor[pa],&p[0],u);
			//??plane->IntersectWithLine(&coor[pc],&coor[pa],&p[3],u);
			pt1.x = static_cast<float>(p[0]);
			pt1.y = static_cast<float>(p[1]);
			pt1.z = static_cast<float>(p[2]);

			pt2.x = static_cast<float>(p[3]);
			pt2.y = static_cast<float>(p[4]);
			pt2.z = static_cast<float>(p[5]);

			p[6]  = coor[pa+3];
			p[7]  = coor[pa+3];
			p[6] += u*(coor[pb+3] - coor[pa+3]);
			p[7] += u*(coor[pc+3] - coor[pa+3]);					
				
			if (result[0] < 0) {
				pt3.x = static_cast<float>(coor[pa++]);
				pt3.y = static_cast<float>(coor[pa++]);
				pt3.z = static_cast<float>(coor[pa++]);
				
				triangle.Add(pt1,p[6]);
				triangle.Add(pt2,p[7]);
				triangle.Add(pt3,coor[pa]);
				
				if(!triangle.Check()) {
					grEls.push_back(triangle);
					++cnt;
				}

			} else {

				pt3.x = static_cast<float>(coor[pb++]);
				pt3.y = static_cast<float>(coor[pb++]);
				pt3.z = static_cast<float>(coor[pb++]);	
							
				pt4.x = static_cast<float>(coor[pc++]);
				pt4.y = static_cast<float>(coor[pc++]);
				pt4.z = static_cast<float>(coor[pc++]);

				triangle.Clear();
				triangle.Add(pt1,p[6]);
				triangle.Add(pt3,coor[pb]);
				triangle.Add(pt2,p[7]);
				
				if(!triangle.Check()) {
					grEls.push_back(triangle);
					++cnt;
				}

				triangle.Clear();
				triangle.Add(pt2,p[7]);
				triangle.Add(pt3,coor[pb]);
				triangle.Add(pt4,coor[pc]);
				
				if(!triangle.Check()) {
					grEls.push_back(triangle);
					++cnt;
				}						

			}

		} else {
			int key_id = 0;
			for (register int k(0); k< 4; ++k) {
				register int pa = (iPrismFaces[j][k] << 2);
				if (plane->CheckLocation(&coor[pa]) < 1)
					key_id |= (1<<k);
			}

			if (key_id == 0) continue;
			int ed_flg = iQuadFlags[key_id];
						
			for (register int k(0);k<4;++k) {
				register int pa = (iPrismFaces[j][k] << 2);
				memcpy(&p[k<<2],&coor[pa],4*sizeof(double));
			}

			for (register int k(0); k<4; ++k) {
				if (ed_flg & (1<<k)) {
					register int m = (k << 2) + 16;
					register int pa = iQuadEdgesIds[k][0] << 2;
					register int pb = iQuadEdgesIds[k][1] << 2;
					//??plane->IntersectWithLine(&p[pa],&p[pb],&p[m],u);
					m  += 3;
					pa += 3;
					pb += 3;
					p[m]  = p[pb];
					p[m] += u*(p[pa] - p[pb]);
				}
			}

			const int *pp = &iQuadEdges[key_id][0];
			triangle.Clear();

			switch(*pp) {
			case 3:
				{
				register int pa = *(++pp);
				pt1.x = static_cast<float>(p[pa++]);
				pt1.y = static_cast<float>(p[pa++]);
				pt1.z = static_cast<float>(p[pa++]);
				
				register int pb = *(++pp);
				pt2.x = static_cast<float>(p[pb++]);
				pt2.y = static_cast<float>(p[pb++]);
				pt2.z = static_cast<float>(p[pb++]);
				
				register int pc = *(++pp);
				pt3.x = static_cast<float>(p[pc++]);
				pt3.y = static_cast<float>(p[pc++]);
				pt3.z = static_cast<float>(p[pc++]);
				
				triangle.Add(pt1,p[pa]);
				triangle.Add(pt1,p[pb]);
				triangle.Add(pt1,p[pc]);
				if(!triangle.Check()) {
					grEls.push_back(triangle);
					++cnt;
				}
			}
				break;
			case 4:
				{
				register int pa = *(++pp) << 2;
				pt1.x = static_cast<float>(p[pa++]);
				pt1.y = static_cast<float>(p[pa++]);
				pt1.z = static_cast<float>(p[pa++]);
				
				register int pb = *(++pp) << 2;
				pt2.x = static_cast<float>(p[pb++]);
				pt2.y = static_cast<float>(p[pb++]);
				pt2.z = static_cast<float>(p[pb++]);
				
				register int pc = *(++pp) << 2;
				pt3.x = static_cast<float>(p[pc++]);
				pt3.y = static_cast<float>(p[pc++]);
				pt3.z = static_cast<float>(p[pc++]);

				register int pd = *(++pp) << 2;
				pt4.x = static_cast<float>(p[pd++]);
				pt4.y = static_cast<float>(p[pd++]);
				pt4.z = static_cast<float>(p[pd++]);
				
				triangle.Add(pt1,p[pa]);
				triangle.Add(pt2,p[pb]);
				triangle.Add(pt3,p[pc]);
				if(!triangle.Check()) {
					grEls.push_back(triangle);
					++cnt;
				}

				triangle.Clear();
				triangle.Add(pt3,p[pc]);
				triangle.Add(pt4,p[pd]);
				triangle.Add(pt1,p[pa]);
				if(!triangle.Check()) {
					grEls.push_back(triangle);
					++cnt;
				}
				}
				break;
			case 5:
				{
				}
				break;
			case 6:
				std::cerr <<"This case shouldn't happend!!\n";
			default:
				break;
			}
			/*for (register int k(1); k<=*p;++k) {
				if (k == *pp) {
					i1 = pp + 1;
					i2 = pp + k;
				} else {
					i1 = pp + k;
					i2 = i1 + 1;
				}

				register int pa = *i1 << 2;
							
				pt1.x = static_cast<float>(p[pa++]);
				pt1.y = static_cast<float>(p[pa++]);
				pt1.z = static_cast<float>(p[pa]);

				pa = *i2 << 2;
				pt2.x = static_cast<float>(p[pa++]);
				pt2.y = static_cast<float>(p[pa++]);
				pt2.z = static_cast<float>(p[pa]);

			
			}*/
/*
	}


	return cnt;
}
*/
/*
int  CutPrism(const double prism_vertices[18],const int nodes[7], FemViewer::CutPlane&  cut_plane,std::vector<myVertex>& out_vertices)
{
	myVertex ed, bg;
	int key_id(0);
    for (int i(0),index = 0;i<6;i++,index += 3) {
    	const double* coor = &prism_vertices[index];
	    if (cut_plane.CheckLocation(*coor++,*coor++,*coor) < 1)
			key_id |= (1<<i);
   	}

    int ed_flg = iPrizmEdgeTbl[key_id];
	if (ed_flg == 0) return 0;

	for (int i(0);i<9;++i) {
		if (ed_flg & (1 << i)) {
			int v0 = iPrismEdges[i][0];
			bg.position.x = prism_vertices[3*v0];
			bg.position.y = prism_vertices[3*v0+1];
			bg.position.z = prism_vertices[v0+2];
		    int v1 = iPrismEdges[i][1];
		    ed.position.x = prism_vertices[3*v1];
		    ed.position.y = prism_vertices[3*v1+1];
		    ed.position.z = prism_vertices[3*v1+2];
		    // 0 <= u <= 1
		    double u;
		    if (cut_plane.IntersectWithLine(ed.position.v,bg.position.v,u)) {
		       	out_vertices.push_back(ed);
		    }
		}
	}
	assert(out_vertices.size()>=3);
	// Sort vertices in CCW order
	myVertex plane_normal = myVertex(cut_plane.GetParams());
	myVertex orig = out_vertices[0];
	std::sort(out_vertices.begin()+1,out_vertices.end(),
			[&](const myVertex& lhs,const myVertex& rhs) -> bool {
				myVertex v = (lhs - orig) * (rhs - orig);
				return FemViewer::fvmath::Dot(v,plane_normal) > 0;
			} );



	return static_cast<int>(out_vertices.size());
}
*/

/* Check for given prism element as a table of coordinates, if is selected(cutted) by a planes
 * Return the index of cutted plane (face) or 0 if cut
 * or -1 if not
 */
int  IsPrismSelected(const double coords[6*4],const double* plane)
{
	int key_id = 0;
	int edges[9] = {0};
	int faces[5] = {0};

	for (int i=0; i<6;i++)
	{
		int pos = Plane::Check(plane,&coords[3*i]);
		if (pos < 1) key_id |= (1<<i);
		if (! pos) {
			int e0 = iPrismVerticesOfEdges[i][0];
			int e1 = iPrismVerticesOfEdges[i][1];
			int e2 = iPrismVerticesOfEdges[i][2];
			edges[e0]++;
			if (edges[e0] == 2) {
				int f0 = iPrismEdgesOfFaces[e0][0];
				int f1 = iPrismEdgesOfFaces[e0][1];
				faces[f0]++;
				faces[f1]++;
			}
			edges[e1]++;
			if (edges[e1] == 2) {
				int f0 = iPrismEdgesOfFaces[e1][0];
				int f1 = iPrismEdgesOfFaces[e1][1];
				faces[f0]++;
				faces[f1]++;
			}
			edges[e2]++;
			if (edges[e2] == 2) {
				int f0 = iPrismEdgesOfFaces[e2][0];
				int f1 = iPrismEdgesOfFaces[e2][1];
				faces[f0]++;
				faces[f1]++;
			}
		}
	}

	int ed_flg = iPrizmEdgeTbl[key_id];
	if (ed_flg == 0) return -1;

	// Loop over faces
	int result  = 0, nf= 0;
	do {
		if (faces[nf] >= 2) { result = nf + 1; break; }
	} while(nf++ < 5);

	return result;
}

/* Check for given terahedron element as a table of coordinates, if is selected(cutted) by a planes
 * Return the index of cutted plane (face) or 0 if it is cut
 * or -1 if not
 */

int IsTetraSelected(const double coords[4*4],const double* plane)
{
	static const int lTVF[4][3] = { {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3} };
	int key_id = 0;
	int faces[4] = {0};
	for (int i=0; i < 4; ++i)
	{
		int pos = Plane::Check(plane,&coords[3*i]);
		if (pos < 1) key_id |= (1< i);
		if (pos == 0) {
			int f0 = lTVF[i][0];
			int f1 = lTVF[i][1];
			int f2 = lTVF[i][2];
			faces[f0]++;
			faces[f1]++;
			faces[f2]++;
		}
	}

	int ed_flg = iPrizmEdgeTbl[key_id];
	if (ed_flg == 0) return -1;

	int result = 0, nf= 0;
	do {
		if (faces[nf] == 3) { result = nf + 1; break; }
	} while(nf++ < 4);

	return result;
}

// Function slices an element ((0)terahedron/(1)prism)
// and returns indices of polygon sorted in CCW order (triangle-fan)
unsigned int  CutElement(
		const double vertices[18],  /* In: element physical coordinates */
		const int nodes[7],			/* In: indices of element's nodes */
		const myPlane*  cut_plane,  /* In: handle to cut-plane */
		std::vector<myVertex>& out_vertices /* Out: vertices obtained from cut */
		)
{
	static const int  *iEdgesTbl[2] = { &iTetraEdgeTbl[0], &iPrizmEdgeTbl[0] };
	static const int (*iEdgesConnect[])[2] = { iTetraEdgeConnect, iPrismEdges };
    static const int   nedges[2] = {6,9};

	myVertex ed, bg, out;
	int key_id(0);

	// Set a type of an element
    const int type(nodes[0] == 6);
    const double* p(cut_plane->GetParams());
    //std::cout << "Cutplane: " << p[0] << " " << p[1] << " " << p[2] << " " << p[3]  << std::endl;
   // mfp_debug("Print element (%d) coordinates: ",type);
    //for (int i=0;i<18;i+=3){
    //	mfp_debug("vertex %d = {%f, %f,%f}\n",i/3,vertices[i],vertices[i+1],vertices[i+2]);
   // }

    out_vertices.clear();

	for (int i=0,index=0; i<nodes[0]; i++,index+=3) {
	   	const double* coor = &vertices[index];
	   	int res;
	    if ((res=cut_plane->CheckLocation(*coor,*(coor+1),*(coor+2))) < 1)
			key_id |= (1<<i);
	    //mfp_debug("location for pt = {%f %f %f} is %d and key is %d",*(coor),*(coor+1),*(coor+2),res,key_id);
	}

	//mfp_debug("Edge key: %d\n",key_id);

	int ed_flg = iEdgesTbl[type][key_id];
	if (ed_flg == 0) return 0;

	for (int i(0);i<nedges[type];++i) {

		if (ed_flg & (1 << i)) {

			int v0 = iEdgesConnect[type][i][0];
			bg.position.x = vertices[3*v0];
			bg.position.y = vertices[3*v0+1];
			bg.position.z = vertices[3*v0+2];
			int v1 = iEdgesConnect[type][i][1];
			ed.position.x = vertices[3*v1];
			ed.position.y = vertices[3*v1+1];
			ed.position.z = vertices[3*v1+2];
			//std::cout << "Edge " << i << " indices in key: " << v0 << " " << v1 <<std::endl;

			// 0 <= u <= 1
			double u;
			//std::cout << "\nbegin: " << bg.position << " end: " << ed.position;
			if (cut_plane->IntersectWithLine(ed.position.v,bg.position.v,u)) {
				//std::cout << "\t intersection: " << ed.position << std::endl;
				//int type = VertexType::EDGE_SHAPE;
				if (fvmath::Compare(u,0.0)) {
					ed.info = nodes[v0+1];
					//type |= VertexType::CORNER_SHAPE;
				}
				else if (fvmath::Compare(u,1.0)) {
					ed.info = nodes[v1+1];
					//type |= VertexType::CORNER_SHAPE;
				}
				else ed.info = 0;

				// Remove duplicates ariesd from coincident edges
				bool found=false;
				unsigned iv=0;
				while (!found && iv < out_vertices.size()) {
					found = (ed.position == out_vertices[iv].position);
					iv++;
				}

				if (!found) out_vertices.push_back(ed);
			}
			//else std::cout << "\t no intersection\n";
		}
	}
	//std::cout <<"Out_lines: " << out_vertices.size() << std::endl;
	assert(out_vertices.size()>=3);

	// Sort vertices in CCW order
	const myVec3d plane_normal = myVec3d(cut_plane->GetParams());
	const myVec3d orig = out_vertices[0].position;

	std::sort(out_vertices.begin(),out_vertices.end(),
			[&](const myVertex& lhs,const myVertex& rhs) -> bool {
				myVec3d v = (lhs.position - orig) * (rhs.position - orig);
				return FemViewer::fvmath::Dot(v,plane_normal) > 0;
			} );

	return out_vertices.size();
}


typedef unsigned int indextype;

bool getSimilarVertexIndex(
	Vertex & packed,
	std::map<Vertex,indextype,comp> & VertexToOutIndex,
	indextype & result)
{
	std::map<Vertex,indextype,comp>::iterator it = VertexToOutIndex.find(packed);
	if ( it == VertexToOutIndex.end() ){
		return false;
	}else{
		result = it->second;
		return true;
	}
}




void indexVBO(std::vector< Vertex >& inout_vertices,
			  std::vector< indextype >& inout_indices)
{
	// Temporary container
	std::vector<Vertex> out_vertices;
	out_vertices.reserve(inout_vertices.size());

	assert(!inout_indices.empty());

	indextype index = 0;
	indextype num = inout_indices.size();
	std::map<Vertex,indextype,comp> VertexToOutIndex;

	// For each input index
	for (indextype i=0; i < num; i++)
	{
		// Get index value
		indextype vt_index = inout_indices[i];
		// Get vertex from input vector
		Vertex v = Vertex(inout_vertices[ vt_index ]);

		// Check whether vertex is present already
		indextype ii;
		bool found = getSimilarVertexIndex(v,VertexToOutIndex,ii);

		if (found) {
			inout_indices[i] = ii;
		}
		else {
			//mfp_debug("Size of vertices: %u",out_vertices.size());
			out_vertices.push_back(v);
			indextype newindex = out_vertices.size() - 1;
			inout_indices[i] = newindex;
			VertexToOutIndex[v] = newindex;
		}
	}
	// Resize
	//std::vector<Vertex>(out_vertices).swap(inout_vertices);
	inout_vertices = out_vertices;
}


template<typename T>
T mark_edges_in_triangle(int s,int u,const int dim)
{
  // Mark corners
  if ((u + s == 0) || (s == dim) || (u == dim)) return T(1);
  // Mark top edge
  else if (u == 0) return T(2);
  // Mark skew edgecust
  else if (u + s == dim) return(3);
  // Mark left edge
  else if (s == 0) return T(4);
  // Mark internal vertices
  return T(0);
}

// Calculate points within triangle by its corresponding
// barycentric coordinates calculated in OXY unit triangle
// return the number of new vertices
unsigned TessellateTriangle(
		const myVertex* a, /* Traingle corners */
		const myVertex* b,
		const myVertex* c,
		int ndiv,		/* The edensity of sampling */
		std::vector<myVertex>& vertices, /* out: buffer of vertices */
		std::vector<unsigned>& indices   /* coressponding indices */
		)
{
	// get a num of vertices as current offset
	unsigned base = vertices.size();
	unsigned size = base;

	//std::cout << "Before inserting: " << vertices.size() << std::endl;
	if (ndiv == 1) {
		vertices.push_back(*a);
		vertices.push_back(*b);
		vertices.push_back(*c);
		indices.push_back(base++);
		indices.push_back(base++);
		indices.push_back(base);
	}
	else {

		// Generate point coords in columns
		const double dlx = 1.0 / (double)ndiv;
        const double dly = dlx;

        // Upgrade ndiv to num of points
        ndiv++;

        myVertex node;
        int type, info;

        for (int x=0; x < ndiv; x++) {

			node.info = 0;

			if (x == 0) {
				type = c->type & VertexType::EDGE_SHAPE ? 1 : 0;
			}
		    else if (x == (ndiv-1)) {
		    	node = *b;
		    	node.type = VertexType::EDGE_SHAPE; //? skoro w b juz jest
		    	vertices.push_back(node); break;
		    }
			else {
				type = a->type & VertexType::EDGE_SHAPE ?  1 : 0;
			}

			for (int y=0, ny=ndiv-x; y < ny; y++) {
				int edge_flg = mark_edges_in_triangle<int>(y,x,ndiv-1);
				if (x == 0) {
					if (y == 0) {
						node = *a;
						node.type = VertexType::EDGE_SHAPE;
						vertices.push_back(node);
						continue;
					}
					else if (y == (ndiv-1)) {
						node = *c;
						node.type = VertexType::EDGE_SHAPE;
						vertices.push_back(node);
						break;
					}
					else {
						type = c->type & VertexType::EDGE_SHAPE ? 1 : 0;
					}
				}
				else if (y == (ndiv-1-x)) {
					type = b->type & VertexType::EDGE_SHAPE ? 1 : 0;
				}
				else if (y > 0) {
					type = 0;
				}


				double beta  = x*dlx;
				double gamma = y*dly;
				double alfa  = 1.0 - beta - gamma;
				assert(gamma >= 0 && gamma <= 1.0);
				node.position = alfa*a->position + beta*b->position + gamma*c->position;
				node.info = 0;
				node.type = type * edge_flg;
				vertices.push_back(node);
			}
		}

        // Generates indices
        uint_t offset = ndiv--;
        // Loop over small triangles within big triangle
        // The number of small triangles is ndiv * ndiv
        for (int x = 0; x < ndiv; x++) {
        	int ny = ndiv-x;
            //std::cout << "offset: " << offset << " ny: " << ny << "\n";
            for (int y = 0; y < ny; y++) {
            	// Odd triangle
            	unsigned id1 = base + y;
            	unsigned id3 = id1 + 1;
            	unsigned id2 = id1 + offset;
            	indices.push_back(id1);
            	indices.push_back(id2);
            	indices.push_back(id3);
            	// Even triangle
            	if ((y+1) < ny) {
            		indices.push_back(id3);
            		indices.push_back(id2);
            		indices.push_back(++id2);
            	}
            }
            offset--;
            base += ny+1;
        }
    }

//	std::cout << "VERTICES:\n";
//	for (const auto & v : vertices) std::cout << v.position << std::endl;
//	std::cout << "TRIANGLE INDICES:\n";
//	for (unsigned i=0 ; i< indices.size(); i+=3)
//		std::cout << "{"<<indices[i]<<", " << indices[i+1] << ", " << indices[i+2] << "}\n";
	return vertices.size() - size;
}



