#include "Fire.h"
#include "MathHelper.h"
#include "Point3D.h"
#include "Legend.h"
#include "GraphicElem.hpp"
#include "fv_assert.h"

#include<cassert>

namespace FemViewer {

const int Fire::indexCase[6][3] = {
		{0, 1, 2}, {1, 2, 0}, {2, 0, 1},
		{2, 1, 0}, {1, 0, 2}, {0, 2, 1}
};




	void Fire::Exe3(GraphElement2<double>& elIn, std::vector<GraphElement2<double> >& elOut)
	{
		std::vector<double> vDValues;
		Point3D<double> p;
		GraphElement2<double> triangle;
		vGEL2d subTri;

		elIn.SortVertex(*this);

		// clear vector
		elOut.clear();

		// get borders from min to max values
		legend.GetBorders(elIn.values[0],elIn.values[2], vDValues);
		
		// check
		if(vDValues.size() == 2) {
			// all vertices are in the same bound
			legend.GetColors(elIn);
			if(elIn.colors.size() != 3) {
			FV_ASSERT(elIn.colors.size() == 3);
			}

			//std::cout<< elIn.AsString();
			if(!elIn.Check()) {
				elOut.push_back(elIn);
			} else std::cout << "tn nie exe3\n" << triangle.AsString();
			return;
		}
	
		// extract pentagram
		int line1, line2;
		ExtractPentagram2(elIn,vDValues,elOut,subTri,line1,line2);
		
		int idx_sub(0);
		// cut sub-triangles - left and right
		if(line1 > -2) CutTriangle(subTri[idx_sub++],vDValues,line1,elOut); // left
		if(line2 >= 0) CutTriangle(subTri[idx_sub],vDValues,line2,elOut); // right
		

		FV_ASSERT(elOut.size() > 0);

		return;
	}

	void Fire::CreateLine(GraphElement2<double>& elIn, std::vector<GraphElement2<double> >& elOut)
	{
		double val, val0, val1, val2;

		Point3D<double> ps, pe, p0, p1, p2;
		GraphElement2<double> line;//(2);

		elIn.SortVertex(*this);
		
		val0 = elIn.values[0];
		val1 = elIn.values[1];
		val2 = elIn.values[2];

		p0 = elIn.points[0];
		p1 = elIn.points[1];
		p2 = elIn.points[2];



		// clear vector
		elOut.clear();

		for(register size_t i=1; i <= legend.GetColors().size() - 1; i++)
		{
			val  = legend.GetWartoscGranicy(i);

			if(fvmath::ValueWithin(val,val0,val1))
			{
				ps = GetContourPoint<double>(p0,p1,val,val0,val1);
				pe = GetContourPoint<double>(p0,p2,val,val0,val2);
				if(ps.Compare(pe)) continue;

			} else if(fvmath::ValueWithin(val,val1,val2))
			{
				ps = GetContourPoint<double>(p1,p2,val,val1,val2);
				pe = GetContourPoint<double>(p0,p2,val,val0,val2);
				if(ps.Compare(pe)) continue;

			} else {
				continue;
			}
			// store
			line.Clear();
			line.Add(ps,val);
			line.Add(pe,val);
			if(!line.Check2()) {
				legend.GetColors(line);
				elOut.push_back(line);
			}
		}
		
		return;

	}

	int Fire::GetSorted(const double *values,const int*& index)
	{
		int result = 0;

		for(int j=0;j<6;++j)
		{
			if(values[indexCase[j][0]]<= values[indexCase[j][1]] && values[indexCase[j][1]] <= values[indexCase[j][2]]) {
				index = indexCase[j];
				result = (j>2) ? -1 : 1;
				break;
			}
		}
		return result;
	}

//	void Fire::ExtractPentagram( const GraphElement& el,
//								const std::vector<double>& values,
//								const int *index,
//								vGEL2& elOut,
//								vGEL2& subTriangles,
//								int& line1, int& line2)
//	{
//		double interp;
//		//const mmt_nodes *v0, *v1, *v2;
//		double Node0[3], Node1[3], Node2[3];
//		GraphElement2<double> triangle;
//		Point3D<double> p12, p13, p23, p31, p1, p2, tmp;
//
//		// clear vectors
//		elOut.clear();
//		subTriangles.clear();
//
//		for(size_t i(0); i < (values.size()-1); ++i)
//		{
//			if((values[i] >= el.values[index[0]] && values[i] < el.values[index[1]])
//				&&
//				(values[i+1] > el.values[index[1]] && values[i+1] <= el.values[index[2]]))
//			{
//				// get vertices from global table nodes
//				//v0 = pMesh->GetNode( el.vertex[index[0]]);
//				//v1 = pMesh->GetNode( el.vertex[index[1]]);
//				//v2 = pMesh->GetNode( el.vertex[index[2]]);
//
//				pMesh->Get_node_coor(pMesh->Id,el.vertex[index[0]],Node0);
//				pMesh->Get_node_coor(pMesh->Id,el.vertex[index[1]],Node1);
//				pMesh->Get_node_coor(pMesh->Id,el.vertex[index[2]],Node2);
//
//				// create points
//				interp = legend.GetInterp(el.values[index[0]],el.values[index[1]],values[i]);
//				//p1.x = v0->x; p1.y = v0->y; p1.z = v0->z;
//				//p2.x = v1->x; p2.y = v1->y; p2.z = v1->z;
//				p1.x = Node0[0]; p1.y = Node0[1]; p1.z = Node0[2];
//				p2.x = Node1[0]; p2.y = Node1[1]; p2.z = Node1[2];
//				tmp = p2;
//				p12.SetInterp(p1,p2,interp);
//
//				interp = legend.GetInterp(el.values[index[0]],el.values[index[2]],values[i]);
//				//p2.x = v2->x; p2.y = v2->y; p2.z = v2->z;
//				p2.x = Node2[0]; p2.y = Node2[1]; p2.z = Node2[2];
//				p13.SetInterp(p1,p2,interp);
//
//				interp = legend.GetInterp(el.values[index[1]],el.values[index[2]],values[i+1]);
//				p23.SetInterp(tmp,p2,interp);
//
//				interp = legend.GetInterp(el.values[index[0]],el.values[index[2]],values[i+1]);
//				p31.SetInterp(p1,p2,interp);
//
//				// clear elemnt
//				triangle.Clear();
//				// add 1 triangle p12,p23,tmp
//				triangle.Add(p12,values[i]);
//				triangle.Add(p23,values[i+1]);
//				triangle.Add(tmp,el.values[index[1]]);
//				legend.GetColors(triangle);
//				FV_ASSERT(triangle.colors.size() == 3); ///Tutaj sie najczesciej wywala ale czemu!!!!!
//				if(!triangle.Check()) elOut.push_back(triangle);
//
//				// clear elemnt
//				triangle.Clear();
//				// add 2 triangle p13,p23,p12
//				triangle.Add(p13,values[i]);
//				triangle.Add(p23,values[i+1]);
//				triangle.Add(p12,values[i]);
//				legend.GetColors(triangle);
//				FV_ASSERT(triangle.colors.size() == 3);
//				if(!triangle.Check()) {
//					elOut.push_back(triangle);
//				} /*else std::cout << "tn nie" << triangle.AsString();*/
//
//
//				// clear elemnt
//				triangle.Clear();
//				// add 3 triangle p13,p31,p23
//				triangle.Add(p13,values[i]);
//				triangle.Add(p31,values[i+1]);
//				triangle.Add(p23,values[i+1]);
//				legend.GetColors(triangle);
//				FV_ASSERT(triangle.colors.size() == 3);
//				if(!triangle.Check()) elOut.push_back(triangle);
//
//				//create subtriangles
//
//				triangle.Clear();
//				// add left sub-traingle
//				triangle.Add(p1, el.values[index[0]]);
//				triangle.Add(p13,values[i]);
//				triangle.Add(p12,values[i]);
//				// check if subtriangle is not a point or a line
//				if(!triangle.Check()) {
//					triangle.SortVertex(*this);
//					subTriangles.push_back(triangle);
//					line1 = i;
//				} else line1 = -2;
//
//				triangle.Clear();
//				// add right sub-triangle
//				triangle.Add(p23,values[i+1]);
//				triangle.Add(p31,values[i+1]);
//				triangle.Add(p2,el.values[index[2]]);
//				// check if subtriangle is not a point or a line
//				if(!triangle.Check()) {
//					triangle.SortVertex(*this);
//					subTriangles.push_back(triangle);
//					line2 = i + 1;
//				} else line2 = -2;
//
//				//line1 = i;
//				//line2 = i + 1;
//
//				return;
//
//
//			} // if
//
//		}// for over values
//
//		// convert GraphElement 2 GraphElement2
//		el.ToGraphElement2(*pMesh,triangle);
//		triangle.SortVertex(*this);
//		subTriangles.push_back(triangle);
//
//		line1 = line2 = -1;
//
//		return;
//
//	}

	void Fire::ExtractPentagram2(GraphElement2<double>& el, const std::vector<double>& values, vGEL2d& elOut, vGEL2d& subTriangles, int& line1, int& line2)
	{
		double interp;
		Point3D<double> p12, p13, p23, p31, p1, p2, tmp;
		GraphElement2<double> triangle;
		
		// clear vectors
		elOut.clear();
		subTriangles.clear();

		for(size_t i(0); i < (values.size()-1); ++i)
		{
			if((values[i] >= el.values[0] && values[i] < el.values[1])
				&&
				(values[i+1] >= el.values[1] && values[i+1] <= el.values[2]))
			{

				p1 = el.points[0];
				tmp = p2 = el.points[1];
	
				// create points
				interp = legend.GetInterp(el.values[0],el.values[1],values[i]);
				tmp = p2;
				p12.SetInterp(p1,p2,interp);
				
				interp = legend.GetInterp(el.values[0],el.values[2],values[i]);
				p2 = el.points[2];
				p13.SetInterp(p1,p2,interp);

				interp = legend.GetInterp(el.values[1],el.values[2],values[i+1]);
				p23.SetInterp(tmp,p2,interp);
	
				interp = legend.GetInterp(el.values[0],el.values[2],values[i+1]);
				p31.SetInterp(p1,p2,interp);

				// clear elemnt
				triangle.Clear();
				// add 1 triangle p12,p23,tmp
				triangle.Add(p12,values[i]);
				triangle.Add(p23,values[i+1]);
				triangle.Add(tmp,el.values[1]);
				legend.GetColors(triangle);
				FV_ASSERT(triangle.colors.size() == 3); ///Tutaj sie najczesciej wywala ale czemu!!!!!
				//std::cout<< triangle.AsString();
				if(!triangle.Check()) {
					elOut.push_back(triangle);
				} /*else std::cout << "tn nie" << triangle.AsString();*/

				// clear elemnt
				triangle.Clear();
				// add 2 triangle p13,p23,p12
				triangle.Add(p13,values[i]);
				triangle.Add(p23,values[i+1]);
				triangle.Add(p12,values[i]);
				legend.GetColors(triangle);
				FV_ASSERT(triangle.colors.size() == 3);
				//std::cout<< triangle.AsString();
				if(!triangle.Check()) {
					elOut.push_back(triangle);
				} /*else std::cout << "tn nie" << triangle.AsString();*/

				// clear elemnt
				triangle.Clear();
				// add 3 triangle p13,p31,p23
				triangle.Add(p13,values[i]);
				triangle.Add(p31,values[i+1]);
				triangle.Add(p23,values[i+1]);
				legend.GetColors(triangle);
				FV_ASSERT(triangle.colors.size() == 3);
				//std::cout<< triangle.AsString();
				if(!triangle.Check()) {
					elOut.push_back(triangle);
				} /*else std::cout << "tn nie\n" << triangle.AsString();*/

				//create subtriangles

				triangle.Clear();
				// add left sub-traingle
				triangle.Add(p1, el.values[0]);
				triangle.Add(p13,values[i]);
				triangle.Add(p12,values[i]);
				// check if subtriangle is not a point or a line
				if(!triangle.Check()) {
					triangle.SortVertex(*this);
					subTriangles.push_back(triangle);
					line1 = i;
				} else line1 = -2;

				triangle.Clear();
				// add right sub-triangle
				triangle.Add(p23,values[i+1]);
				triangle.Add(p31,values[i+1]);
				triangle.Add(p2,el.values[2]);
				// check if subtriangle is not a point or a line
				if(!triangle.Check()) {
					triangle.SortVertex(*this);
					subTriangles.push_back(triangle);
					line2 = i + 1;
				} else line2 = -2;

				return;


			} // if

		}// for over values

		el.SortVertex(*this);
		subTriangles.push_back(el);

		line1 = line2 = -1;

		return;
	}

	void Fire::CutTriangle(GraphElement2<double>& elIn, const std::vector<double>& values, int line, vGEL2d& elOut)
	{
		double interp;
		Point3D<double> p23, p13;
		GraphElement2<double> triangle;

		if(FullTriangleFallage(elIn,values,elOut))
			return;

		for(register int i = 0; i < static_cast<int>(values.size()); ++i)
		{
			if(i == line) continue;

			if(values[i] >= elIn.values[1] && values[i] <= elIn.values[2])
			{
				interp = legend.GetInterp(elIn.values[1], elIn.values[2], values[i]);
				p23.SetInterp(elIn.points[1],elIn.points[2],interp);

				interp = legend.GetInterp(elIn.values[0], elIn.values[2], values[i]);
				p13.SetInterp(elIn.points[0],elIn.points[2],interp);

				triangle.Clear();
				triangle.Add(elIn.points[0],elIn.values[0]);
				triangle.Add(p13,values[i]);
				triangle.Add(p23,values[i]);
				legend.GetColors(triangle);
				if ( triangle.colors.size() != 3) {
				FV_ASSERT(triangle.colors.size() == 3);
				}
				if(!triangle.Check()) {
					elOut.push_back(triangle);
				} /*else std::cout << "tn nie\n" << triangle.AsString();*/

				triangle.Clear();
				triangle.Add(elIn.points[0],elIn.values[0]);
				triangle.Add(p23,values[i]);
				triangle.Add(elIn.points[1],elIn.values[1]);
				legend.GetColors(triangle);
				FV_ASSERT(triangle.colors.size() == 3);
				if(!triangle.Check()) {
					elOut.push_back(triangle);
				} /*else std::cout << "tn nie\n" << triangle.AsString();*/

				triangle.Clear();
				triangle.Add(p23,values[i]);
				triangle.Add(p13,values[i]);
				triangle.Add(elIn.points[2],elIn.values[2]);
				//FV_ASSERT(!triangle.Check());
				triangle.SortVertex(*this);

				CutTriangle(triangle,values,i,elOut);
				return;

			}// if for right 

			if((values[i] >= elIn.values[0]) && (values[i] <= elIn.values[1])
			&& (i == values.size() -1 || /*values[i+1] >= elIn.values[1] ||*/ values[i+1] >= elIn.values[2])) 
			{
				interp = legend.GetInterp(elIn.values[0], elIn.values[2], values[i]);
				p23.SetInterp(elIn.points[0],elIn.points[2],interp);

				interp = legend.GetInterp(elIn.values[0], elIn.values[1], values[i]);
				p13.SetInterp(elIn.points[0],elIn.points[1],interp);

				triangle.Clear();
				triangle.Add(elIn.points[2],elIn.values[2]);
				triangle.Add(p23,values[i]);
				triangle.Add(p13,values[i]);
				//FV_ASSERT(!triangle.Check());
				legend.GetColors(triangle);
				FV_ASSERT(triangle.colors.size() == 3);
				if(!triangle.Check()) {
					elOut.push_back(triangle);
				} /*else std::cout << "tn nie\n" << triangle.AsString();*/

				triangle.Clear();
				triangle.Add(p13,values[i]);
				triangle.Add(elIn.points[1],elIn.values[1]);
				triangle.Add(elIn.points[2],elIn.values[2]);
				//FV_ASSERT(!triangle.Check());
				legend.GetColors(triangle);
				FV_ASSERT(triangle.colors.size() == 3);
				if(!triangle.Check()) {
					elOut.push_back(triangle);
				} /*else std::cout << "tn nie\n" << triangle.AsString();*/

				triangle.Clear();
				triangle.Add(elIn.points[0],elIn.values[0]);
				triangle.Add(p13,values[i]);
				triangle.Add(p23,values[i]);
				//FV_ASSERT(!triangle.Check());
				triangle.SortVertex(*this);			

				CutTriangle(triangle,values,i,elOut);
				return;

			}// if for lef
		}// for
	}

	bool Fire::FullTriangleFallage(GraphElement2<double>& elIn, const std::vector<double>& values, vGEL2d& elOut)
	{

		for(int i=0; i<values.size(); ++i)
		{
			if(values[i] >= elIn.values[2] && i == 0)
			{
				legend.GetColors(elIn);
				FV_ASSERT(elIn.colors.size() == 3);
				elOut.push_back(elIn);
				return true;
			} 
			else if(values[i] <= elIn.values[0] && i < (values.size() -1) && values[i+1] >= elIn.values[2])
			{
				legend.GetColors(elIn);
				if(elIn.colors.size() != 3) {
					FV_ASSERT(elIn.colors.size() == 3);
				}
				elOut.push_back(elIn);
				return true;
			}
			else if(values[i] <= elIn.values[0] && i == (values.size() -1))
			{
				legend.GetColors(elIn);
				FV_ASSERT(elIn.colors.size() == 3);
				elOut.push_back(elIn);
				return true;
			}
		}// for
		return false;
	}






} //end namespace FEMSv 
