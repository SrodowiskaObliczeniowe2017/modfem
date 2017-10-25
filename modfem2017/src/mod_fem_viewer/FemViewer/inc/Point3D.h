#ifndef _POINT3D_H_
#define _POINT3D_H_

#include "MathHelper.h"
#include <sstream>
#include <string>

namespace FemViewer {

	template <class T>
	class Point3D
	{
		public:
			T x, y, z;


		public:
			//! Konstruktor domy�lny wyplenia wartosciami domyslnymi dla danego typu dnaych
			Point3D(void) : x(T()), y(T()), z(T()){};

			//! Konstruktor inicjuj�cy trzema warto�ciami.
			Point3D(T a, T b, T c):x(a), y(b), z(c){};

			//! Konstruktor kopiuj�cy zawarto�c z istniej�cego ju� obiektu
			Point3D(const Point3D<T> & src) : x(src.x), y(src.y), z(src.z){};

			//! Operator przypisania
			Point3D<T>& operator=(const Point3D<T>& rhs)
			{
				x = rhs.x; y = rhs.y; z = rhs.z;
				return *this;
			}

			void Set(T a, T b, T c)
			{
				x = a;
				y = b;
				z = c;
			};

			//! Metoda por�wnyj�ca warto�ci punktu z wartosciami podanymi jako parametr
			/*!
				Zwraca tru jesli sa zgodne, false je�li si� r�znia
			*/
			bool Compare(const Point3D<T> & p) const
			{
				//if ( (x == p.x) && (y == p.y) && (z == p.z) )
				if
				(
				( fvmath::Compare(x, p.x) )
					&&
					( fvmath::Compare(y, p.y) )
					&&
					( fvmath::Compare(z, p.z) )
				)
				{
					return true;
				}
				else
					return false;
			};


			void SetInterp(Point3D<T> a, Point3D<T> b, double interp)
			{

				//lv.red = legendValues[indexMin].red + (legendValues[indexMax].red - legendValues[indexMin].red) * interp;
				this->x = a.x + (b.x-a.x) * interp;
				this->y = a.y + (b.y-a.y) * interp;
				this->z = a.z + (b.z-a.z) * interp;

			};

			//! Zwraca wartosc <0-1> interpolacji punktu p0 le��cego miedzy p1 a p2
			static double GetInterp(Point3D<T> p1, Point3D<T> p2, Point3D<T> p0)
			{
				double l, la;

				l = p1.Length(p2);

				la = p1.Length(p0);

				return la/l;
			};


			double Length(const Point3D<T> & p)
			{
				T vx, vy, vz;
				vx = p.x - x;
				vy = p.y - y;
				vz = p.z - z;

				return sqrtf(  vx*vx + vy*vy + vz*vz  );
			};

			std::string AsString(void) const
			{
				std::stringstream oss;

				oss << x << ", " << y << ", " << z;

				return oss.str();
			}
	};


} // end namespace FEMSv

#endif /* _POINT3D_H_
*/
