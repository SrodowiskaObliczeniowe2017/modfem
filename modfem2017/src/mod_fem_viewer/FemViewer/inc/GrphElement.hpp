#ifndef _GRAPH_ELEMENT_H_
#define _GRAPH_ELEMENT_H_

#include<string.h>

namespace FemViewer
{


	class GraphElement
	{
		public:

			int type;//dane typu figury
			int vertex[8];//indeksy wsp�ednych punkt�w
			double values[8];//warto�ci w punktach


			//float norm[3];
			//wektor normalny do figury

			void clear(void)
			{
				memset(vertex,0,sizeof(int)   << 3);
				memset(values,0,sizeof(double) << 3);
			};


			double GetInterpValue(double a, double b, double interp){return a + (b-a) * interp;};

	};
} // end namespace FemViewer

#endif /* _GRAPH_ELEMENT_H_
*/
