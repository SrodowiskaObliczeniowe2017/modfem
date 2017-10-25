#ifndef _TILE_H
#define _TILE_H

namespace FemViewer {

	class Tile
	{
	public:
		/* Constructors */
		Tile(const int iImgWidth,
			 const int iImgHeight,
			 const float fAspectR);

		Tile(const int iImgWidth,
			 const int iImgHeight,
		     const int iXmin,
		     const int iYmin,
		     const int iXmax,
		     const int iYmax);
		
		/* Destructor 
		**/
		~Tile();

		/* Initialize frustrum matrix 
		**/
		void InitFrustrumMatrix(const float fFOV,
							    const float fZNear,
							    const float fZFar,
							    float mout[]) const;

		/* Initialize frustrum matrix from array of floats
		**/
		//void InitFrustrumMatrix(const float mt[]) const;

		/* Initialize viewport
		**/
		void InitOrtho2DMatrix() const;

		/* Initialize ortho matrix
		**/
		void InitiOrtho3DMatrix(const float fFOV, const float fZNear,
			const float fZFar, float mout[]) const;

		/* Initialize viewport
		**/
		void InitViewport() const;

	private:


	#ifdef FV_DEBUG
		void invariants() const;
	#endif
		/* Private members */
		float fImageWidth;
		float fImageHeight;
		float fAspectRatio;
		float fXmin;
		float fXmax;
		float fYmin;
		float fYmax;
	};

} // end namespace FemViewer

#endif // TILE_H
