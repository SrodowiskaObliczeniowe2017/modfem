#include "Matrix.h"

namespace FemViewer {
template<>
bool Matrix<float>::IsSingular() const
{
	float det = 
		matrix[0]*(matrix[5]*matrix[10] - matrix[6]*matrix[9]) +
		matrix[1]*(matrix[6]*matrix[8] - matrix[4]*matrix[10]) +
		matrix[2]*(matrix[4]*matrix[9] - matrix[5]*matrix[8]);
	return(det == 0.0f);
}

template<>
void Matrix<float>::rotateAbout(const float* coord, const float angle)
{
	Vec3D vec(coord);
	const float co = cos(angle);
	const float si = sin(angle);
	const float vxx = vec._x() * vec._x();
	const float vyy = vec._y() * vec._y();
	const float vzz = vec._z() * vec._z();

	const float lRot00  = vxx+co*(1-vxx);
    const float lRot01  = vec._x()*vec._y()*(1-co)-vec._z()*si;
    const float lRot02  = vec._z()*vec._x()*(1-co)+vec._y()*si;
    const float lRot10  = vec._x()*vec._y()*(1-co)+vec._z()*si;
    const float lRot11  = vyy+co*(1-vyy);
    const float lRot12  = vec._y()*vec._z()*(1-co)-vec._x()*si;
    const float lRot20  = vec._z()*vec._x()*(1-co)-vec._y()*si;
    const float lRot21  = vec._y()*vec._z()*(1-co)+vec._x()*si;
    const float lRot22  = vzz+co*(1-vzz);

    const float lVals00 = matrix[0];
    const float lVals01 = matrix[4];
    const float lVals02 = matrix[8];
    const float lVals10 = matrix[1];
    const float lVals11 = matrix[5];
    const float lVals12 = matrix[9];
    const float lVals20 = matrix[2];
    const float lVals21 = matrix[6];
    const float lVals22 = matrix[10];

    matrix[0] = lVals00*lRot00 + lVals01*lRot10 + lVals02*lRot20;
    matrix[4] = lVals00*lRot01 + lVals01*lRot11 + lVals02*lRot21;
    matrix[8] = lVals00*lRot02 + lVals01*lRot12 + lVals02*lRot22;
    matrix[1] = lVals10*lRot00 + lVals11*lRot10 + lVals12*lRot20;
    matrix[5] = lVals10*lRot01 + lVals11*lRot11 + lVals12*lRot21;
    matrix[9] = lVals10*lRot02 + lVals11*lRot12 + lVals12*lRot22;
    matrix[2] = lVals20*lRot00 + lVals21*lRot10 + lVals22*lRot20;
    matrix[6] = lVals20*lRot01 + lVals21*lRot11 + lVals22*lRot21;
    matrix[10]= lVals20*lRot02 + lVals21*lRot12 + lVals22*lRot22;


}

template<>
void Matrix<float>::scale(const Vec3D& vec)
{
  const float x = vec._x();
  const float y = vec._y();
  const float z = vec._z();

  matrix[0] *= x; matrix[1] *= x; matrix[2] *= x;
  matrix[4] *= y; matrix[5] *= y; matrix[6] *= y;
  matrix[8] *= z; matrix[9] *= z; matrix[10] *= z;
}

template<>
void Matrix<float>::translate(const Vec3D& vec)
{
	const float x = vec._x();
	const float y = vec._y();
	const float z = vec._z();

	matrix[12] += matrix[0]*x + matrix[4]*y + matrix[8]*z;
	matrix[13] += matrix[1]*x + matrix[5]*y + matrix[9]*z;
	matrix[14] += matrix[2]*x + matrix[6]*y + matrix[10]*z; 
}

template<>
Matrix<float> Matrix<float>::Invers() const
{
	if (this->isIdentity) return *this;

	int i, j, k;
	Matrix s, t (*this);

	// Forward elimination
	for (i = 0; i < 3 ; i++) {
		int pivot = i;
		float pivotsize = t[4*i + i];

       	if (pivotsize < 0)
       		pivotsize = -pivotsize;

       	for (j = i + 1; j < 4; j++) {
       		float tmp = t[4*i + j];

			if (tmp < 0)
				tmp = -tmp;

			if (tmp > pivotsize) {
				pivot = j;
				pivotsize = tmp;
			}
		}

		if (pivotsize == 0) return Matrix();


		if (pivot != i) {
			for (j = 0; j < 4; j++) {
				float tmp;

				tmp = t[4*j + i];
				t[4*j + i] = t[4*j + pivot];
				t[4*j + pivot] = tmp;

				tmp = s[4*j + i];
				s[4*j + i] = s[4*j + pivot];
				s[4*j + pivot] = tmp;
			}
		}

		for (j = i + 1; j < 4; j++) {
			float f = t[4*j + i] / t[4*i + i];
			for (k = 0; k < 4; k++) {
				t[4*k + j] -= f * t[4*k + i];
				s[4*k + j] -= f * s[4*k + i];
			}
		}
	}

	// Backward substitution
	for (i = 3; i >= 0; --i) {
		float f;
		if ((f = t[4*i + i]) == 0) {
			return Matrix();
		}

		for (j = 0; j < 4; j++) {
			t[4*j + i] /= f;
			s[4*j + i] /= f;
		}
		for (j = 0; j < i; j++) {
			f = t[4*i + j];

			for (k = 0; k < 4; k++) {
				t[4*k + j] -= f * t[4*k + i];
				s[4*k + j] -= f * s[4*k + i];
			}
		}
	}

	return s;
}

} // end namespace FemViewer
