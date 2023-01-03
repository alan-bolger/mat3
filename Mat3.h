#ifndef MAT3_H
#define MAT3_H

#include "Vec3.h"

constexpr auto PI = 3.14159265358979323846;

class Mat3
{
public:
	double m11;
	double m21;
	double m31;
	double m12;
	double m22;
	double m32;
	double m13;	
	double m23;	
	double m33;

	Mat3()
	{ 
		m11 = 0.0;
		m21 = 0.0;
		m31 = 0.0;
		m12 = 0.0;
		m22 = 0.0;
		m32 = 0.0;
		m13 = 0.0;
		m23 = 0.0;
		m33 = 0.0;
	};

	Mat3(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33)
	{ 
		m11 = static_cast<double>(a11);
		m21 = static_cast<double>(a21);
		m31 = static_cast<double>(a31);
		m12 = static_cast<double>(a12);
		m22 = static_cast<double>(a22);
		m32 = static_cast<double>(a32);
		m13 = static_cast<double>(a13);
		m23 = static_cast<double>(a23);
		m33 = static_cast<double>(a33);
	};

	Mat3(Vec3 r1, Vec3 r2, Vec3 r3) 
	{
		m11 = r1.x;
		m21 = r2.x;
		m31 = r3.x;
		m12 = r1.y;
		m22 = r2.y;
		m32 = r3.y;
		m13 = r1.z;
		m23 = r2.z;
		m33 = r3.z;
	};

	~Mat3()	{ };

	Mat3 transpose() const
	{
		return Mat3(m11, m21, m31, m12, m22, m32, m13, m23, m33);
	};

	double determinant() const
	{
		return m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m31 * m23) + m13 * (m21 * m32 - m22 * m31);
	};

	Mat3 inverse() const
	{
		double det = determinant();

		if (det != 0)
		{
			Mat3 trans = transpose();

			// 1. Get matrix of cofactors from transpose of original matrix
			// 2. Get adjugate matrix by changing the sign of every second element
			// 3. Divide each term by the determinant
			double a11 = (trans.m22 * trans.m33 - trans.m23 * trans.m32) * 1 / det;
			double a12 = ((trans.m21 * trans.m33 - trans.m23 * trans.m31) * -1) * 1 / det;
			double a13 = (trans.m21 * trans.m32 - trans.m22 * trans.m31) * 1 / det;
			double a21 = ((trans.m12 * trans.m33 - trans.m13 * trans.m32) * -1) * 1 / det;
			double a22 = (trans.m11 * trans.m33 - trans.m13 * trans.m31) * 1 / det;
			double a23 = ((trans.m11 * trans.m32 - trans.m12 * trans.m31) * -1) * 1 / det;
			double a31 = (trans.m12 * trans.m23 - trans.m13 * trans.m22) * 1 / det;
			double a32 = ((trans.m11 * trans.m23 - trans.m13 * trans.m21) * -1) * 1 / det;
			double a33 = (trans.m11 * trans.m22 - trans.m12 * trans.m21) * 1 / det;

			return Mat3(a11, a12, a13, a21, a22, a23, a31, a32, a33);
		}
		else
		{
			// If there's no inverse, then return the original matrix
			return Mat3(m11, m12, m13, m21, m22, m23, m31, m32, m33);
		}
	};

	Vec3 row(const int row) const
	{
		if (row >= 0 && row < 3)
		{
			if (row == 0)
			{
				return Vec3(m11, m12, m13);
			}
			else if (row == 1)
			{
				return Vec3(m21, m22, m23);
			}
			else
			{
				return Vec3(m31, m32, m33);
			}
		}
		else
		{
			return Vec3(0.0, 0.0, 0.0);
		}
	};

	Vec3 col(const int col) const
	{
		if (col >= 0 && col < 3)
		{
			if (col == 0)
			{
				return Vec3(m11, m21, m31);
			}
			else if (col == 1)
			{
				return Vec3(m12, m22, m32);
			}
			else
			{
				return Vec3(m13, m23, m33);
			}
		}
		else
		{
			return Vec3(0.0, 0.0, 0.0);
		}
	};

	static Mat3 rotationX(const double radians)
	{
		return Mat3(1.0, 0.0, 0.0, 0.0, std::cos(radians), std::sin(radians * -1.0), 0.0, std::sin(radians), std::cos(radians));
	};

	static Mat3 rotationY(const double radians)
	{
		return Mat3(std::cos(radians), 0.0, std::sin(radians), 0.0, 1.0, 0.0, std::sin(radians * -1.0), 0.0, std::cos(radians));
	};

	static Mat3 rotationZ(const double radians)
	{
		return Mat3(std::cos(radians), std::sin(radians * -1.0), 0.0, std::sin(radians), std::cos(radians), 0.0, 0.0, 0.0, 1.0);
	};

	static Mat3 translate(const Vec3 displacement)
	{
		return Mat3(1.0, 0.0, displacement.x, 0.0, 1.0, displacement.y, 0.0, 0.0, 1.0);
	};

	static Mat3 scale(const double factor)
	{
		return Mat3(factor, 0.0, 0.0, 0.0, factor, 0.0, 0.0, 0.0, 1.0);
	};

	std::string toString() const
	{
		return "[" + std::to_string(m11) + "," + std::to_string(m12) + "," + std::to_string(m13) +
			"]\n[" + std::to_string(m21) + "," + std::to_string(m22) + "," + std::to_string(m23) +
			"]\n[" + std::to_string(m31) + "," + std::to_string(m32) + "," + std::to_string(m33) + "]";
	};

	bool operator== (const Mat3 mat) const
	{
		return m11 == mat.m11 && m12 == mat.m12 && m13 == mat.m13 && m21 == mat.m21 && m22 == mat.m22 && m23 == mat.m23 && m31 == mat.m31 && m32 == mat.m32 && m33 == mat.m33;
	};

	bool operator!= (const Mat3 mat) const
	{
		return m11 != mat.m11 || m12 != mat.m12 || m13 != mat.m13 || m21 != mat.m21 || m22 != mat.m22 || m23 != mat.m23 || m31 != mat.m31 || m32 != mat.m32 || m33 != mat.m33;
	};

	Mat3 operator+ (const Mat3 mat) const
	{
		return Mat3(m11 + mat.m11, m12 + mat.m12, m13 + mat.m13, m21 + mat.m21, m22 + mat.m22, m23 + mat.m23, m31 + mat.m31, m32 + mat.m32, m33 + mat.m33);
	};

	Mat3 operator- (const Mat3 mat) const
	{
		return Mat3(m11 - mat.m11, m12 - mat.m12, m13 - mat.m13, m21 - mat.m21, m22 - mat.m22, m23 - mat.m23, m31 - mat.m31, m32 - mat.m32, m33 - mat.m33);
	};

	Mat3 operator* (const Mat3 mat) const
	{
		double a11 = m11 * mat.m11 + m12 * mat.m21 + m13 * mat.m31;
		double a12 = m11 * mat.m12 + m12 * mat.m22 + m13 * mat.m32;
		double a13 = m11 * mat.m13 + m12 * mat.m23 + m13 * mat.m33;
		double a21 = m21 * mat.m11 + m22 * mat.m21 + m23 * mat.m31;
		double a22 = m21 * mat.m12 + m22 * mat.m22 + m23 * mat.m32;
		double a23 = m21 * mat.m13 + m22 * mat.m23 + m23 * mat.m33;
		double a31 = m31 * mat.m11 + m32 * mat.m21 + m33 * mat.m31;
		double a32 = m31 * mat.m12 + m32 * mat.m22 + m33 * mat.m32;
		double a33 = m31 * mat.m13 + m32 * mat.m23 + m33 * mat.m33;

		return Mat3(a11, a12, a13, a21, a22, a23, a31, a32, a33);
	};

	Vec3 operator* (const Vec3 vec) const
	{
		return Vec3(m11 * vec.x + m12 * vec.y + m13 * vec.z, m21 * vec.x + m22 * vec.y + m23 * vec.z, m31 * vec.x + m32 * vec.y + m33 * vec.z);
	};

	Mat3 operator* (const double scalar) const
	{
		return Mat3(m11 * scalar, m12 * scalar, m13 * scalar, m21 * scalar, m22 * scalar, m23 * scalar, m31 * scalar, m32 * scalar, m33 * scalar);
	};
};

#endif //!MAT3_H