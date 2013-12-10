// SWD_Maths.h: interface for the SWD_Maths class.
//
/*--------------------------------------------------------------
*
*\brief   about mesh computation on maths
*\author   李欧州
*\date     20-9-2008
*
*
*
------------------------------------------------------------------*/
//history
/*-----------------------------------------------------------------
*brief  modify for 3DQModelRetrieval
*author 李欧州
*date   17-02-2009 
-----------------------------------------------------------------*/
//////////////////////////////////////////////////////////////////////
#ifndef SWD_MATHS_H
#define SWD_MATHS_H
#include"myconfig.h"
class Vertex3d;
class Face;
class Mesh;
class SpheCoord;

class SWD_Maths  
{
public:
	static bool low_compare1(const SpheCoord & m1, const SpheCoord & m2);
	static bool up_compare1(const SpheCoord & m2,const SpheCoord & m1);
	static bool low_compare2(const SpheCoord & m1, const SpheCoord & m2);
	static bool up_compare2(const SpheCoord & m2,const SpheCoord & m1);

	static double AreaOneTriangle(Vertex3d v1,Vertex3d v2,Vertex3d v3);
	static Vertex3d* CenGravOneTriangle(Vertex3d v1, Vertex3d v2, Vertex3d v3);
	//计算模型的重心，和存储每个三角面片的面积和重心，和总面积
    static void ComputBarMesh(double* areaT, Vertex3d * veGArray,Mesh* mesh,Vertex3d* barMesh
				   , double& sAllTrian);
    //计算CovMat
    static void ComputCovMat(double* areaT, Vertex3d * veGArray,Mesh* mesh, double CovMat[3][3]
	 ,Vertex3d* barMesh,double& sAllTrian);
	
    static void ConvertCartesianToSpherical(Vertex3d& Pos, float& rLong, float& rLat);
	static bool less_second(const SpheCoord & m1, const SpheCoord & m2);
    static float max3(float a, float b, float c);
    static float min3(float a, float b, float c);
    static  void SortEigen(double eigenvalues[3], double eigenvectors[3][3]);
    static void ComputFaceNormal(Face& face);
    static void Jacobi_Cyclic_Method( double eigenvalues[], 
		double *eigenvectors, double *A, int n);

    static void ComputFirstMoments(double& fVector, float matrix[ImHeight][ImWidth]
		, int hFlag, int wFlag, int nLevel);
	static void ComputSecondMoments(double& fVector, double fVectorfirst, float matrix[ImHeight][ImWidth]
		, int hFlag, int wFlag, int nLevel);
	static void ComputThirdMoments(double& fVector, double fVectorfirst, float matrix[ImHeight][ImWidth]
		, int hFlag, int wFlag, int nLevel);

	static void ComputSWTFirstMoments(vector<double>& vect, vector<double>& descriptorVect);
    static void ComputSWTSecondMoments(vector<double>& vect, vector<double>& descriptorVect);
	static void ComputSWTThirdMoments(vector<double>& vect, vector<double>& descriptorVect);


	
private:
	static void ComputHelperCovMat(Vertex3d vert, double arrayF[3][3]);
	static void ComputHelperMat(double CovMat[3][3], double arrayF[3][3]);
	
    
};

#endif 
