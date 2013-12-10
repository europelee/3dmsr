// SWD_PoseNormalization.h: interface for the SWD_PoseNormalization class.
//
/*--------------------------------------------------------------
*
*\brief   3D model pre-process, pose normalization
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
#ifndef SWD_POSENORMALIZATION_H
#define SWD_POSENORMALIZATION_H

class Mesh;
class Vertex3d;

class SWD_PoseNormalization  
{
public:
	void ComputFlipInvia(int F[3][3], Mesh* mesh,double* areaT,double& sAllTrian);
	void ComputFlipInviaCPCA(int F[3][3], Mesh* mesh,double* areaT,double& sAllTrian);
	void ComputRotaInvia(Mesh* mesh, double eigenvectors[3][3]);
	void ComputScaleInvia(double& scaleV, Mesh*  mesh,double* areaT,double&  sAllTrian, Vertex3d *veGArray);
    void ComputTranslatInvia(Mesh* mesh,Vertex3d* barMesh);
    void PoseNormalization(Mesh* mesh,double&  scaleV, int F[3][3]);
    void ComputHelperFlip(double& fx, float x[3], double* areaT, int m); 
	SWD_PoseNormalization();
	virtual ~SWD_PoseNormalization();

};

#endif // !defined(AFX_SWD_POSENORMALIZATION_H__E64AFAB7_DCBC_4FDB_BA1D_2553405BB9D0__INCLUDED_)
