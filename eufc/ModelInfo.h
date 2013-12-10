// ModelInfo.h: interface for the CModelInfo class.
//
/*-----------------------------------------------------------------
*\ brief   aout model info 
* \ author   李欧州	
*\ date      25-9-2008
*
*
*
*
*------------------------------------------------------------------*/
//history
/*-----------------------------------------------------------------
*brief  modify for 3DQModelRetrieval
*author 李欧州
*date   17-02-2009 
-----------------------------------------------------------------*/
//////////////////////////////////////////////////////////////////////

#ifndef MODELINFO_H
#define MODELINFO_H

class Vertex3d;
class CModelInfo  
{
public:
	
	Vertex3d* barMesh;
	double sAllTrian;
	double *areaT;
	Vertex3d *veGArray;
	double CovMat[3][3];
	double eigenvalues[3];
	double eigenvectors[3][3];
	int F[3][3];
	double scaleV;
	CModelInfo();
	virtual ~CModelInfo();

};

#endif 
