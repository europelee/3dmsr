// ModelInfo.cpp: implementation of the CModelInfo class.
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

#include "myconfig.h"
#include "ModelInfo.h"
#include "Vertex3d.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


CModelInfo::CModelInfo()
{
  barMesh = new Vertex3d();
   sAllTrian = 0;
  //c++ 引用，参考贝尔实验室，stroustrup
//c++ 程序设计语言 翻译版 
//double *areaT = new double[mesh->m_nfaces];//存储每个三角形面积
//Vertex3d *veGArray = new Vertex3d[mesh->m_nfaces];//存储每个三角形重心,是否分配new??,debug
  for(int i= 0 ; i<3 ; i++)
	  for(int j= 0; j<3; j++)
	  {
	  
        CovMat[i][j] = 0;//协方差矩阵
        eigenvectors[i][j] = 0; 
         F[i][j] =0;
	  }
for(int i= 0 ; i<3 ; i++)
 eigenvalues[i] = 0;

 
  scaleV = 0;
}

CModelInfo::~CModelInfo()
{
  delete barMesh;
  barMesh = NULL;
  delete []areaT;
  areaT = NULL;
  delete []veGArray;
  veGArray = NULL;
}
