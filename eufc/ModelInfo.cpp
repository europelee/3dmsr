// ModelInfo.cpp: implementation of the CModelInfo class.
//
/*-----------------------------------------------------------------
*\ brief   aout model info 
* \ author   ��ŷ��	
*\ date      25-9-2008
*
*
*
*
*------------------------------------------------------------------*/
//history
/*-----------------------------------------------------------------
*brief  modify for 3DQModelRetrieval
*author ��ŷ��
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
  //c++ ���ã��ο�����ʵ���ң�stroustrup
//c++ ����������� ����� 
//double *areaT = new double[mesh->m_nfaces];//�洢ÿ�����������
//Vertex3d *veGArray = new Vertex3d[mesh->m_nfaces];//�洢ÿ������������,�Ƿ����new??,debug
  for(int i= 0 ; i<3 ; i++)
	  for(int j= 0; j<3; j++)
	  {
	  
        CovMat[i][j] = 0;//Э�������
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
