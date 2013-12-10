// DepthMap.h: feature for 3d shape rerieval.
//
//////////////////////////////////////////////////////////////////////

/*--------------------------------------------------------------
*
*\brief   about mesh depthmap 
*\author   ��ŷ��
*\date     18-4-2009
*
*
*
------------------------------------------------------------------*/
#ifndef DEPTHMAPP_H
#define DEPTHMAPP_H

#ifdef EUDIM_DLLAPI
#else
#define EUDIM_DLLAPI _declspec(dllimport)
#endif

#ifndef MYCONFIG_H //����Ƶ���򿪣����ر��ļ�����
#include "myconfig.h"
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////
//#define DepthMapSize 128
/////////////////////////////////////////////////////////////////////////////////////////////////////

class Mesh;
class Face;
class Vertex3d;

class EUDIM_DLLAPI EU_DepthMapp
{
    protected:
         float wEBB;
		 float cxEBB;
		 float cyEBB;
		 float czEBB;
        //static const 
			int planeNum;//6
       // static const 
			int mapSize;//256
      //  const 
			//float w;//1
        float stepLength;
    public:
		float vranicfft[438];
        EU_DepthMapp();
        virtual ~EU_DepthMapp();
		float*** depMapValue;//�������Χ�У�ÿ�����϶�Ӧ��image
        void ComputeDepMap(const char *filename);
		void SaveFeature2File();//Ĭ�ϵ�c:\\myvranicdimdll.txt,Ϊ����Ŀ��
    private:
        void HelpComputeDep(Face& face);
        void HelpLineBresenham(float tempdata[DmSize][DmSize],int index[3][2]);
		void HelpLineBresenham2(float tempdata[DmSize][DmSize],int index[3][2]);
		void Bresenham_FromGame(float tempdata[DmSize][DmSize],int index[3][2],
			int& x0,int & y0,int& x1, int& y1);
        void FillTriangle(float tempdata[DmSize][DmSize],int& x1,int& y1,
                              int& x2,int& y2,int index[3][2]);
        void HelpInterpolation(float coeff[3],int index[3][2],int& xi, int& yi);            
        void CompareZBuffer(float** dep2d,float tempdata[DmSize][DmSize],
           int& x1,int& y1, int& x2,int& y2 );
        
};    

#endif
