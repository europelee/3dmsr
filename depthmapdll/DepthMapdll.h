// DepthMap.h: feature for 3d shape rerieval.
//
//////////////////////////////////////////////////////////////////////

/*--------------------------------------------------------------
*
*\brief   about mesh depthmap 
*\author   李欧州
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

#ifndef MYCONFIG_H //避免频繁打开，即关闭文件操作
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
		float*** depMapValue;//立方体包围盒，每个面上对应的image
        void ComputeDepMap(const char *filename);
		void SaveFeature2File();//默认到c:\\myvranicdimdll.txt,为测试目的
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
