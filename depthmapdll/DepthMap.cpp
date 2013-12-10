//HelpLineBresenham2 come from "3d game engine design" pdf
//HelpLineBresenham 唐圣哲 计算机图形学
//24-04-2009 00:03 by europelee
#define EUDIM_DLLAPI  _declspec(dllexport)

#include "myconfig.h"
#include"DepthMapdll.h"
#include "Mesh.h"
#include "Face.h"
#include "Vertex3d.h"
#include "fft2.h"
#include "ModelRead.h"
#include "SWD_PoseNormalization.h"
#include "ModelInfo.h"
#include "SWD_Maths.h"
//const int EU_DepthMap::planeNum = 6;
//const int EU_DepthMap::mapSize = DepthMapSize;
//const int EU_DepthMap::w = 1;
//float vranicfft[438];
EU_DepthMapp::EU_DepthMapp()
{
	planeNum = 6;
	mapSize = DmSize;
	// w = 1;

  //stepLength = 2*w/mapSize;
  depMapValue = new float**[planeNum];
  for(int i=0;i<planeNum;i++)
  {
	  int j,k;
      depMapValue[i] = new float*[mapSize];
      for(j=0;j<mapSize;j++)
	  {
         depMapValue[i][j] = new float[mapSize];
      for(k=0; k<mapSize;k++)
          depMapValue[i][j][k] =0.0;
	  }
  }    
       
} 

 EU_DepthMapp::~EU_DepthMapp()
 {
     for(int i=0;i<planeNum;i++)
      for(int j=0;j<mapSize;j++)
      delete[] (depMapValue[i][j]);
      
     for(int i=0;i<planeNum;i++)
        delete[] depMapValue[i];
        
     delete[] depMapValue;
     depMapValue = 0; 
     
 } 
 void EU_DepthMapp::FillTriangle(float tempdata[DmSize][DmSize],int& x1,int& y1,
                              int& x2,int& y2,int index[3][2]) 
 {  
     int k;int j; int left=-1,right=-1; float coeff[3]={0};
	 //assert(x1-1>=0&&x2+1<=127);
    for(int i=y1-1;i<=y2+1;i++)//may be bug happen 
    {
		if(i<0)
			i = 0;
		if(i>DmSize-1)
			break;
      for(j=x1-1,k=x2+1;j<=x2+1,k>=x1-1;j++,k--)
      {
		  if(j<0)
			  j = 0;
		  if(j>DmSize-1)
			  break;
		  if(k>DmSize-1)
			  k = DmSize-1;
		  if(k<0)
			   break;

          if(fabs(tempdata[i][j]-0.0)>0.0000000001&&left<0)
            left = j;
          if(fabs(tempdata[i][k]-0.0)>0.0000000001&&right<0)
            right = k;
           if(left>0&&right>0)
             break;             
      } 
      for(int t=left;t<=right;t++)
      {
          if(fabs(tempdata[i][t]-0.0)<0.0000000001)
          {
             
          HelpInterpolation(coeff,index,i,t);
          tempdata[i][t]=coeff[0]*tempdata[index[0][0]][index[0][1]] +
           coeff[1]*tempdata[index[1][0]][index[1][1]]+
           coeff[2]*tempdata[index[2][0]][index[2][1]]; 
          }    
      }    
	  left = right =-1;//care!!!! I find it,so depthmap over,01:37 24.04.2009
    }       
	 //int k;int j; int left=-1,right=-1; float coeff[3]={0};
	 ////assert(x1-1>=0&&x2+1<=127);
	 //for(int i=0;i<=127;i++)//may be bug happen 
	 //{
		// if(i<0)
		//	 i = 0;
		// if(i>127)
		//	 break;
		// for(j=0,k=127;j<=127,k>=0;j++,k--)
		// {
		//	 if(j<0)
		//		 j = 0;
		//	 if(j>127)
		//		 break;
		//	 if(k>127)
		//		 k = 127;
		//	 if(k<0)
		//		 break;

		//	 if(fabs(tempdata[i][j]-0.0)>0.0000000001&&left<0)
		//		 left = j;
		//	 if(fabs(tempdata[i][k]-0.0)>0.0000000001&&right<0)
		//		 right = k;
		//	 if(left>0&&right>0)
		//		 break;             
		// } 
		// for(int t=left;t<=right;t++)
		// {
		//	 if(fabs(tempdata[i][t]-0.0)<0.0000000001)
		//	 {

		//		 HelpInterpolation(coeff,index,i,t);
		//		 tempdata[i][t]=coeff[0]*tempdata[index[0][0]][index[0][1]] +
		//			 coeff[1]*tempdata[index[1][0]][index[1][1]]+
		//			 coeff[2]*tempdata[index[2][0]][index[2][1]]; 
		//	 }    
		// }  
		// left = right = -1;
	 //}        
 }
 void EU_DepthMapp::CompareZBuffer(float** dep2d,float tempdata[DmSize][DmSize],
           int& x1,int& y1, int& x2,int& y2 )
{
    for(int i=y1-1;i<=y2+1;i++)//may be bug happen 
    { 
		if(i<0)
			i = 0;
		if(i>DmSize-1)
			break;
        for(int j=x1-1;j<=x2+1;j++)
		{
			if(j<0)
				j = 0;
			if(j>DmSize-1)
				break;
          if(fabs(tempdata[i][j]-0.0)>0.0000000001)
           dep2d[i][j] = dep2d[i][j]>tempdata[i][j]?dep2d[i][j]:tempdata[i][j];
		}
    }  
	//for(int i=0;i<=127;i++)//may be bug happen 
	//{ 
	//	if(i<0)
	//		i = 0;
	//	if(i>127)
	//		break;
	//	for(int j=0;j<=127;j++)
	//	{
	//		if(j<0)
	//			j = 0;
	//		if(j>127)
	//			break;
	//		if(fabs(tempdata[i][j]-0.0)>0.0000000001)
	//			dep2d[i][j] = dep2d[i][j]>tempdata[i][j]?dep2d[i][j]:tempdata[i][j];
	//	}
	//}  

    
}    
 void EU_DepthMapp::HelpInterpolation(float coeff[3],int index[3][2],int& xi, int& yi)
 {
   coeff[0] = abs(index[1][0]*index[2][1]-index[2][0]*index[1][1]+yi*index[2][0]
    -xi*index[2][1]+xi*index[1][1]-yi*index[1][0]);
   coeff[1] = abs(index[0][0]*index[2][1]-index[2][0]*index[0][1]+yi*index[2][0]
    -xi*index[2][1]+xi*index[0][1]-yi*index[0][0]);
   coeff[2] = abs(index[0][0]*index[1][1]-index[1][0]*index[0][1]+yi*index[1][0]
    -xi*index[1][1]+xi*index[0][1]-yi*index[0][0]); 
   float s = coeff[0]+coeff[1]+coeff[2];
   coeff[0] = coeff[0]/s;
   coeff[1] = coeff[1]/s;
   coeff[2] = coeff[2]/s; 
     
 }
 void EU_DepthMapp::Bresenham_FromGame(float tempdata[DmSize][DmSize],int index[3][2],
	 int& x0,int & y0,int& x1, int& y1)
 {
   //start point of line
	 int x =x0,y=y0;
	 //direction of line
	 int dx = x1-x0,dy = y1-y0;
	 //increment or decrement depending on direction of line
	 int sx,sy;
	 if (dx>0)
        sx = 1;
	 else
		 if (dx<0)
		 {
			 sx = -1;
			 dx = -dx;
		 }
		 else
			 sx = 0;
	 if (dy>0)
		 sy = 1;
	 else
		 if (dy<0)
		 {
			 sy = -1;
			 dy = -dy;
		 }
		 else
			 sy = 0;
	 int ax = 2*dx,ay = 2*dy;
	 if (dy<=dx)
	 {
		 //single step in x-direction
		 for (int decy=ay-dx; /**/; x += sx, decy += ay)
		 {
			 float coeff[3]={0};
			 HelpInterpolation(coeff,index,y,x);
			 tempdata[y][x]=coeff[0]*tempdata[index[0][0]][index[0][1]] +
				 coeff[1]*tempdata[index[1][0]][index[1][1]]+
				 coeff[2]*tempdata[index[2][0]][index[2][1]]; 
			 /////////////////////////////////////////////////////////////////////////////////////
			 //take bresenham step
			 if (x == x1)
			   break;
			 if (decy >=0)
			 {
				 decy -= ax;
				 y += sy;
			 }
		 }//for (int decy=ay-dx;x += sx; decy += ay) 
	 }// if (dy<=dx)
	 else
	 {
		 //single step in y-direction
		 for (int decx = ax-dy;  ;y +=sy, decx +=ax)
		 {
			 float coeff[3]={0};
			 HelpInterpolation(coeff,index,y,x);
			 tempdata[y][x]=coeff[0]*tempdata[index[0][0]][index[0][1]] +
				 coeff[1]*tempdata[index[1][0]][index[1][1]]+
				 coeff[2]*tempdata[index[2][0]][index[2][1]]; 
			 /////////////////////////////////////////////////////////////////////////////////////
			 //take bresenham step
			 if (y == y1)
				 break;
			 if (decx >=0)
			 {
				 decx -= ay;
				 x += sx;
			 }
		 }
	 }

 }
 void EU_DepthMapp::HelpLineBresenham2(float tempdata[DmSize][DmSize],int index[3][2])
 {
	 for(int i=0;i<3;i++)
	 {
		 int j = (i+1)%3;
		 if(index[i][0]<index[j][0])
         Bresenham_FromGame(tempdata,index,index[i][1],index[i][0],index[j][1],index[j][0]);
		 else
         Bresenham_FromGame(tempdata,index,index[j][1],index[j][0],index[i][1],index[i][0]);
	 }
 }
void EU_DepthMapp::HelpLineBresenham(float tempdata[DmSize][DmSize],int index[3][2])
{
    for(int i=0;i<3;i++)
    {
        int j = (i+1)%3;
        int x = index[i][1];//col
        int y = index[i][0];//row
        int s1,s2;
        int deltax = abs(index[j][1]-index[i][1]);
        int deltay = abs(index[j][0]-index[i][0]);
        if(index[j][1]-index[i][1]>=0) 
           s1=1;
         else
            s1=-1;
         if(index[j][0]-index[i][0]>=0)     
           s2=1;
         else
            s2 = -1;
         int interchange;   
         if(deltay>deltax)
           {
             int temp = deltax;
             deltax = deltay;
             deltay = temp;
             interchange = 1;  
           }
          else
            interchange = 0;
        int f = 2*deltay-deltax;
        for(int k=1;k<deltax;k++)
        {
          float coeff[3]={0};
          HelpInterpolation(coeff,index,y,x);
          tempdata[y][x]=coeff[0]*tempdata[index[0][0]][index[0][1]] +
           coeff[1]*tempdata[index[1][0]][index[1][1]]+
           coeff[2]*tempdata[index[2][0]][index[2][1]]; 
          /////////////////////////////////////////////////////////
           if(f>=0)
           {
             if(interchange==1)
              x = x+s1;
             else
              y = y+s2;
              f = f-2*deltax;  
           }
           if(interchange==1)
             y= y+s2;
           else
             x= x+s1;
           f = f+2*deltay;                 
        }                
    }    
    
    
}      
void EU_DepthMapp::HelpComputeDep(Face& face)
{   
   int index[3][2] = {0}; 
   //float tempdepth[3]={0};
   //存放min,max
   int min_1d=mapSize-1, max_1d=0,min_2d=mapSize-1, max_2d=0; 
   float tempImage[DmSize][DmSize] = {0};
   float tempmax[3]= {0};
   //x_min clip plane     
    for(int i=0;i<3;i++)
    {                
     int a,b;
     a = mapSize*(((*(face.m_verts[i])).m_y + wEBB-cyEBB )/ (2*wEBB));
     b = mapSize*(((*(face.m_verts[i])).m_z + wEBB-czEBB)/ (2*wEBB));
	 if(a<0)
		  a = 0;
	 if (a>mapSize-1)
	 {
		 a = mapSize-1;
	 }
	 if(b<0)
		 b = 0;
	 if (b>mapSize-1)
	 {
		 b = mapSize-1;
	 }
     //assert(a>=0&&a<=mapSize-1);
     //assert(b>=0&&b<=mapSize-1);
     index[i][0]=b;//row
     if(b<min_2d)
        min_2d = b;
     if(b>max_2d)
        max_2d = b;
     index[i][1]=a;//col
      if(a<min_1d)
        min_1d = a;
     if(a>max_1d)
        max_1d = a; 
     tempImage[b][a] = (wEBB+cxEBB- (*(face.m_verts[i])).m_x)/(2*wEBB);
     //max
     tempmax[i] = (wEBB+(*(face.m_verts[i])).m_x-cxEBB)/(2*wEBB);
    
    } 
    //三角形填充
   ////////////////////////////////////////////////////////////////////////////////            
     HelpLineBresenham2(tempImage,index);//三角形边的光栅化
     //内部填充
     FillTriangle(tempImage,min_1d,min_2d,max_1d,
                              max_2d,index);  
      
    CompareZBuffer(depMapValue[0],tempImage,
           min_1d,min_2d,max_1d,max_2d );
     
      for(int i=0;i<mapSize;i++)           
        for(int j=0;j<mapSize;j++) 
             tempImage[i][j]=0;
       
       for(int i=0;i<3;i++)
          tempImage[index[i][0]][index[i][1]] = tempmax[i];
          
         HelpLineBresenham2(tempImage,index);//三角形边的光栅化
     //内部填充
     FillTriangle(tempImage,min_1d,min_2d,max_1d,
                              max_2d,index);  
      
    CompareZBuffer(depMapValue[1],tempImage,
           min_1d,min_2d,max_1d,max_2d );               
    ////////////////////////////////////////////////////////////////////////////  
      //y_min clip plane
      //y_max clip plane
	for(int i=0;i<mapSize;i++)           
		for(int j=0;j<mapSize;j++) 
			tempImage[i][j]=0;
	for(int i=0;i<3;i++)
	{
		tempmax[i] = 0;
		for(int j=0;j<2;j++)
			index[i][j] = 0;
	}
    min_1d=mapSize-1, max_1d=0,min_2d=mapSize-1, max_2d=0; 
    
      for(int i=0;i<3;i++)
    {                   
     int a,b;
     a = mapSize*(((*(face.m_verts[i])).m_z + wEBB-czEBB )/ (2*wEBB));
     b = mapSize*(((*(face.m_verts[i])).m_x + wEBB-cxEBB )/ (2*wEBB));
	 if(a<0)
		 a = 0;
	 if (a>mapSize-1)
	 {
		 a = mapSize-1;
	 }
	 if(b<0)
		 b = 0;
	 if (b>mapSize-1)
	 {
		 b = mapSize-1;
	 }
     index[i][0]=b;
	 if(b<min_2d)
		 min_2d = b;
	 if(b>max_2d)
		 max_2d = b;
     index[i][1]=a;
	 if(a<min_1d)
		 min_1d = a;
	 if(a>max_1d)
		 max_1d = a; 
	 tempImage[b][a] = (wEBB- (*(face.m_verts[i])).m_y+cyEBB)/(2*wEBB);
	 //max
	 tempmax[i] = (wEBB+(*(face.m_verts[i])).m_y-cyEBB)/(2*wEBB);

     //float temp = (w- (*(face.m_verts[i])).m_y)/(2*w);
     //depMapValue[2][b][a] =  depMapValue[2][b][a]> temp? depMapValue[2][b][a]:temp;
     ////x_max clip plane
     //temp = (w+(*(face.m_verts[i])).m_y)/(2*w);
     //depMapValue[3][b][a] =  depMapValue[3][b][a]> temp? depMapValue[3][b][a]:temp;
    } 
    //三角形填充
	  ////////////////////////////////////////////////////////////////////////////////            
	  HelpLineBresenham2(tempImage,index);//三角形边的光栅化
	  //内部填充
	  FillTriangle(tempImage,min_1d,min_2d,max_1d,
		  max_2d,index);  

	  CompareZBuffer(depMapValue[2],tempImage,
		  min_1d,min_2d,max_1d,max_2d );

	  for(int i=0;i<mapSize;i++)           
		  for(int j=0;j<mapSize;j++) 
			  tempImage[i][j]=0;

	  for(int i=0;i<3;i++)
		  tempImage[index[i][0]][index[i][1]] = tempmax[i];

	  HelpLineBresenham2(tempImage,index);//三角形边的光栅化
	  //内部填充
	  FillTriangle(tempImage,min_1d,min_2d,max_1d,
		  max_2d,index);  

	  CompareZBuffer(depMapValue[3],tempImage,
		  min_1d,min_2d,max_1d,max_2d );               
	  ////////////////////////////////////////////////////////////////////////////  
   ///////////////////////////////////////////////////////////////////////////////   
      //z_min clip plane
      //z_max clip plane 
	  for(int i=0;i<mapSize;i++)           
		  for(int j=0;j<mapSize;j++) 
			  tempImage[i][j]=0;
	  for(int i=0;i<3;i++)
	  {
		  tempmax[i] = 0;
		  for(int j=0;j<2;j++)
			  index[i][j] = 0;
	  }
	  min_1d=mapSize-1, max_1d=0,min_2d=mapSize-1, max_2d=0; 

	  for(int i=0;i<3;i++)
	  {                   
		  int a,b;
		  a = mapSize*(((*(face.m_verts[i])).m_x + wEBB-cxEBB )/ (2*wEBB));
		  b = mapSize*(((*(face.m_verts[i])).m_y + wEBB-cyEBB )/ (2*wEBB));
		  if(a<0)
			  a = 0;
		  if (a>mapSize-1)
		  {
			  a = mapSize-1;
		  }
		  if(b<0)
			  b = 0;
		  if (b>mapSize-1)
		  {
			  b = mapSize-1;
		  }
		  index[i][0]=b;
		  if(b<min_2d)
			  min_2d = b;
		  if(b>max_2d)
			  max_2d = b;
		  index[i][1]=a;
		  if(a<min_1d)
			  min_1d = a;
		  if(a>max_1d)
			  max_1d = a; 
		  tempImage[b][a] = (wEBB- (*(face.m_verts[i])).m_z+czEBB)/(2*wEBB);
		  //max
		  tempmax[i] = (wEBB+(*(face.m_verts[i])).m_z-czEBB)/(2*wEBB);

		  	  } 
	  //三角形填充
	  ////////////////////////////////////////////////////////////////////////////////            
	  HelpLineBresenham2(tempImage,index);//三角形边的光栅化
	  //内部填充
	  FillTriangle(tempImage,min_1d,min_2d,max_1d,
		  max_2d,index);  

	  CompareZBuffer(depMapValue[4],tempImage,
		  min_1d,min_2d,max_1d,max_2d );

	  for(int i=0;i<mapSize;i++)           
		  for(int j=0;j<mapSize;j++) 
			  tempImage[i][j]=0;

	  for(int i=0;i<3;i++)
		  tempImage[index[i][0]][index[i][1]] = tempmax[i];

	  HelpLineBresenham2(tempImage,index);//三角形边的光栅化
	  //内部填充
	  FillTriangle(tempImage,min_1d,min_2d,max_1d,
		  max_2d,index);  
	  CompareZBuffer(depMapValue[5],tempImage,
		  min_1d,min_2d,max_1d,max_2d );      
}       
void EU_DepthMapp::ComputeDepMap(const char *filename)
{ 
	Mesh * mesh = ModelRead::ReadOffFile(filename);
	if(mesh==NULL)
		return;
	///////////////////////////////////////////////cpca
	CModelInfo modelInfo;
	SWD_PoseNormalization posNormal;	
	modelInfo.areaT = new double[(mesh)->m_nfaces];
	modelInfo.veGArray = new Vertex3d[(mesh)->m_nfaces];
	SWD_Maths::ComputBarMesh(modelInfo.areaT,modelInfo.veGArray,mesh,
		modelInfo.barMesh,modelInfo.sAllTrian);	
	posNormal.ComputTranslatInvia(mesh,modelInfo.barMesh);
	SWD_Maths::ComputCovMat(modelInfo.areaT, modelInfo.veGArray,mesh
		,modelInfo.CovMat, modelInfo.barMesh, modelInfo.sAllTrian);
	SWD_Maths::Jacobi_Cyclic_Method(modelInfo.eigenvalues, (double *)(modelInfo.eigenvectors)
		,(double *)(modelInfo.CovMat), 3);
	SWD_Maths::SortEigen(modelInfo.eigenvalues, modelInfo.eigenvectors);	
	posNormal.ComputRotaInvia(mesh, modelInfo.eigenvectors);
	posNormal.ComputFlipInviaCPCA(modelInfo.F,mesh,
		modelInfo.areaT, modelInfo.sAllTrian);
	posNormal.ComputScaleInvia(modelInfo.scaleV,mesh,modelInfo.areaT
		,modelInfo.sAllTrian, modelInfo.veGArray);
	posNormal.PoseNormalization(mesh, modelInfo.scaleV,modelInfo.F);	
	////////////////////////////////////////////////cpca
	float x1=30000,x2=-30000,y1=30000,y2=-30000,z1=30000,z2=-30000;
	 for(int i=0;i<mesh->m_nverts;i++)
	 {
        if(x1>mesh->m_verts[i].m_x)
			x1 = mesh->m_verts[i].m_x;
		if(x2<mesh->m_verts[i].m_x)
			x2 = mesh->m_verts[i].m_x;
		if(y1>mesh->m_verts[i].m_y)
			y1 = mesh->m_verts[i].m_y;
		if(y2<mesh->m_verts[i].m_y)
			y2 = mesh->m_verts[i].m_y;
		if(z1>mesh->m_verts[i].m_z)
			z1=mesh->m_verts[i].m_z;
		if(z2<mesh->m_verts[i].m_z)
			z2=mesh->m_verts[i].m_z;
	 }
	 wEBB=MAXMAX((x2-x1)/2,(y2-y1)/2,(z2-z1)/2);
	 cxEBB=(x1+x2)/2;cyEBB=(y1+y2)/2;czEBB=(z1+z2)/2;

	stepLength = 2*wEBB/mapSize;
    for(int m=0; m<mesh->m_nfaces; m++) 
    {
       HelpComputeDep(mesh->m_faces[m]);      
    }//for(int m=0; m<mesh->m_nfaces; m++)
    
	/////////////////////////////////////////////////////////test
	for(int i=0;i<DmSize;i++)
		for(int j=0;j<DmSize/2;j++)
		{
			float temp = depMapValue[0][i][j];
			depMapValue[0][i][j] = depMapValue[0][i][DmSize-1-j];
			depMapValue[0][i][DmSize-1-j] = temp;

			temp = depMapValue[3][i][j];
			depMapValue[3][i][j] = depMapValue[3][i][DmSize-1-j];
			depMapValue[3][i][DmSize-1-j] = temp;
            
			temp = depMapValue[5][i][j];
			depMapValue[5][i][j] = depMapValue[5][i][DmSize-1-j];
			depMapValue[5][i][DmSize-1-j] = temp;

		}
		
		/////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////
      
		_complex** fftinput = new _complex*[DmSize];//[DmSize][DmSize];
        for (int i=0;i<DmSize;i++)
			fftinput[i] = new _complex[DmSize];
		float** tempmat = new float*[DmSize];
		for (int i=0;i<DmSize;i++)
		{
			tempmat[i] = new float[DmSize];
		}
		for(int k=0;k<6;k++)
		{
			int a,b;
		for (int i=0;i<DmSize;i++)
			for (int j=0;j<DmSize;j++)
			{  
				a = (i+DmSize/2)%DmSize;
				b = ( j+DmSize/2)%DmSize;
                fftinput[a][b].x = depMapValue[k][i][j];
				fftinput[a][b].y = 0;
			}    
			FFT2D(fftinput,DmSize,DmSize,1);
			
			for (int i=0;i<DmSize;i++)
				for (int j=0;j<DmSize;j++)
				{
					a = (i+DmSize/2)%DmSize;
					b = ( j+DmSize/2)%DmSize;
					tempmat[a][b] = sqrt(fftinput[i][j].x*fftinput[i][j].x+fftinput[i][j].y*fftinput[i][j].y);
				}
				/////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////////////////////////
				int i2=0;
                for(int i=8;i>=0;i--)
				{
					int tt ;	
					  if(i==8)
					  {
                            for(tt=0;tt<=i;tt++)
					       {
							   vranicfft[i2+k*73] = tempmat[DmSize/2][DmSize/2-tt];	
                               		i2++;
					       }
					  }
					  else
					  {
                         for(tt=-i;tt<=i;tt++)
						 {
                             vranicfft[i2+k*73] = tempmat[DmSize/2-i+8][DmSize/2+tt];
							 i2++;
						 }
					  }
			
				}
		}
			for (int i=0;i<DmSize;i++)
			{
			delete [] fftinput[i] ;
			}
			delete []fftinput;
			for (int i=0;i<DmSize;i++)
			{
				delete [] tempmat[i] ;
			}
			delete []tempmat;
			if(mesh)
			{
				delete mesh;
				mesh = NULL;
			}
			//////////////////////////////////////////////////////////////
			//FILE* pff;
			//pff = fopen("c:\\myvranicdimdll.txt","w");
   //        for(int i=0;i<438;i++)//438 need insteaded
		 //  {
			//   fprintf(pff,"%.10f\n",vranicfft[i]);
		 //  }
		 //  fclose(pff);
		//////////////////////////////////////////////////////////////////////////////////////////////////////
    
}      
void EU_DepthMapp::SaveFeature2File()
{
	FILE* pff;
	pff = fopen("c:\\myvranicdimdll.txt","w");
	for(int i=0;i<438;i++)//438 need insteaded
	{
		fprintf(pff,"%.10f\n",vranicfft[i]);
	}
	fclose(pff);
}