// SWD_PoseNormalization.cpp: implementation of the SWD_PoseNormalization class.
//
//////////////////////////////////////////////////////////////////////

#include "myconfig.h"
#include "SWD_PoseNormalization.h"
#include "Mesh.h"
#include "Face.h"
#include "Vertex3d.h"
//#include "SWD_Maths.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
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
/*
*brief: bug!!,in my_retri.dll,SWD_PoseNormalization.cpp,func,PoseNormalization(Mesh* mesh,double&  scaleV, int F[3][3])
*verts coordinates are divided by scaleV two times!!!!,so head models from lab,retrieval results bad.
*solved!!Then update database!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
* author: europelee 
*date: 03-06-2009 18:04
*/
SWD_PoseNormalization::SWD_PoseNormalization()
{

}

SWD_PoseNormalization::~SWD_PoseNormalization()
{

}
void SWD_PoseNormalization::ComputHelperFlip(double& fx, float x[3], double* areaT, int m)
{
	sort(x,&(x[3]));
	//float minV = SWD_Maths::min3(x[0],x[1],x[2]);
	//float maxV = SWD_Maths::max3(x[0],x[1],x[2]);
	
	if( x[0] >= 0)
	{
		fx = fx + (pow(x[0],2) + pow(x[1],2)
			+ pow(x[2],2) + x[0]*x[1] + x[0]*x[2] + x[1]*x[2] )*areaT[m];
	}
	else
	{
		if( x[2] < 0)
		{
			fx = fx + (-(pow(x[0],2) + pow(x[1],2)
				+ pow(x[2],2) + x[0]*x[1] + x[0]*x[2] + x[1]*x[2] ))*areaT[m];
		} 	  
		else
		{
			if( x[0]*x[1] <= 0)
			{
				float fi = (pow(x[0],2) + pow(x[1],2)
					+ pow(x[2],2) + x[0]*x[1] + x[0]*x[2] + x[1]*x[2] )
					- 2* pow(x[0],4)/((x[1]-x[0])*(x[2]-x[0]));
				
				fx = fx + fi*areaT[m];
				
			}
			else
				if( x[0]*x[1] > 0 ) 
				{
					float temp = x[2]; x[2] = x[0] ; x[0] = temp;
					float fi = (-(pow(x[0],2) + pow(x[1],2)
						+ pow(x[2],2) + x[0]*x[1] + x[0]*x[2] + x[1]*x[2] ))
						+ 2* pow(x[0],4)/((x[1]-x[0])*(x[2]-x[0]));
					
					fx = fx + fi*areaT[m];
				}
		}
	}
	
	
}
void SWD_PoseNormalization::ComputFlipInviaCPCA(int F[3][3], Mesh* mesh,double* areaT,double& sAllTrian)
{   
	double fx = 0;
	double fy = 0;
	double fz = 0;
	for(int m=0; m<mesh->m_nfaces; m++)
	{
      Face* fac=&(mesh->m_faces[m]);
	  float x[3] = {fac->m_verts[0]->m_x,fac->m_verts[1]->m_x,fac->m_verts[2]->m_x};
	  float y[3] = {fac->m_verts[0]->m_y,fac->m_verts[1]->m_y,fac->m_verts[2]->m_y};
	  float z[3] = {fac->m_verts[0]->m_z,fac->m_verts[1]->m_z,fac->m_verts[2]->m_z};
      ComputHelperFlip(fx, x, areaT, m);
	  ComputHelperFlip(fy, y, areaT, m);
	  ComputHelperFlip(fz, z, areaT, m);
	}
    fx = fx / (6*sAllTrian);
	fy = fy / (6*sAllTrian);
	fz = fz / (6*sAllTrian);

	if(fx == 0)
    F[0][0] = 0;
	else
    F[0][0] = fx / fabs(fx);
	assert(F[0][0] == 0 || F[0][0] == 1 || F[0][0] == -1);
	if(fy == 0)
		F[1][1] = 0;
	else
		F[1][1] = fy / fabs(fy);
	assert(F[1][1] == 0 || F[1][1] == 1 || F[1][1] == -1);
	if(fz == 0)
		F[2][2] = 0;
	else
		F[2][2] = fz / fabs(fz);
	assert(F[2][2] == 0 || F[2][2] == 1 || F[2][2] == -1);

}
void SWD_PoseNormalization::ComputFlipInvia(int F[3][3], Mesh* mesh,double* areaT,double& sAllTrian)
 {
	 double fx = 0;
	double fy = 0;
	double fz = 0;

for(int m=0; m<mesh->m_nfaces; m++)
{
  Face* fac=&(mesh->m_faces[m]);
  double tempV = 0;
  int signV = 0;
/////////////////test////////////////////////

//  if (m == 1) 
//	  {
//		  printf("\n!!!!!%.10e,%.10e,%.10e,!!!!\n",(*(fac->m_verts[0])).m_x,
//			  (*(fac->m_verts[0])).m_y,(*(fac->m_verts[0])).m_z);
//	  } 


   ///////////////////////////
   tempV = (fac->m_verts[0])->m_x + (fac->m_verts[1])->m_x + (fac->m_verts[2])->m_x;
   if (tempV > 0)
   {
	   signV = 1;
   }

   if (tempV == 0)
   {
	   signV = 0;
   }

   if (tempV < 0)
   {
	   signV = -1;
   } 
  fx += signV * areaT[m] * ((fac->m_verts[0])->m_x + (fac->m_verts[1])->m_x + (fac->m_verts[2])->m_x) *
           ((fac->m_verts[0])->m_x + (fac->m_verts[1])->m_x + (fac->m_verts[2])->m_x) / 9;       
   ////////////////////////////////
  
  ///////////////////////////
   tempV = (fac->m_verts[0])->m_y + (fac->m_verts[1])->m_y + (fac->m_verts[2])->m_y;
   if (tempV > 0)
   {
	   signV = 1;
   }

   if (tempV == 0)
   {
	   signV = 0;
   }

   if (tempV < 0)
   {
	   signV = -1;
   } 
  fy += signV * areaT[m] * ((fac->m_verts[0])->m_y + (fac->m_verts[1])->m_y + (fac->m_verts[2])->m_y) *
           ((fac->m_verts[0])->m_y + (fac->m_verts[1])->m_y + (fac->m_verts[2])->m_y) / 9;       
   ////////////////////////////////

  /////////////////////////////////

  ///////////////////////////
   tempV = (fac->m_verts[0])->m_z + (fac->m_verts[1])->m_z + (fac->m_verts[2])->m_z;
   if (tempV > 0)
   {
	   signV = 1;
   }

   if (tempV == 0)
   {
	   signV = 0;
   }

   if (tempV < 0)
   {
	   signV = -1;
   } 
  fz += signV * areaT[m] * ((fac->m_verts[0])->m_z + (fac->m_verts[1])->m_z + (fac->m_verts[2])->m_z) *
           ((fac->m_verts[0])->m_z + (fac->m_verts[1])->m_z + (fac->m_verts[2])->m_z) / 9;       
   ////////////////////////////////

	
}

fx = fx / sAllTrian;
fy = fy / sAllTrian;
fz = fz / sAllTrian;
////////////////////////////
//printf("\nfx = %.10e\n",fx);
//printf("\nfy = %.10e\n",fy);
//printf("\nfz = %.10e\n",fz);
///////////////////////////


  if (fx > 0)
   {
	   F[0][0] = 1;
   }

   if (fx  == 0)
   {
	   F[0][0] = 0;
   }

   if (fx < 0)
   {
	   F[0][0]  = -1;
   } 
   
   if (fy > 0)
   {
	   F[1][1] = 1;
   }

   if (fy == 0)
   {
	   F[1][1] = 0;
   }

   if (fy < 0)
   {
	   F[1][1]  = -1;
   } 

   if (fz > 0)
   {
	   F[2][2] = 1;
   }

   if (fz  == 0)
   {
	   F[2][2] = 0;
   }

   if (fz < 0)
   {
	   F[2][2]  = -1;
   } 
 }

 void SWD_PoseNormalization::ComputRotaInvia(Mesh* mesh, double eigenvectors[3][3])
 {  
	 //float scaleV1 = 0, scaleV2 = 0;
	for(int i=0;i<mesh->m_nverts;i++)
	{  
	
	Vertex3d* temp3d = new Vertex3d((mesh->m_verts[i]).m_x,
		                            (mesh->m_verts[i]).m_y,(mesh->m_verts[i]).m_z);//运算符重载
	
	
      (mesh->m_verts[i]).m_x = eigenvectors[0][0]*(temp3d->m_x)
		                       +eigenvectors[0][1]*(temp3d->m_y)
							   +eigenvectors[0][2]*(temp3d->m_z);

	  (mesh->m_verts[i]).m_y = eigenvectors[1][0]*(temp3d->m_x)
		                       +eigenvectors[1][1]*(temp3d->m_y)
							   +eigenvectors[1][2]*(temp3d->m_z);

	  (mesh->m_verts[i]).m_z = eigenvectors[2][0]*(temp3d->m_x)
		                       +eigenvectors[2][1]*(temp3d->m_y)
							   +eigenvectors[2][2]*(temp3d->m_z);
	  
	  delete temp3d;
	  temp3d = NULL;

  } 

 }

 void SWD_PoseNormalization::ComputScaleInvia(double& scaleV, Mesh*  mesh,double* areaT,
	                                          double&  sAllTrian, Vertex3d *veGArray)
 {

//////////////////hamid method/////////////////////////////////////
 for(int m=0; m<mesh->m_nverts; m++) 
 {
   Vertex3d* tempVex;
  tempVex = &(mesh->m_verts[m]);
   scaleV += ( tempVex->m_x *  tempVex->m_x + tempVex->m_y *  tempVex->m_y 
	          + tempVex->m_z *  tempVex->m_z );

	 
 }
  
  scaleV = scaleV / (mesh->m_nverts);
  scaleV = 2*sqrt(scaleV);
  
 }

 void SWD_PoseNormalization::ComputTranslatInvia(Mesh* mesh,Vertex3d* barMesh)
{
   for(int i=0;i<mesh->m_nverts;i++)
	{  
	
	Vertex3d* temp3d =((mesh->m_verts[i]) - (*barMesh));//运算符重载
	 mesh->m_verts[i].m_x = temp3d->m_x;
     mesh->m_verts[i].m_y = temp3d->m_y;
	 mesh->m_verts[i].m_z = temp3d->m_z;
	 delete temp3d;
	 temp3d = NULL;
   }
}

  void SWD_PoseNormalization::PoseNormalization(Mesh* mesh,double&  scaleV, int F[3][3])
 {   
	
	 for(int i=0;i<mesh->m_nverts;i++)
 {
	 
		

//******************	翻转不变性暂时搁置**************************
      (mesh->m_verts[i]).m_x = (F[0][0]*((mesh->m_verts[i]).m_x)
		                       +F[0][1]*((mesh->m_verts[i]).m_y)
							   +F[0][2]*((mesh->m_verts[i]).m_z))
							   /scaleV
							   ;

	  (mesh->m_verts[i]).m_y = (F[1][0]*((mesh->m_verts[i]).m_x)
		                       +F[1][1]*((mesh->m_verts[i]).m_y)
							   +F[1][2]*((mesh->m_verts[i]).m_z))
							   /scaleV
							   ;

	  (mesh->m_verts[i]).m_z = (F[2][0]*((mesh->m_verts[i]).m_x)
		                       +F[2][1]*((mesh->m_verts[i]).m_y)
							   +F[2][2]*((mesh->m_verts[i]).m_z))
							   /scaleV
							  ;
	  
} 

 }
