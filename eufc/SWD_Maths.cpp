// SWD_Maths.cpp: implementation of the SWD_Maths class.
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

#include "myconfig.h"
#include "SWD_Maths.h"
#include "Mesh.h"
#include "Vertex3d.h"
#include "SpheCoord.h"
#include "Face.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
bool SWD_Maths::low_compare1(const SpheCoord & m1, const SpheCoord & m2) {
        return (m1.m_long) < m2.m_long;
}

bool SWD_Maths::up_compare1(const SpheCoord & m2,const SpheCoord & m1) {
        return (m1.m_long) > m2.m_long;
}

bool SWD_Maths::low_compare2(const SpheCoord & m1, const SpheCoord & m2) {
        return (m1.m_lat) < m2.m_lat;
}

bool SWD_Maths::up_compare2(const SpheCoord & m2,const SpheCoord & m1) {
        return (m1.m_lat) > m2.m_lat;
}

bool SWD_Maths::less_second(const SpheCoord & m1, const SpheCoord & m2) {
        return (m1.m_lat) < (m2.m_lat);
}

float SWD_Maths::max3(float a, float b, float c)
{
	float d = (a>b)?a:b;
	return (c>d)?c:d;
}

float SWD_Maths::min3(float a, float b, float c)
{
	float d = (a<b)?a:b;
	return (c<d)?c:d;
}

//area one triangle
double SWD_Maths::AreaOneTriangle(Vertex3d v1,Vertex3d v2,Vertex3d v3)
{ 
  //默认是三角网格
  double vec1[3];
  double vec2[3];
  
  //计算向量
  vec1[0] = v2.m_x - v1.m_x;
  vec1[1] = v2.m_y - v1.m_y;
  vec1[2] = v2.m_z - v1.m_z;
  
  vec2[0] = v3.m_x - v1.m_x;
  vec2[1] = v3.m_y - v1.m_y;
  vec2[2] = v3.m_z - v1.m_z;

  double a = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  double b = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  double c = vec1[0]*vec2[1] - vec1[1]*vec2[0];

  return 0.5*sqrt(a*a+b*b+c*c);

}

//计算三角型重心
Vertex3d* SWD_Maths::CenGravOneTriangle(Vertex3d v1, Vertex3d v2, Vertex3d v3)
{
	Vertex3d* vertex = new Vertex3d();
	
	vertex->m_x = (v1.m_x + v2.m_x + v3.m_x)/3;
    vertex->m_y = (v1.m_y + v2.m_y + v3.m_y)/3;
    vertex->m_z = (v1.m_z + v2.m_z + v3.m_z)/3;

	return vertex;
}

//计算模型的重心，和存储每个三角面片的面积和重心，和总面积
void SWD_Maths::ComputBarMesh(double* areaT, Vertex3d * veGArray,Mesh* mesh,Vertex3d* barMesh
				   , double& sAllTrian)
{
	for(int m=0; m<mesh->m_nfaces; m++)
  {   

	  Face* fac=&(mesh->m_faces[m]);
     
	  double sTriangle = AreaOneTriangle(*(fac->m_verts[0]), *(fac->m_verts[1]), *(fac->m_verts[2]));
      
	  areaT[m] = sTriangle;

	  Vertex3d* veG = CenGravOneTriangle(*(fac->m_verts[0]), *(fac->m_verts[1]), *(fac->m_verts[2]));
      
	  (veGArray[m]).m_x = veG->m_x;
	  (veGArray[m]).m_y = veG->m_y;
	  (veGArray[m]).m_z = veG->m_z;
      
	  

	  barMesh->m_x += sTriangle * (veG->m_x);
      barMesh->m_y += sTriangle * (veG->m_y);
	  barMesh->m_z += sTriangle * (veG->m_z);
      
      delete veG;

	  (sAllTrian) += sTriangle;

  }
   
  barMesh->m_x = barMesh->m_x/sAllTrian;
  barMesh->m_y = barMesh->m_y/sAllTrian;
  barMesh->m_z = barMesh->m_z/sAllTrian;

  
//	printf("%f,%f,%f,\n",barMesh->m_x,barMesh->m_y,barMesh->m_z);
  
//////////////////////////////////////////////////////////////
}

void SWD_Maths::ComputHelperCovMat(Vertex3d vert, double arrayF[3][3])
{
	  arrayF[0][0] =  (vert.m_x) * (vert.m_x);
	  arrayF[0][1] =  (vert.m_x) * (vert.m_y);
	  arrayF[0][2] =  (vert.m_x) * (vert.m_z);

	  arrayF[1][0] =  (vert.m_y) * (vert.m_x);
	  arrayF[1][1] =  (vert.m_y) * (vert.m_y);
	  arrayF[1][2] =  (vert.m_y) * (vert.m_z);

	  arrayF[2][0] =  (vert.m_z) * (vert.m_x);
	  arrayF[2][1] =  (vert.m_z) * (vert.m_y);
	  arrayF[2][2] =  (vert.m_z) * (vert.m_z);
}

void SWD_Maths::ComputHelperMat(double CovMat[3][3], double arrayF[3][3])
{
	  CovMat[0][0] += arrayF[0][0];
	  CovMat[0][1] += arrayF[0][1];
	  CovMat[0][2] += arrayF[0][2];

	  CovMat[1][0] += arrayF[1][0];
	  CovMat[1][1] += arrayF[1][1];
	  CovMat[1][2] += arrayF[1][2];

	  CovMat[2][0] += arrayF[2][0];
	  CovMat[2][1] += arrayF[2][1];
	  CovMat[2][2] += arrayF[2][2];
}





//计算CovMat
 void SWD_Maths::ComputCovMat(double* areaT, Vertex3d * veGArray,Mesh* mesh, double CovMat[3][3]
	 ,Vertex3d* barMesh,double& sAllTrian)
 {   
	 double arrayF[3][3] ={0};
	 //计算CovMat
  for(int m=0; m<mesh->m_nfaces; m++)
  {   
	  double CovMat2[3][3] = {0};

	  Face* fac=&(mesh->m_faces[m]);

	  ComputHelperCovMat(*(fac->m_verts[0]), arrayF );
      ComputHelperMat(CovMat2, arrayF);

      ComputHelperCovMat(*(fac->m_verts[1]), arrayF );
      ComputHelperMat(CovMat2, arrayF);

	  ComputHelperCovMat(*(fac->m_verts[2]), arrayF );
      ComputHelperMat(CovMat2, arrayF);
	  
	  Vertex3d* temp3d =(veGArray[m] - (*barMesh));//运算符重载
	  ComputHelperCovMat(*temp3d, arrayF);
      delete temp3d;
	  temp3d = NULL;
	  
	  for (int i =0 ; i<3; i++)
	  for (int k =0 ; k<3; k++)
		  arrayF[i][k] = 9 * arrayF[i][k];
	  
      ComputHelperMat(CovMat2, arrayF);
      
	  for (int i =0 ; i<3; i++)
	  for (int k =0 ; k<3; k++)
	  CovMat[i][k] += areaT[m] * CovMat2[i][k];
	  
  }
  


  for (int i =0 ; i<3; i++)
	  for (int k =0 ; k<3; k++)
		  CovMat[i][k] = CovMat[i][k]/(12 * sAllTrian);
//printf("covmat is ... \n");
//  for ( i =0 ; i<3; i++)
//  {
//	  for (int k =0 ; k<3; k++)
//	  {
//		  printf("%.10e ,",CovMat[i][k]);//.10e 科学计数法
//	  }
//	  printf("\n");
//  }
 }

 void SWD_Maths::ConvertCartesianToSpherical(Vertex3d& Pos, float& rLong, float& rLat)
{
	
	float posLen = sqrt( Pos.m_x* Pos.m_x +  Pos.m_y* Pos.m_y + Pos.m_z* Pos.m_z);
	posLen = Pos.m_z / posLen ;
    
	    rLat =  asin( posLen );

	float rVal =  sqrt( Pos.m_x* Pos.m_x +  Pos.m_y* Pos.m_y );
	if( rVal==0 )
		rLong = 0;
	else
	{
		rVal = Pos.m_x / rVal;
		
		if( rVal < -1 )
			rVal = -1;
		if( rVal > 1 )
			rVal = 1;
		rLong =  acos( rVal );
		if( Pos.m_y<0  )
			rLong = SWD_TWOPI- rLong;
	}
}

 void SWD_Maths::SortEigen(double eigenvalues[3], double eigenvectors[3][3])
 {
	 double med_eigenValue;
double med_eigenVector;

for (int i =0 ; i<3; i++)
	  for (int k =i ; k<3; k++)
	  {   
		  double tempvector;

		  tempvector = eigenvectors[k][i];
		  eigenvectors[k][i] = eigenvectors[i][k];
          eigenvectors[i][k] = tempvector;
	  }

for (int i =0 ; i<2; i++)
	  for (int ii = i+1; ii<3; ii++)
	  {   
		  if (eigenvalues[i] < eigenvalues[ii]) 
		  {
			  med_eigenValue = eigenvalues[i];
			  eigenvalues[i] = eigenvalues[ii];
			  eigenvalues[ii] = med_eigenValue;
			  for (int k=0;k<3;k++) 
			  {
				  med_eigenVector = eigenvectors[i][k];
				  eigenvectors[i][k] = eigenvectors[ii][k];
				  eigenvectors[ii][k] = med_eigenVector;
			  }
		  }
	  }
       //euclidean unit length
	  for(int i=0 ; i<2; i++)
	  {   
		  double tt = sqrt(eigenvectors[i][0]*eigenvectors[i][0]
			  +eigenvectors[i][1]*eigenvectors[i][1]
			  +eigenvectors[i][2]*eigenvectors[i][2]);
		  for(int j= 0; j<2; j++)
		  {
			  eigenvectors[i][j] = eigenvectors[i][j]/tt;
		  }
	}
//	  for ( i =0 ; i<3; i++)
//        printf("VALUE IS %.10e ,",eigenvalues[i]);
//       printf("\n");

//	  for ( i =0 ; i<3; i++)
//  {
//	  for (int k =0 ; k<3; k++)
//	  {
//		  printf("%.10e ,",eigenvectors[i][k]);//.10e 科学计数法
//	  }
//	  printf("\n");
//  }
 }

 void SWD_Maths::ComputFaceNormal(Face& face)
 {
	
      Vertex3d* v1 = new Vertex3d();
	  Vertex3d* v2 = new Vertex3d();
      
	  v1->m_x = (face.m_verts[1])->m_x - (face.m_verts[0])->m_x;
	  v1->m_y = (face.m_verts[1])->m_y - (face.m_verts[0])->m_y;
	  v1->m_z = (face.m_verts[1])->m_z - (face.m_verts[0])->m_z;

	  v2->m_x = (face.m_verts[2])->m_x - (face.m_verts[0])->m_x;
	  v2->m_y = (face.m_verts[2])->m_y - (face.m_verts[0])->m_y;
	  v2->m_z = (face.m_verts[2])->m_z - (face.m_verts[0])->m_z;

      face.m_normal[0] = v1->m_y*v2->m_z - v1->m_z*v2->m_y;
      face.m_normal[1] = v1->m_z*v2->m_x - v1->m_x*v2->m_z;
	  face.m_normal[2] = v1->m_x*v2->m_y - v1->m_y*v2->m_x;
	  delete v1;
	  v1 = NULL;
	  delete v2;
	  v2 = NULL;
	  float squared_normal_length = 0.0;
      squared_normal_length += face.m_normal[0]*face.m_normal[0];
      squared_normal_length += face.m_normal[1]*face.m_normal[1];
      squared_normal_length += face.m_normal[2]*face.m_normal[2];
      float normal_length = sqrt(squared_normal_length);
     // if (normal_length > 1.0E-6) {
        face.m_normal[0] /= normal_length;
        face.m_normal[1] /= normal_length;
        face.m_normal[2] /= normal_length;
 }

 void SWD_Maths::Jacobi_Cyclic_Method(double eigenvalues[], double *eigenvectors,
                                                    double *A, int n)
{
   int row, i, j, k, m;
   double *pAk, *pAm, *p_r, *p_e;
   double threshold_norm;
   double threshold;
   double tan_phi, sin_phi, cos_phi, tan2_phi, sin2_phi, cos2_phi;
   double sin_2phi, cos_2phi, cot_2phi;
   double dum1;
   double dum2;
   double dum3;
   double r;
   double max;

                  // Take care of trivial cases

   if ( n < 1) return;
   if ( n == 1) {
      eigenvalues[0] = *A;
      *eigenvectors = 1.0;
      return;
   }

          // Initialize the eigenvalues to the identity matrix.

   for (p_e = eigenvectors, i = 0; i < n; i++)
      for (j = 0; j < n; p_e++, j++)
         if (i == j) *p_e = 1.0; else *p_e = 0.0;
  
            // Calculate the threshold and threshold_norm.
 
   for (threshold = 0.0, pAk = A, i = 0; i < ( n - 1 ); pAk += n, i++) 
      for (j = i + 1; j < n; j++) threshold += *(pAk + j) * *(pAk + j);
   threshold = sqrt(threshold + threshold);
   threshold_norm = threshold * DBL_EPSILON;
   max = threshold + 1.0;
   while (threshold > threshold_norm) {
      threshold /= 10.0;
      if (max < threshold) continue;
      max = 0.0;
      for (pAk = A, k = 0; k < (n-1); pAk += n, k++) {
         for (pAm = pAk + n, m = k + 1; m < n; pAm += n, m++) {
            if ( fabs(*(pAk + m)) < threshold ) continue;

                 // Calculate the sin and cos of the rotation angle which
                 // annihilates A[k][m].

            cot_2phi = 0.5 * ( *(pAk + k) - *(pAm + m) ) / *(pAk + m);
            dum1 = sqrt( cot_2phi * cot_2phi + 1.0);
            if (cot_2phi < 0.0) dum1 = -dum1;
            tan_phi = -cot_2phi + dum1;
            tan2_phi = tan_phi * tan_phi;
            sin2_phi = tan2_phi / (1.0 + tan2_phi);
            cos2_phi = 1.0 - sin2_phi;
            sin_phi = sqrt(sin2_phi);
            if (tan_phi < 0.0) sin_phi = - sin_phi;
            cos_phi = sqrt(cos2_phi); 
            sin_2phi = 2.0 * sin_phi * cos_phi;
            cos_2phi = cos2_phi - sin2_phi;

                     // Rotate columns k and m for both the matrix A 
                     //     and the matrix of eigenvectors.

            p_r = A;
            dum1 = *(pAk + k);
            dum2 = *(pAm + m);
            dum3 = *(pAk + m);
            *(pAk + k) = dum1 * cos2_phi + dum2 * sin2_phi + dum3 * sin_2phi;
            *(pAm + m) = dum1 * sin2_phi + dum2 * cos2_phi - dum3 * sin_2phi;
            *(pAk + m) = 0.0;
            *(pAm + k) = 0.0;
            for (i = 0; i < n; p_r += n, i++) {
               if ( (i == k) || (i == m) ) continue;
               if ( i < k ) dum1 = *(p_r + k); else dum1 = *(pAk + i);
               if ( i < m ) dum2 = *(p_r + m); else dum2 = *(pAm + i);
               dum3 = dum1 * cos_phi + dum2 * sin_phi;
               if ( i < k ) *(p_r + k) = dum3; else *(pAk + i) = dum3;
               dum3 = - dum1 * sin_phi + dum2 * cos_phi;
               if ( i < m ) *(p_r + m) = dum3; else *(pAm + i) = dum3;
            }
            for (p_e = eigenvectors, i = 0; i < n; p_e += n, i++) {
               dum1 = *(p_e + k);
               dum2 = *(p_e + m);
               *(p_e + k) = dum1 * cos_phi + dum2 * sin_phi;
               *(p_e + m) = - dum1 * sin_phi + dum2 * cos_phi;
            }
         }
         for (i = 0; i < n; i++)
            if ( i == k ) continue;
            else if ( max < fabs(*(pAk + i))) max = fabs(*(pAk + i));
      }
   }
   for (pAk = A, k = 0; k < n; pAk += n, k++) eigenvalues[k] = *(pAk + k); 
}


 void SWD_Maths::ComputFirstMoments(double& fVector, float matrix[ImHeight][ImWidth]
		,  int hFlag, int wFlag, int nLevel)
 {
   int i,j;
   int mNum = pow(2.0,nLevel);
   int mNumW = mNum * 2;

   for(i = (mNum * hFlag) ; i< (mNum * hFlag)+mNum; i++)
		for(j = (mNumW * wFlag )  ; j<(mNumW * wFlag )+mNumW ; j++)
		{          
          fVector += fabs(matrix[i][j]);
		}

	fVector = fVector/(mNum*mNumW) ;
 }

void SWD_Maths::ComputSecondMoments(double& fVector, double fVectorfirst, float matrix[ImHeight][ImWidth]
		, int hFlag, int wFlag, int nLevel)
{
 double temp = 0;
 int i,j;
 int mNum = pow(2.0,nLevel);
 int mNumW = mNum * 2;

	for(i = (mNum * hFlag) ; i<(mNum * hFlag)+mNum; i++)
		for(j = (mNumW * wFlag )  ; j< (mNumW * wFlag)+mNumW ; j++)
		{ 
			 temp = (fabs(matrix[i][j]) - fVectorfirst)*
                            (fabs(matrix[i][j]) - fVectorfirst);
          fVector += temp;
		}

	fVector = fVector / (mNum*mNumW-1) ;
	fVector = sqrt(fVector) ;
}
 void SWD_Maths::ComputThirdMoments(double& fVector, double fVectorfirst, float matrix[ImHeight][ImWidth]
		, int hFlag, int wFlag, int nLevel)
 {
   double f1,f2;
	f1 = 0 ;
	f2 = 0 ;
	double temp = 0;
    int i,j;
    int mNum = pow(2.0,nLevel);
    int mNumW = mNum * 2;
 
	for(i = (mNum * hFlag) ; i<(mNum * hFlag)+mNum; i++)
		for(j = (mNumW * wFlag )  ; j< (mNumW * wFlag)+mNumW ; j++)
		{ 
          temp = pow(fabs(matrix[i][j]) - fVectorfirst, 3);
          f1 += temp;

		  temp = pow(fabs(matrix[i][j]) - fVectorfirst, 2);
		  f2 += temp;

		}
	
	f2 = pow(f2, 3);
	f2 = sqrt(f2);

	fVector = f1/f2;

	fVector = fVector*(mNum*mNumW)*sqrt((double)(mNum*mNumW - 1))
		   /(mNum*mNumW - 2);
 }


 /*----------------------------------------------------------------
 *\ brief  about SWT 1 ,2 ,3 MOMENTS
 *\author  李欧州
 *\date    10-8-2008
 *----------------------------------------------------------*/
 void SWD_Maths::ComputSWTFirstMoments(vector<double>& vect, vector<double>& descriptorVect)
 {   
	 double temp = 0;
	 for(int i = 0; i<vect.size(); i++)
	 {
       temp += fabs(vect[i]);
	 }
	  
	  temp = temp/vect.size();

	  descriptorVect.push_back(temp);

 }

 void SWD_Maths::ComputSWTSecondMoments(vector<double>& vect, vector<double>& descriptorVect)
 {   
	 double temp = 0;
	 int vsize = vect.size();
	 double dt = *(descriptorVect.rbegin());
	 for(int i = 0; i<vsize; i++)
	 {
       double secondtemp = fabs(vect[i]) - dt;
       temp += (secondtemp*secondtemp);
	 }
	  
	  temp = temp/(vsize-1);
	  temp = sqrt(temp);

	  descriptorVect.push_back(temp);

 }

 void SWD_Maths::ComputSWTThirdMoments(vector<double>& vect, vector<double>& descriptorVect)
 {   
	 double temp = 0;
	 double f1 = 0;
	 double f2 = 0;
	 int vsize = vect.size();
	 double dt = *(descriptorVect.rbegin()+1);
	 for(int i = 0; i<vsize; i++)
	 {
       double secondtemp = fabs(vect[i]) - dt;
       f1 += pow(secondtemp,3);
	   f2 += pow(secondtemp,2);
	 }
	  
	 f2 = pow(f2, 3);
	 f2 = sqrt(f2);
	 
	 f1 = f1/f2;
	 temp = f1 * vsize * sqrt((double)(vsize-1))/(vsize-2);
	  descriptorVect.push_back(temp);

 }
