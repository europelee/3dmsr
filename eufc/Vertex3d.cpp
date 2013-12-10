// Vertex3d.cpp: implementation of the Vertex3d class.
//
/*--------------------------------------------------------------
*
*\brief   Vertex3d(also as vector3d)
*\author   李欧州
*\date     20-9-2008
*
*
*
------------------------------------------------------------------*/
//////////////////////////////////////////////////////////////////////
//history
/*-----------------------------------------------------------------
*brief  modify for 3DQModelRetrieval
*author 李欧州
*date   17-02-2009 
-----------------------------------------------------------------*/
#include "myconfig.h"
//#include "3DMDggpile.h"
#include "Vertex3d.h"
/*
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
*/
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
Vertex3d::Vertex3d()
{m_x=0;m_y=0;m_z=0;}

Vertex3d::Vertex3d(float c1, float c2, float c3)
	{m_x=c1;m_y=c2;m_z=c3;}

void Vertex3d::VertNormal()
{   
	float normv = 0;
	normv = sqrt((this->m_x)*(this->m_x)+(this->m_y)*(this->m_y)
				 +(this->m_z)*(this->m_z));
            this->m_x = this->m_x/normv;
			this->m_y = this->m_y/normv;
			 this->m_z = this->m_z/normv;
}

                        // required for DBL_EPSILON
Vertex3d* Vertex3d::operator+ (Vertex3d add)
{
  Vertex3d* val2 = new Vertex3d();

	val2->m_x = this->m_x + add.m_x;
	val2->m_y = this->m_y + add.m_y;
	val2->m_z = this->m_z + add.m_z;

	return val2;
}
Vertex3d* Vertex3d::operator- (Vertex3d sub)
{
    Vertex3d* val = new Vertex3d();

	val->m_x = this->m_x - sub.m_x;
	val->m_y = this->m_y - sub.m_y;
	val->m_z = this->m_z - sub.m_z;

	return val;
}

float Vertex3d::operator* (Vertex3d v)
{
    
	return (this->m_x * v.m_x + this->m_y * v.m_y + this->m_z * v.m_z);
}

Vertex3d Vertex3d::operator^ (Vertex3d v)
{   

	return Vertex3d(this->m_y * v.m_z - this->m_z * v.m_y,
         this->m_z * v.m_x - this->m_x * v.m_z,
		 this->m_x * v.m_y - this->m_y * v.m_x
		);
}

Vertex3d::~Vertex3d()
{

}
