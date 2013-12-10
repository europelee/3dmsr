// Vertex3d.h: interface for the Vertex3d class.
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
//history
/*-----------------------------------------------------------------
*brief  modify for 3DQModelRetrieval
*author 李欧州
*date   17-02-2009 
-----------------------------------------------------------------*/
//////////////////////////////////////////////////////////////////////

/*#if !defined(AFX_VERTEX3D_H__01C253AB_6DFE_4B63_9D91_A3949EA11C66__INCLUDED_)
#define AFX_VERTEX3D_H__01C253AB_6DFE_4B63_9D91_A3949EA11C66__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
*/
#ifndef VERTEX3D_H
#define VERTEX3D_H

class Vertex3d  
{
public:
	float m_x;
	float m_y;
	float m_z;
	Vertex3d();
	Vertex3d(float c1, float c2, float c3);


	Vertex3d* operator- (Vertex3d sub);
	Vertex3d* operator+ (Vertex3d add);
	float operator* (Vertex3d v);//inner product
	Vertex3d operator^ (Vertex3d v);// cross
	float vectLength() {
		return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
	}
	void VertNormal();
	virtual ~Vertex3d();

};
#endif
//#endif // !defined(AFX_VERTEX3D_H__01C253AB_6DFE_4B63_9D91_A3949EA11C66__INCLUDED_)
