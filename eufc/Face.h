// Face.h: interface for the Face class.
//
/*--------------------------------------------------------------
*
*\brief   Face(only for triangle face )
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

#ifndef FACE_H
#define FACE_H

class Vertex3d;

class Face  
{
public:
	int m_nverts;
  Vertex3d **m_verts;
  float m_normal[3];
  int m_verIndex[3];
  //static int objectCount;
  Face();	
 virtual ~Face();


};

#endif 
