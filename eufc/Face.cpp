// Face.cpp: implementation of the Face class.
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

#include "myconfig.h"
#include "Face.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Face::Face()
  {   
	 // objectCount++;
	  m_nverts=0;
	  m_normal[0]=0;
	  m_normal[1]=0;
	  m_normal[2]=0;
	  m_verIndex[0] = 0;
	   m_verIndex[1] = 0;
	    m_verIndex[2] = 0;
  }

Face::~Face()
{
   delete []m_verts;
	  m_verts = NULL;
	  
}
