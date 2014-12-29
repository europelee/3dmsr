// SpheCoord.cpp: implementation of the SpheCoord class.
//
/*--------------------------------------------------------------
*
*\brief   about spherical coordinate
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
#include "SpheCoord.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


SpheCoord::	SpheCoord()
{m_lat=0;m_long=0;m_rowIndex=0;m_colIndex=0;}
SpheCoord::SpheCoord(float f1,float f2,int r,int c)
{
  m_lat=f1;m_long=f2;m_rowIndex=r;m_colIndex=c;  
}    
bool SpheCoord::operator < (const SpheCoord &m)const {
		  if (m_long < m.m_long) 
		         return true;
		  else
			  if(m_long == m.m_long && m_lat < m.m_lat)
                return true;
			  else
				  return false;
				  
        }
SpheCoord::~SpheCoord()
{

}
