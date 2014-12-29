// SpheCoord.h: interface for the SpheCoord class.
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
/*-----------------------------------------------------------------
*brief  modify for 3DQModelRetrieval,add  a constructor with 4 params
*author 李欧州
*date   21-02-2009 
-----------------------------------------------------------------*/
//////////////////////////////////////////////////////////////////////

#ifndef SPHECOORD_H
#define SPHECOORD_H

class SpheCoord  
{
public:
	float m_lat;
	float m_long;
	int m_rowIndex;
	int m_colIndex;
	bool operator < (const SpheCoord &m)const;
	SpheCoord();
	SpheCoord(float f1,float f2,int r,int c);
	virtual ~SpheCoord();

};

#endif 
