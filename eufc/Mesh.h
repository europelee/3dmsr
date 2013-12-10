// Mesh.h: interface for the Mesh class.
//
/*--------------------------------------------------------------
*
*\brief   Mesh(3d model, i.e 3d triangle mesh  )
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

#ifndef MESH_H
#define MESH_H
class MeshTrianSph  
{
public:
	float m_minLong ;
	float m_minLat ;
	float m_maxLong ;
	float m_maxLat ;
	bool m_crossX;
	MeshTrianSph()
	{
	 m_crossX = false;
	 m_minLong = 0;
	 m_minLat = 0;
	 m_maxLong = 0;
	 m_maxLat = 0;
	}
	virtual ~MeshTrianSph(){}

};
class Vertex3d;
class Face;

class Mesh  
{
public:
	int m_nverts;
	int m_nfaces;
	Vertex3d *m_verts;
	Face *m_faces;
	void MeshFaceNormal();
	Mesh* MeshClone();
	Mesh();
	~Mesh();

};

#endif 
