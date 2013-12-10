// Mesh.cpp: implementation of the Mesh class.
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

#include "myconfig.h"
#include "Mesh.h"
#include "Face.h"
#include "Vertex3d.h"
#include "SWD_Maths.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Mesh::Mesh()
{

m_nfaces=0;
m_nverts=0;


}

Mesh::~Mesh()
{

delete []m_verts; 
m_verts = NULL; 
m_nverts = 0;
delete []m_faces; 
m_faces = NULL;
m_nfaces = 0;

}

void Mesh::MeshFaceNormal()
{
  	for(int m=0; m < m_nfaces; m++)
	{  
	   Face& face = m_faces[m];
	   SWD_Maths::ComputFaceNormal(face);
	}

}
Mesh* Mesh::MeshClone()
{
   Mesh* m_meshBack = new Mesh();
   m_meshBack->m_nfaces = this->m_nfaces;
   m_meshBack->m_nverts = this->m_nverts;
   m_meshBack->m_verts = new Vertex3d [m_meshBack->m_nverts];
   m_meshBack->m_faces = new Face [m_meshBack->m_nfaces];

   int i;
   for(i=0; i<this->m_nverts; i++)
   {
	  m_meshBack->m_verts[i].m_x = this->m_verts[i].m_x;
	  m_meshBack->m_verts[i].m_y = this->m_verts[i].m_y;
	  m_meshBack->m_verts[i].m_z = this->m_verts[i].m_z;
   }
   for(i=0; i<this->m_nfaces; i++)
   {
	   m_meshBack->m_faces[i].m_nverts = 3;
	   m_meshBack->m_faces[i].m_verts = new Vertex3d *[3];

	   for(int j = 0 ; j<3 ;j++)
	   {
		   //vtkIdType index = pveIndexface[j];
		   //mesh_dec->m_faces[i].m_verts[j] = &(mesh_dec->m_verts[index]);
		   int index = this->m_faces[i].m_verIndex[j];
            m_meshBack->m_faces[i].m_verts[j] = &(m_meshBack->m_verts[index]);
			m_meshBack->m_faces[i].m_verIndex[j] = index;
	   }
   }

   return m_meshBack;

}
