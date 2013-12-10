// ModelRead.cpp: implementation of the ModelRead class.
//
//////////////////////////////////////////////////////////////////////

#include "myconfig.h"
#include "ModelRead.h"
#include "Mesh.h"
#include "Face.h"
#include "Vertex3d.h"



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ModelRead::ModelRead()
{
}

ModelRead::~ModelRead()
{

}
Mesh* ModelRead::ReadOffFile(const char *filename)
{
	FILE *fp;
	
  if (!(fp = fopen(filename, "r"))) 
  {
    //fprintf(stderr, "Unable to open file %s\n", filename);
	  //AfxMessageBox("Unable to open 3d model file");
    return NULL;
  }

  Mesh *mesh = new Mesh();
  if (!mesh) 
  {
    //fprintf(stderr, "Unable to allocate memory for file %s\n", filename);
	   //AfxMessageBox("Unable to allocate memory for file");
    fclose(fp);
    return NULL;
  }

  // Read file
  int nverts = 0;
  int nfaces = 0;
  int nedges = 0;
  int line_count = 0;
  char buffer[1024];

  while (fgets(buffer, 1023, fp))
  {
	  
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) 
		bufferp++;

	// Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Check section
    if (nverts == 0) 
	{
      // Read header 
      if (!strstr(bufferp, "OFF")) 
	  {
        // Read mesh counts
        if ((sscanf(bufferp, "%d%d%d", &nverts, &nfaces, &nedges) != 3) || (nverts == 0)) 
		{
          //fprintf(stderr, "Syntax error reading header on line %d in file %s\n", line_count, filename);
            //AfxMessageBox("Syntax error reading header on line ");
			fclose(fp);
          return NULL;
		}
		
		// Allocate memory for mesh
        mesh->m_verts = new Vertex3d [nverts];
//        ASSERT(mesh->m_verts);
        mesh->m_faces = new Face [nfaces];
   //     ASSERT(mesh->m_faces);
	  }
     }
	  else 
		  if (mesh->m_nverts < nverts)
		  {
              // Read vertex coordinates
            Vertex3d& vert = mesh->m_verts[mesh->m_nverts++];
	  
            if (sscanf(bufferp, "%f%f%f", &(vert.m_x), &(vert.m_y), &(vert.m_z)) != 3) 
			{
              //fprintf(stderr, "Syntax error with vertex coordinates on line %d in file %s\n", line_count, filename);
               //AfxMessageBox("Syntax error with vertex coordinates on line");
				fclose(fp);
               return NULL;
			}
			int aaa = 1;
		  }

		  else 
			  if (mesh->m_nfaces < nfaces) 
			  {
                 // Get next face
                Face& face = mesh->m_faces[mesh->m_nfaces++];

                // Read number of vertices in face 
                bufferp = strtok(bufferp, " \t");
                if (bufferp) 
					face.m_nverts = atoi(bufferp);
                else 
				{
                  //fprintf(stderr, "Syntax error with face on line %d in file %s\n", line_count, filename);
                   //AfxMessageBox("Syntax error with face on line");
					fclose(fp);
                  return NULL;
				}

				 // Allocate memory for face vertices
                 face.m_verts = new Vertex3d *[face.m_nverts];
//                 ASSERT(face.m_verts);

                 // Read vertex indices for face
                 for (int i = 0; i < face.m_nverts; i++) 
				 {
                    bufferp = strtok(NULL, " \t");//
                    if (bufferp) 
					{
						face.m_verts[i] = &(mesh->m_verts[atoi(bufferp)]);
						face.m_verIndex[i] = atoi(bufferp);
					}
                    else 
					{
                     //fprintf(stderr, "Syntax error with face on line %d in file %s\n", line_count, filename);
                     // AfxMessageBox("Syntax error with face on line");
						fclose(fp);
                     return NULL;
					}
				 }
				 
				  // Compute normal for face
                   face.m_normal[0] = face.m_normal[1] = face.m_normal[2] = 0;
///////////////////////////my code////////////////////////////////////////////////////////
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
	  /////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Vertex3d *v1 = face.m_verts[face.m_nverts-1];
//
//                   for (i = 0; i < face.m_nverts; i++) 
//				   {
//                    Vertex3d *v2 = face.m_verts[i];
//                    face.m_normal[0] += (v1->m_y - v2->m_y) * (v1->m_z + v2->m_z);
//                    face.m_normal[1] += (v1->m_z - v2->m_z) * (v1->m_x + v2->m_x);
//                    face.m_normal[2] += (v1->m_x - v2->m_x) * (v1->m_y + v2->m_y);
//                    v1 = v2;
//				   }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
                  // Normalize normal
                  float squared_normal_length = 0.0;
                  squared_normal_length += face.m_normal[0]*face.m_normal[0];
                  squared_normal_length += face.m_normal[1]*face.m_normal[1];
                  squared_normal_length += face.m_normal[2]*face.m_normal[2];
                  float normal_length = sqrt(squared_normal_length);
				  
                  if (normal_length > 1.0E-6) 
				  {
                    face.m_normal[0] /= normal_length;
                    face.m_normal[1] /= normal_length;
                    face.m_normal[2] /= normal_length;
				  }
			  }
              else 
			  {
               // Should never get here
                 //fprintf(stderr, "Found extra text starting at line %d in file %s\n", line_count, filename);
                  //AfxMessageBox("Found extra text starting at line");
				  break;
			  }
		
} //while (fgets(buffer, 1023, fp))
	// Check whether read all faces
  if (nfaces != mesh->m_nfaces)
  {
    // AfxMessageBox("Expected faces, but read only faces in file");
	  //fprintf(stderr, "Expected %d faces, but read only %d faces in file %s\n", nfaces, mesh->m_nfaces, filename);
     return NULL; 
  }

  // Close file
  fclose(fp);

  // Return mesh 
  return mesh;			

}


Mesh* ModelRead::ReadObjFile(const char *filename)
{
   //Mesh *mesh = new Mesh();
   //vtkOBJReader *pvtkObjReader = vtkOBJReader::New();
   //pvtkObjReader->SetFileName(filename);
   //pvtkObjReader->Update();
  
   //vtkPolyData* temp = pvtkObjReader->GetOutput();

   //mesh->m_nfaces = temp->GetNumberOfPolys();
   //mesh->m_nverts = temp->GetNumberOfPoints();
   //mesh->m_verts = new Vertex3d [mesh->m_nverts];
   //mesh->m_faces = new Face [mesh->m_nfaces];

   //double x2[3];int i;
   //for(i=0; i<mesh->m_nverts; i++)
   //{
	  // temp->GetPoint(i,x2);
	  // mesh->m_verts[i].m_x = x2[0];
	  // mesh->m_verts[i].m_y = x2[1];
	  // mesh->m_verts[i].m_z = x2[2];
   //}

   //vtkIdType vnum_face;
   //int*  pveIndexface ;//= new int[3];不必为其申请内存，看getnextcell源码，即清楚
   ////vtkIdType * pvI = veIndex_face;
   //vtkCellArray * pPoly = temp->GetPolys();


   //for(i=0; i<mesh->m_nfaces; i++)
   //{
	  // int flag = pPoly->GetNextCell(vnum_face,pveIndexface);

	  // if(flag == 0)
	  // {
		 //  CString str ;
		 //  str.Format("error for dec  mesh %d,%d",i,flag);
		 //  AfxMessageBox(str);
		 //  delete mesh;
		 //  pvtkObjReader->Delete();
		 //  pvtkObjReader = NULL;
		 //  return NULL;

	  // }
	  // if(vnum_face != 3 )
	  // {   
		 //  CString str ;
		 //  str.Format("error for dec not a triangle mesh %d,%d",i,vnum_face);
		 //  AfxMessageBox(str);
   //        delete mesh;
		 //  pvtkObjReader->Delete();
		 //  pvtkObjReader = NULL;
		 //  return NULL;
	  // }
	  // mesh->m_faces[i].m_nverts = vnum_face;
	  // mesh->m_faces[i].m_verts = new Vertex3d *[vnum_face];

	  // for(int j = 0 ; j<3 ;j++)
	  // {
		 //  vtkIdType index = pveIndexface[j];
		 //  mesh->m_faces[i].m_verts[j] = &(mesh->m_verts[index]);
   //        mesh->m_faces[i].m_verIndex[j] = index;
	  // }
   //}

   //pvtkObjReader->Delete();
   //pvtkObjReader = NULL;

   //return mesh;
	return NULL;
}
Mesh* ModelRead::Read3DSFile(const char *filename)
{

	return NULL;
}