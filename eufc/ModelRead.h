// ModelRead.h: interface for the ModelRead class.
//
/*-------------------------------------------------------------------
*
*\brief   about some kinds of 3d model file read
*\author  ¿Ó≈∑÷›
*\date    22-9-2008
*
*
-----------------------------------------------------*/
//////////////////////////////////////////////////////////////////////

#ifndef MODELREAD_H
#define MODELREAD_H

class Mesh;

class ModelRead  
{
public:
	static Mesh* ReadOffFile(const char *filename);
	static Mesh* ReadObjFile(const char *filename);
	static Mesh* Read3DSFile(const char *filename);
	ModelRead();
	virtual ~ModelRead();

};

#endif // !defined(AFX_MODELREAD_H__643CB214_A9B7_48EF_B75A_32E4D20D51FD__INCLUDED_)
