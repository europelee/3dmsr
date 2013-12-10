/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkOFFReader.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOFFReader - read  .off files
// .SECTION Description
// vtkOFFReader is a source object that reads .off
// files. The output of this source object is polygonal data.

//vtkOFFReader from vtkOBJReader
//by europelee 
//2009-7-27

#ifndef __vtkOFFReader_h
#define __vtkOFFReader_h

#ifdef EUVTKOFF_DLLAPI
#else
#define EUVTKOFF_DLLAPI _declspec(dllimport)
#endif

//class vtkPolyDataAlgorithm;
#include <vtkPolyDataAlgorithm.h>
class EUVTKOFF_DLLAPI vtkOFFReader : public vtkPolyDataAlgorithm 
{
public:
  static vtkOFFReader *New();
  vtkTypeRevisionMacro(vtkOFFReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of .off file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  vtkOFFReader();
  ~vtkOFFReader();
  
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  char *FileName;
private:
  vtkOFFReader(const vtkOFFReader&);  // Not implemented.
  void operator=(const vtkOFFReader&);  // Not implemented.
};

#endif
