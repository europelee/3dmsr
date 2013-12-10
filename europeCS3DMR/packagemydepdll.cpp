/********************************************************
*brief: package my depmapdll into a .exe file
*author: europelee
*date: 2009-11-22 16:53
*********************************************************/
#include <afx.h> 
#include<stdio.h>
#include <stdlib.h>
#include "DepthMapdll.h"

void SaveFeature2File(const char *fileName,float *vranicfft,int vsize)
{
	FILE* pff;
	pff = fopen(fileName,"ab");
	//for(int i=0;i<vsize;i++)//438 need insteaded
	//{
	//	fprintf(pff,"%.10f\n",vranicfft[i]);
	//}
	fwrite(vranicfft,sizeof(float),vsize,pff);
	fclose(pff);
}
//void TestWriteFile(const char *fileName)
//{
//  FILE* pff;
//  pff = fopen(fileName,"rb");
//  float temp[438]={0};
//  fread(temp,sizeof(float),438,pff);
//  for(int i=0;i<438;i++)
//	  printf("%.10f  * ",temp[i]);
//}
void AllFeatureExt(CString dirPath,const char *fileName)
{
	CFileFind finder;
	CString strWildCard(dirPath);
	strWildCard+=_T("\\*.*");
	BOOL bWorking=finder.FindFile(strWildCard);

	while(bWorking)
	{
		bWorking=finder.FindNextFile();

		if(finder.IsDots() /*|| finder.IsDirectory()*/)
			continue; 
		else
			if(finder.IsDirectory())
			{
				CString filePath = finder.GetFilePath();
				AllFeatureExt(filePath,fileName);
			}
			//continue;
			else
			{
				CString filePath = finder.GetFilePath();
				char temp[500];
				memset(temp,0,500*sizeof(char));
				USES_CONVERSION;
				strcpy(temp,(char*)(W2A((LPCTSTR)filePath)));
				// strncpy(temp,filePath,filePath.GetLength());
				
				EU_DepthMapp euobj;
				euobj.ComputeDepMap(temp);
				SaveFeature2File(fileName,euobj.vranicfft,438);
			}
	}//while(bWorking)
   finder.Close();
}
int main(int argc, char *argv[])
{
     char *InDirectory, *OutFile;  
    if (argc < 3) {
    printf("\n\
3d model feature extraction for 3d shape retrieval.\n\
            by Europe Lee.\n\
          in November 2009.\n\
Usage: \n\
  %s InDirectory OutFile \n\
    InDirectory    : only support .off file\n\
    OutFile          : output feature into the file\n\
",argv[0]);
   // system("pause");
   exit(1);
  }
  
  
  InDirectory = argv[1];
  OutFile = argv[2];
  CString str = InDirectory;
  
  AllFeatureExt(str,OutFile);
  //TestWriteFile(OutFile);
    return 0;
}    
