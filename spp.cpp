#include "RinexNavRead.h"
#include "RinexObsRead.h"
#include "stdlib.h"
#include "ComputeSatPosition.h"
#include "DataManage.h"
#include <iostream>

using namespace std;

int main()
{
 char strn[100],stro[100];
cout<<"###################################################"<<endl;
cout<<"###                                             ###"<<endl;
cout<<"###             ���㶨λ���� V1.0               ###"<<endl;
cout<<"###                                             ###"<<endl;
cout<<"###################################################"<<endl;
//cout<<"������N�ļ�����";
//gets(strn);
//cout<<"������O�ļ�����";
//gets(stro);

string path("AVCA3000.13o");
GMO worun;
worun=ReadRinexObsFile(path);
string pathnav("brdc3000.13n");
vector<GMNREC>  nav;
 nav=ReadRinexNavFile(pathnav);
vector<PPONERESULT>  all;
all=PP(worun,nav);
 ResultOut(all);
	return 0;
}

