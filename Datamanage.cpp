#include "DataManage.h"

void ResultOut(vector<PPONERESULT> m_result)
{
	FILE *fp;
	if((fp=fopen("Result.dat","w"))==NULL)
	{
	    printf("ʧ��\n");
	    exit(0);
	}
    
	char dis1[100]="#���㶨λ�������ļ�\n";
    char year[10]="#�� ,";
	char month[10]="��,";
	char day[10]="��,";
	char hour[10]="ʱ,";
	char minute[10]="��,";
	char second[10]="��, ";

	char X[20]="          X(m)   ";
	char Y[20]="      Y(m)     ";
	char Z[20]="     Z(m)     ";
	char T[20]="clk_bias(ns)   ";

	char SX[20]="     sigmaX(m)   ";
	char SY[20]="   sigmaY(m)   ";
	char SZ[20]="    sigmaZ(m)    ";

	char PDOP[10]=" PDOP  ";
	char num[10]="������";
	char hang[2]="\n"; 
	 
	fputs(dis1,fp);
	//fputs(dis2,fp);

	fputs(year,fp);
	fputs(month,fp);
	fputs(day,fp);
	fputs(hour,fp);
	fputs(minute,fp);
	fputs(second,fp);
	fputs(X,fp);fputs(SX,fp);
	fputs(Y,fp);fputs(SY,fp);
	fputs(Z,fp);fputs(SZ,fp);
	fputs(T,fp);
	fputs(PDOP,fp);fputs(num,fp);
	fputs(hang,fp);

	int len=m_result.size();
	cout<<"�ļ���д��..."<<endl;

	for(int i=0;i<len;i++)
	{
		fprintf(fp,"%4d,%2d,%2d,%2d,%2d,%6.3lf,",m_result[i].epoch.year,
			m_result[i].epoch.month,m_result[i].epoch.day,m_result[i].epoch.hour,
			m_result[i].epoch.minute,m_result[i].epoch.second);
		fprintf(fp,"%14.3lf,%14.10lf,%14.3lf,%14.10lf,%14.3lf,%14.10lf, %14.11lf",
			m_result[i].crd.x,m_result[i].crdsigma.x,
			m_result[i].crd.y,m_result[i].crdsigma.y,
			m_result[i].crd.z,m_result[i].crdsigma.z,
			m_result[i].clk_bias);
		fprintf(fp,"%8.3lf",m_result[i].PDOP);
		fprintf(fp,"%3d\n",m_result[i].sat_num);
       
	
	}
	

	 cout<<"�ļ��������!"<<endl;
	fclose(fp);

}