#include "ComputeSatPosition.h"
#include "math.h"




GMNREC GetBestGMNREC(vector<GMNREC> navRecord,
				   int nPRN,PCOMMONTIME pctEpoch)
{
    bool flag=false;
	bool satPRN=false;
	int n=0;
	
	JULIANDAY sat_toe,sat_epoch;
    
	CommonTimeToJulianDay (pctEpoch,&sat_epoch);

	for(int i=0;i<(int)navRecord.size();i++)
	{
		if(navRecord[i].PRN==nPRN)
		{
			satPRN=true;
			GPSTimeToJulianDay (&navRecord[i].TOE, &sat_toe);

			if(GetTimeDelta (&sat_epoch,&sat_toe)<0)
			{
				return navRecord[i];
                flag=true;
				break;
			}
			else if(fabs(GetTimeDelta (&sat_epoch,&sat_toe))*_DAY_IN_SECOND<60*60)
			{
				return navRecord[i];
				flag=true;
				break;
				
			}
			n=i;

		}
		
	
	}
	if(!satPRN)
	{
		cout<<"the satPRN didn't exsit!"<<endl;
      
	}
	if(!flag)
         return navRecord[n];

	

}

double EofMe(double M,double e,double tol)
{
	double E0,E1;


	E0=M;
	
	E1=E0-(E0-e*sin(E0)-M)/(1-e*cos(E0));
	
	while(fabs(E1-E0)>tol)
	{
		E0=E1;
		//cout<<"test  1"<<E0<<endl;
		E1=E0-(E0-e*sin(E0)-M)/(1-e*cos(E0));
		//cout<<"M="<<M<<endl;

		//E1=M+e*sin(E0);
		//cout<<"test  2"<<E1<<endl;
		//cout<<"test  3"<<E1-E0<<endl;
     }
	return E1;

}

double Get_atan(double z,double y)
{
   double x;
   if (z==0) x=PI/2;
   else{
	if (y==0) x=PI;
	else{
	      x=atan(fabs(y/z));
	      if ((y>0)&&(z<0)) x=PI-x;
	      else if ((y<0)&&(z<0)) x=PI+x;
		   else if((y<0)&&(z>0)) x=2*PI-x;
	     }
       }
   return x;
}


void GetUtilParameter(vector<GMNREC> navRecord,
				   int nPRN,PCOMMONTIME pctEpoch,PUtilParam pParam)
{
	GMNREC  theBestGMN;

	theBestGMN=GetBestGMNREC(navRecord,nPRN,pctEpoch);

	//��������ƽ�����ٶ�
	double n0=sqrt(GM)/ pow(theBestGMN.SqrtA,3);//3-11

	//��������������ο���Ԫ��ʱ��

    JULIANDAY toe,epoch;

	

	CommonTimeToJulianDay (pctEpoch,&epoch);


	GPSTimeToJulianDay (&theBestGMN.TOE, &toe);




	double tk;
	tk=GetTimeDelta (&epoch,&toe)*_DAY_IN_SECOND;
	if(tk>302400)
		tk-=604800;
	else if(tk<-302400)
		tk+=604800;
	else
		tk=tk;

	//��ƽ�����ٶȽ��и���
	double n=n0+theBestGMN.deltn;//3-12

	(*pParam).n=n;

	//����ƽ�����
	double M=theBestGMN.M0+n*tk;//3-13

	//��ƫ�����
    double E;
	E=EofMe(M,theBestGMN.e,1e-10);//3-14

	(*pParam).E=E;

	//����������
	double vk,cosvk,sinvk;
	cosvk=(cos(E)-theBestGMN.e)/(1-theBestGMN.e*cos(E));
	sinvk=(sqrt(1-theBestGMN.e*theBestGMN.e)*sin(E))/(1-theBestGMN.e*cos(E));

//	vk=atan2(sinvk,cosvk);
//	if(vk<0)
//	 vk+=2*PI;

	vk=Get_atan(cosvk,sinvk);//3-17
	 
	(*pParam).vk=vk;
	 //���������Ǿ�
	 double u0;
	 u0=theBestGMN.omiga+vk;
     (*pParam).u0=u0;

	 //�����㶯������;3-18
	     //1.���������Ǿ�ĸ�����
	 double detU=theBestGMN.Cus*sin(2*u0)+theBestGMN.Cuc*cos(2*u0);
	     //2.�����򾶵ĸ�����
	 double detR=theBestGMN.Crs*sin(2*u0)+theBestGMN.Crc*cos(2*u0);
	     //3.��������Ǹ�����
	 double detI=theBestGMN.Cis*sin(2*u0)+theBestGMN.Cic*cos(2*u0);


	 //������������Ǿ�;3-19
	 double uk=u0+detU;
	 (*pParam).uk=uk;

	 //���㾭����������
	 double r=theBestGMN.SqrtA*theBestGMN.SqrtA*(1-theBestGMN.e*cos(E))+detR;
	 (*pParam).r=r;

	 //���㾭�������Ĺ�����
	 double i=theBestGMN.i0+detI+theBestGMN.iDot*tk;
	 (*pParam).i=i;


	 //����˲ʱ�����㾭��;3-23
	 double L=theBestGMN.omiga0+(theBestGMN.omigaDot-we)*tk
		 -we*(theBestGMN.TOE.tow.sn+theBestGMN.TOE.tow.tos);
	 (*pParam).L=L;



}

void GetOrbNClk(vector<GMNREC> navRecord,int nPRN,
PCOMMONTIME pctEpoch, PCRDCARTESIAN pcrdOrb)
{
    
    UtilParam pParam;
    GetUtilParameter(navRecord,nPRN,pctEpoch,&pParam);

	double n=pParam.n;
	double E=pParam.E;
	double vk=pParam.vk;
	double uk=pParam.uk;
	double r=pParam.r;
	double i=pParam.i;
	double L=pParam.L;
	double u0=pParam.u0;

	//���������ڹ��ƽ���ϵ�λ��;3-20
	 double x1=r*cos(uk);
	 double y1=r*sin(uk);
	 double z1=0;

	 //�����ڵع�����ϵ�µ�λ��;3-24
	 pcrdOrb->x=cos(L)*x1-cos(i)*sin(L)*y1;
	 pcrdOrb->y=sin(L)*x1+cos(i)*cos(L)*y1;
	 pcrdOrb->z=sin(i)*y1;

	 
}

void GetSVClkBias(vector<GMNREC> navRecord,int nPRN,
PCOMMONTIME pctEpoch,double* pdSVClkBias,double *detj)
{
	UtilParam pParam;
    GetUtilParameter(navRecord,nPRN,pctEpoch,&pParam);

	GMNREC  theBestGMN;

	theBestGMN=GetBestGMNREC(navRecord,nPRN,pctEpoch);

	double E=pParam.E;
	
    JULIANDAY/* toe,*/epoch;
	CommonTimeToJulianDay (pctEpoch,&epoch);
	//GPSTimeToJulianDay (&theBestGMN.TOE, &toe);

	 //����C/A���źŷ���ʱ�̵ĸ���
	 double dettr=F*theBestGMN.e*theBestGMN.SqrtA*sin(E);  //4-14

	 JULIANDAY toc;
	 GPSTimeToJulianDay (&theBestGMN.TOC, &toc);

	 double dettoc=GetTimeDelta (&epoch,&toc)*_DAY_IN_SECOND;
     //4-21
	 *pdSVClkBias=theBestGMN.a0+theBestGMN.a1*dettoc
		 +theBestGMN.a2*dettoc*dettoc+dettr;

	 *detj=theBestGMN.a0+theBestGMN.a1*dettoc
		 +theBestGMN.a2*dettoc*dettoc;
	 
}


