#include "Pp.h"
#include "BiasCorrection.h"
#include "time.h"
#include "CoordCovert.h"
#include"matrix.h"

#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix<double> Matrix;
#else
typedef matrix Matrix;
#endif

#ifndef _NO_EXCEPTION
#  define TRYBEGIN()	try {
#  define CATCHERROR()	} catch (const STD::exception& e) { \
						cerr << "Error: " << e.what() << endl; }
#else
#  define TRYBEGIN()
#  define CATCHERROR()
#endif




bool PPOne(GMOREC gmoRecord,/*�۲�ֵ��Ϣ��¼*/
		   GMO gmo,/*�۲�ֵ�ļ�*/
		   vector<GMNREC> navRecord,/*���������ļ�*/
		   PPPONERESULT presult,
		   Matrix& Pk,
		   int iflag=1,
		   double X0=1000000.0,
		   double Y0=1000000.0,
		   double Z0=1000000.0
			)
{
	
	double toprelay,ionolay,p0;//�������ӳ�
	int ncount=0;//��������
	////////////////////////////////////////////////////////////////////////////////
	int m_value=-1;//�õ��Ĺ۲��ļ�������α��C1,P1,P2
	int m_p1=-1;
	int m_p2=-1;

    
	//ȷ��ʹ��α�����ڵ�λ��
	for(int type=0;type<gmo.hdr.MeasureTypeNum;type++)
	{
		if(gmo.hdr.ObsType[type].substr(0,2).compare(0,2,"C1")==0)
		{
			m_value=type;
	
		}
		 if(gmo.hdr.ObsType[type].substr(0,2).compare(0,2,"P1")==0)
		{
			m_p1=type;
		
		}
		 if(gmo.hdr.ObsType[type].substr(0,2).compare(0,2,"P2")==0)
		{
			m_p2=type;
		
		}
		
	}

	//ȷ����ʹ�õ�������Ŀ
    int bz=0,vp=0,vh=0,t=0;
	int sat_value[MAXNUM]={0};
	int sat_pos[MAXNUM]={0};
	int PN;


	///////////////////////////////////////////////////////
	//
	double DetTi;
	double dx,dy,dz,dt;
	if(iflag==1&&gmo.hdr.AppX!=0)
	{
		X0=gmo.hdr.AppX;
		Y0=gmo.hdr.AppY;
		Z0=gmo.hdr.AppZ;
	}
		

        
	DetTi=0;//��ʼ��ֵ

        
	double m_dop=0;
	double X,Y,Z;
	double dett2,dett1;
	double tbias,detj,R=20000000;

	COMMONTIME   tj;//�źŷ���ʱ��
	tj.year=gmoRecord.epochtime.year;
	tj.month=gmoRecord.epochtime.month;
	tj.day=gmoRecord.epochtime.day;
	tj.hour=gmoRecord.epochtime.hour;
	tj.minute=gmoRecord.epochtime.minute;
	tj.second=gmoRecord.epochtime.second;

		 
	JULIANDAY tjJ;//����ʱ��
	JULIANDAY epo;//��Ԫʱ��

		///////////////////////////////
		 //��������������ӳ�
	CRDCARTESIAN crdSite;//��վ����
	CRDCARTESIAN crdSat;//��������
   
	CRDCARTESIAN pcrdOrb1;

	crdSite.x=X0;
	crdSite.y=Y0;
	crdSite.z=Z0;

	

	int pre=0;//�ж�α��۲�ֵ����P1,P2����C1,Ϊ0ʱ��ʾ��P1,P2.

	if(m_p1!=-1&&m_p2!=-1)//�����P1,P2�Ͳ���˫Ƶ����ģ���������������ӳ�
	{   
		for(int j=0;j<gmoRecord.satsum;j++)
		{
			PN=atoi(gmoRecord.PRN_list[j].substr(1,2).c_str());

			p0=gmoRecord.obsValue[j][m_value];
			
			dett1=p0/c;//����ʱ��

			
			GetSVClkBias(navRecord,PN,&tj,&tbias,&detj);

			CommonTimeToJulianDay(&gmoRecord.epochtime, &epo);

			SetTimeDelta (&tjJ, &epo, (-tbias-dett1));
			JulianDayToCommonTime (&tjJ, &tj);
         


			//���������ڵѿ�������ϵ�е�λ��
			GetOrbNClk(navRecord,PN,&tj,&pcrdOrb1);
			//��������C/A���źŷ���ʱ�̵ĸ���
			//GetSVClkBias(navRecord,PN,&tj,&tbias,&detj);

			//������ת�ĸ���,4-120
			X=cos(dett1*we)*pcrdOrb1.x+sin(dett1*we)*pcrdOrb1.y;
			Y=-sin(dett1*we)*pcrdOrb1.x+cos(dett1*we)*pcrdOrb1.y;
			Z=pcrdOrb1.z;

			crdSat.x=X;
			crdSat.y=Y;
			crdSat.z=Z;

			CRDTOPOCENTRIC  pct;

			CartesianToTopocentric (&pct,&crdSat,&crdSite,a,flattening);
		
			CRDTOPOCENTRICPOLAR site;

			TopocentricToTopocentricPolar (&site,&pct);
			double E=site.elevation;

			if(E<(15.0*PI/180.0))
			{
				t++;
				continue;
			}


			for(int k=0;k<navRecord.size();k++)
			{
				
				if(PN==navRecord[k].PRN&&gmoRecord.obsValue[j][m_p1]!=0&&gmoRecord.obsValue[j][m_p2]!=0)
				{
					sat_value[vp]=PN;
					vp++;
					sat_pos[vh]=j;
					vh++;
					bz=1;

				
					break;
				}

			}
			if(bz==0)
				t++;
			else
				bz=0;
		}
	}
	else//���û��P1��P2����ֱ����C1����α�࣬���Ե����Ӱ��
	{   
		pre=1;
		for(int j=0;j<gmoRecord.satsum;j++)
		{
			PN=atoi(gmoRecord.PRN_list[j].substr(1,2).c_str());//��i������PRN��
			p0=gmoRecord.obsValue[j][m_value];//C1α��ֵ
			
			dett1=p0/c;//����ʱ��

			GetSVClkBias(navRecord,PN,&tj,&tbias,&detj);
			CommonTimeToJulianDay(&gmoRecord.epochtime, &epo);


			SetTimeDelta (&tjJ, &epo, (-tbias-dett1));
			JulianDayToCommonTime (&tjJ, &tj);

			//���������ڵѿ�������ϵ�е�λ��
			GetOrbNClk(navRecord,PN,&tj,&pcrdOrb1);
			//��������C/A���źŷ���ʱ�̵ĸ���
			//GetSVClkBias(navRecord,PN,&tj,&tbias,&detj);

			//������ת�ĸ���,4-120
			X=cos(dett1*we)*pcrdOrb1.x+sin(dett1*we)*pcrdOrb1.y;
			Y=-sin(dett1*we)*pcrdOrb1.x+cos(dett1*we)*pcrdOrb1.y;
			Z=pcrdOrb1.z;

			crdSat.x=X;
			crdSat.y=Y;
			crdSat.z=Z;

			CRDTOPOCENTRIC  pct;

			CartesianToTopocentric (&pct,&crdSat,&crdSite,a,flattening);
		
			CRDTOPOCENTRICPOLAR site;

			TopocentricToTopocentricPolar (&site,&pct);
			double E=site.elevation;

			if(E<(15.0*PI/180.0))
			
			{
				t++;
				continue;
			}
			for(int k=0;k<navRecord.size();k++)
			{
					//cout<<"  test 1 "<<endl;
				if(PN==navRecord[k].PRN&&gmoRecord.obsValue[j][m_value]!=0)  //��k��N�ļ���¼
				{
					sat_value[vp]=PN;
					vp++;
					sat_pos[vh]=j;
					vh++;
					bz=1;
					
					break;
				}

			}
			if(bz==0)
				t++;
			else
				bz=0;
				//cout<<"  test 3 "<<endl;
		}

		//cout<<"  test 4 "<<endl; 
	}
		
    int sat_valuesum=gmoRecord.satsum-t;//��������������
	int len=sat_valuesum;

	//cout<<"  test 5 "<<endl;

	if(sat_valuesum>=4)
	{

		Matrix B(len,4),W(len,1),Result(4,1),Q(4,4);
		Matrix Me(4,4),Qk(4,4),Rk(len,len),Kk(4,len);


		for(int i=0;i<=3;i++)
		{
			for(int j=0;j<=3;j++)
			{
				if(i==j)
				{
					Me(i,j)=1;
					Qk(i,j)=10.0;
				}
					else
				{
					Me(i,j)=0;
					Qk(i,j)=0;
				}	

			}
		}
		Qk(3,3)=1;
		
		for( i=0;i<=len-1;i++)
		{
			
		}

		int PRN;//��Ч���ǵ����Ǻ�		
		//double p0;//α��
	    int pos;//��Ч�����ڹ۲������е�λ��




		////////////////////////////////////////////////////
		do
		{
			 ncount++;
			 crdSite.x=X0;
			 crdSite.y=Y0;
			 crdSite.z=Z0;
			 for(int i=0;i<len;i++)
			 {
				
				 PRN=sat_value[i];
				 //��ʼ�������Ӳ�
				 GetSVClkBias(navRecord,PRN,&gmoRecord.epochtime,&tbias,&detj);

				 pos=sat_pos[i];
				
				 if(pre==0)//��P1,P2;4-61
					 p0=gmoRecord.obsValue[pos][m_p1]*2.54573
						 -gmoRecord.obsValue[pos][m_p2]*1.54573;
				 else//��C1
					 p0=gmoRecord.obsValue[pos][m_value];
			
				 dett2=p0/c;//����ʱ��
           
				 do
				 {
					
					 dett1=dett2;
					 //ʱ�����
					 CommonTimeToJulianDay(&gmoRecord.epochtime, &epo);
					 SetTimeDelta (&tjJ, &epo, (-tbias-dett1));
					 JulianDayToCommonTime (&tjJ, &tj);
                

					 //���������ڵѿ�������ϵ�е�λ��
					 GetOrbNClk(navRecord,PRN,&tj,&pcrdOrb1);
					 //��������C/A���źŷ���ʱ�̵ĸ���
					 GetSVClkBias(navRecord,PRN,&tj,&tbias,&detj);

					//������ת�ĸ���,4-120
					 X=cos(dett1*we)*pcrdOrb1.x+sin(dett1*we)*pcrdOrb1.y;
					 Y=-sin(dett1*we)*pcrdOrb1.x+cos(dett1*we)*pcrdOrb1.y;
					 Z=pcrdOrb1.z;

					 crdSat.x=X;
					 crdSat.y=Y;
					 crdSat.z=Z;


					 R=sqrt((X-X0)*(X-X0)+(Y-Y0)*(Y-Y0)+(Z-Z0)*(Z-Z0));//�������
			  
					 dett2=R/c;//����ʱ��

				}while(fabs(dett2-dett1)>1e-12);

				 B(i,0)=(X0-X)/R;
				 B(i,1)=(Y0-Y)/R;
				 B(i,2)=(Z0-Z)/R;
				 B(i,3)=1;
         
				 toprelay=GetTropDelay(&crdSite,&crdSat);//�������ӳ�



				 CRDTOPOCENTRIC  pct;
				 CartesianToTopocentric (&pct,&crdSat,&crdSite,a,flattening);		
				 CRDTOPOCENTRICPOLAR site;
				 TopocentricToTopocentricPolar (&site,&pct);
				 double E=site.elevation;
				 
				 for(int j=0;j<=len-1;j++)
				 {
					if(i==j)
					{
						Rk(i,j)=E*2.0/PI;
					}
						else
					{
						Rk(i,j)=0;
					}	

				 }


				 W(i,0)=p0-R+c*detj-c*DetTi +toprelay;
				 
				 //cout<<R<<endl;
			 
			 }
     
			if(iflag==1)
			{
				/*for(i=0;i<=len-1;i++)
				{
					for(int j=0;j<=3;j++)
					{
						cout<<B(i,j)<<"  ";
					}
					cout<<endl;
				}*/
				Result=!(~B*B)*~B*W;
				Pk=~B*B;
				Q=!Pk;
			}
			else if(iflag==2)
			{
				/*cout<<"len="<<len<<endl;
				for(i=0;i<=len-1;i++)
				{
					for(int j=0;j<=3;j++)
					{
						cout<<B(i,j)<<"  ";
					}
					cout<<endl;
				}*/
				Pk=Pk+Qk;

				/*for(i=0;i<=3;i++)
				{
					for(int j=0;j<=3;j++)
					{
						cout<<Pk(i,j)<<"  ";
					}
					cout<<endl;
				}*/

				//B=B*Pk;

				/*for(i=0;i<=len-1;i++)
				{
					for(int j=0;j<=3;j++)
					{
						cout<<B(i,j)<<"  ";
					}
					cout<<endl;
				}*/

				Kk=Pk*~B*!(B*Pk*~B+Rk);
				
				/*for(i=0;i<=3;i++)
				{
					for(int j=0;j<=len-1;j++)
					{
						cout<<Kk(i,j)<<"  ";
					}
					cout<<endl;
				}*/

				Result=Kk*W;
				Pk=(Me-Kk*B)*Pk;
				Q=!Pk;
			}

			dx=Result(0,0);
			dy=Result(1,0);
			dz=Result(2,0);
			dt=Result(3,0);
			X0+=dx;
			Y0+=dy;
			Z0+=dz;
			DetTi+=dt/c;

			if(ncount>10)//����10 �����ϵĵ�����Ϊʧ��
			{	
				return false;
				break;
			}
		}while(fabs(dx)>0.1||fabs(dy)>0.1||fabs(dz)>0.1);

		for(int q=0;q<3;q++)
			m_dop+=Q(q,q);

				 
      presult->sat_num=len;
      presult->clk_bias=DetTi*1e9;

	  presult->crd.x=X0;
	  presult->crd.y=Y0;
	  presult->crd.z=Z0;

	  presult->epoch.year=gmoRecord.epochtime.year;
	  presult->epoch.month=gmoRecord.epochtime.month;
	  presult->epoch.day=gmoRecord.epochtime.day;
	  presult->epoch.hour=gmoRecord.epochtime.hour;
	  presult->epoch.minute=gmoRecord.epochtime.minute;
	  presult->epoch.second=gmoRecord.epochtime.second;

	  presult->PDOP=sqrt(m_dop);

	  presult->crdsigma.x=Pk(0,0);
	  presult->crdsigma.y=Pk(1,1);
	  presult->crdsigma.z=Pk(2,2);
      
	  
	   return true;
	
	}

	
	else
	{
		return false;
	}
  

}


vector<PPONERESULT> PP(GMO gmo,/*�۲�ֵ�ļ�*/
		               vector<GMNREC> navRecord/*���������ļ�*/
					   )
{
	vector<PPONERESULT> m_result;
	double crdX0,crdY0,crdZ0;
	PPONERESULT m_oneresult;
    Matrix Pk;
	int iflag;


	for(int i=0;i<=3;i++)
	{
		for(int j=0;j<=3;j++)
		{
			if(i==j)
				Pk(i,j)=1;
			else
				Pk(i,j)=0;

		}
	}

	int len=gmo.obs.size();
	if(len!=0)
	{
		cout<<"��ѡ��ƽ���:1��С����,2�������˲�;"<<endl;
		cin>>iflag;
		cout<<"���ݴ�����..."<<endl;
		for(int j=0;j<len;j++)
		{
			if(PPOne(gmo.obs[j], gmo, navRecord,&m_oneresult,Pk))
			{
		
				m_result.push_back(m_oneresult);

				crdX0=m_result[0].crd.x;
				crdY0=m_result[0].crd.y;
				crdZ0=m_result[0].crd.z;
				break;
			}
		}
		

		//cout<<"pre_m_oneresult.crd.x="<<m_oneresult.crd.x<<endl;



		int k=1;
		for(int i=j+1;i<len-1;i++)
		{
			if(PPOne(gmo.obs[i], gmo, navRecord,&m_oneresult,Pk,iflag,crdX0,crdY0,crdZ0))
			{
				//pre_m_oneresult=m_oneresult; 
				
				m_result.push_back(m_oneresult);
				crdX0=m_result[k].crd.x;
				crdY0=m_result[k].crd.y;
				crdZ0=m_result[k].crd.z;
				cout<<"���ڴ�����Ԫ��"<<k<<endl;
				k++;
			}
		}
		cout<<"���ݴ�����ϣ�"<<endl;
	}
	else
	{
		cout<<"û�з���Ҫ������ݣ�"<<endl;
	}

	return m_result;
}

