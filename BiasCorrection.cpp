#include "BiasCorrection.h"
#include "math.h"
#include "set.h"

double GetTropDelay(PCRDCARTESIAN pcrdSite,PCRDCARTESIAN pcrdSat)
{
	CRDGEODETIC crdSite;
	CartesianToGeodetic (&crdSite, pcrdSite,a,flattening);

	double H=crdSite.height;

	CRDTOPOCENTRIC  pct;

	CartesianToTopocentric (&pct,pcrdSat,pcrdSite,a,flattening);
    
	CRDTOPOCENTRICPOLAR site;
    TopocentricToTopocentricPolar (&site,&pct);

	double E=site.elevation;
	double delta=0;
	//Hopfield,4-90
	if(fabs(H)<40000)
	{
	
	double T=T0-0.0065*(H-H0)+273.16;
	double P=P0*pow((1-0.0000226*(H-H0)),5.225);
	double RH=RH0*exp(-0.0006396*(H-H0));

	double e=RH*exp(-37.2465+0.213166*T-0.000256908*T*T);
	double hw=11000;
	double hd=40136+148.72*(T-273.16);

	double Kw=(155.2e-7*4810*e*(hw-H))/(T*T);
	double Kd=(155.2e-7*P*(hd-H))/T;
	 delta=Kd/(sin(sqrt(E*E+6.25)))+Kw/(sin(sqrt(E*E+2.25)));

	}
	else
     delta=0;
	return delta;
}


//double GetIonolay()

/*double GetIonolay(vector<GMNREC> navRecord,int nPRN,
		PCOMMONTIME pctEpoch,
		PCRDCARTESIAN pcrdSat,
		PCRDCARTESIAN pcrdSite,
		double* pionolay)
{
	CRDGEODETIC cgsite;
	CRDTOPOCENTRIC  pct;
	CRDTOPOCENTRICPOLAR site;
	
	UtilParam pParam;
	GetUtilParameter(navRecord,nPRN,pctEpoch,&pParam);
	    
	GMNREC  theBestGMN;
	theBestGMN=GetBestGMNREC(navRecord,nPRN,pctEpoch);	

	CartesianToTopocentric (&pct,&crdSat,&crdSite,a,flattening);

	TopocentricToTopocentricPolar (&site,&pct);
	double E=site.elevation;

	double alpha=site.azimuth;
	
	double EA=(445.0*PI/180.0)/(E+20.0*PI/180.0)-4.0*PI/180.0;//4-53

	CartesianToGeodetic (*cgsite,pcrdSite,a,flattening);

	double PHIIP=cgsit.latitude+EA*cos(alpha);
	double LAMDA=cgsit.longtitude+EA*sin(alpha)/cos(cgsit.latitude);
	double phim=PHIP+(10.07*PI/180.0)*cos(LAMDA-288.04*PI/180.0);

	double A=h
    
	 
}*/
