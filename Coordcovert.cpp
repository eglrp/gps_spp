#include "CoordCovert.h"
#include "math.h"

void CartesianToGeodetic (PCRDGEODETIC pcg, PCRDCARTESIAN pcc,
double dSemiMajorAxis, double dFlattening)
{
	double e2;//��һƫ���ʵ�ƽ��
	e2=2*dFlattening-dFlattening*dFlattening;
	//e2=atan(fabspcc->y/pcc->x);
	pcg->longitude=atan(fabs(pcc->y/pcc->x));
	double W,N,N1=0,B,B1;
	B1=atan(pcc->z/sqrt(pcc->x*pcc->x+pcc->y*pcc->y));
	while(1)
	{   
		W=sqrt(1-e2*sin(B1)*sin(B1));
		N1=dSemiMajorAxis/W;
		B=atan((pcc->z+N1*e2*sin(B1))/sqrt(pcc->x*pcc->x+pcc->y*pcc->y));

			if(fabs(B-B1)<delta)
				break;
			else
				B1=B;
	}

	pcg->latitude=B;
	N=dSemiMajorAxis/sqrt(1-e2*sin(pcg->latitude)*sin(pcg->latitude));
	pcg->height=sqrt(pcc->x*pcc->x+pcc->y*pcc->y)/cos(B)-N;
    

}

//�ɴ������ת��Ϊ�ѿ�������
void GeodeticToCartesian (PCRDCARTESIAN pcc, PCRDGEODETIC pcg,
double dSemiMajorAxis, double dFlattening)
{   
	double e2;//��һƫ���ʵ�ƽ��
	double N;//î��Ȧ�뾶
	e2=2*dFlattening-dFlattening*dFlattening;
	N=dSemiMajorAxis/sqrt(1-e2*sin(pcg->latitude)*sin(pcg->latitude));

	pcc->x=(N+pcg->height)*cos(pcg->latitude)*cos(pcg->longitude);
	pcc->y=(N+pcg->height)*cos(pcg->latitude)*sin(pcg->longitude);
	pcc->z=(N*(1-e2)+pcg->height)*sin(pcg->latitude);

}

//�ɵѿ�������ת��Ϊվ�ĵ�ƽ����
void CartesianToTopocentric (PCRDTOPOCENTRIC pct,
PCRDCARTESIAN pcc,
PCRDCARTESIAN pccCenter,
double dSemiMajorAxis,
double dFlattening)
{
	double dx,dy,dz;
	dx=pcc->x-pccCenter->x;
	dy=pcc->y-pccCenter->y;
	dz=pcc->z-pccCenter->z;

	PCRDGEODETIC pd;
	pd=(PCRDGEODETIC)malloc(sizeof(CRDGEODETIC));

    CartesianToGeodetic (pd,pccCenter,dSemiMajorAxis,dFlattening);

	pct->northing=-sin(pd->latitude)*cos(pd->longitude)*dx
		-sin(pd->latitude)*sin(pd->longitude)*dy
		+cos(pd->latitude)*dz;
	pct->easting=-sin(pd->longitude)*dx
		+cos(pd->longitude)*dy;
	pct->upping=cos(pd->latitude)*cos(pd->longitude)*dx
		+cos(pd->latitude)*sin(pd->longitude)*dy
		+sin(pd->latitude)*dz;
	free(pd);

}

//��վ�ĵ�ƽֱ������ת��Ϊվ�ĵ�ƽ������
void TopocentricToTopocentricPolar (PCRDTOPOCENTRICPOLAR pctp,
PCRDTOPOCENTRIC pct)
{   

	pctp->range=sqrt(pct->northing*pct->northing+pct->easting*pct->easting+pct->upping*pct->upping);
	pctp->azimuth=atan(pct->easting/pct->northing);
	pctp->elevation=asin(fabs(pct->upping/pctp->range));


}

//��վ�ĵ�ƽ������ת��Ϊվ�ĵ�ƽֱ������
void TopocentricPolarToTopocentric (PCRDTOPOCENTRIC pct,
PCRDTOPOCENTRICPOLAR pctp)
{
	pct->northing=pctp->range*cos(pctp->elevation)*cos(pctp->azimuth);
	pct->easting=pctp->range*cos(pctp->elevation)*sin(pctp->azimuth);
	pct->upping=pctp->range*sin(pctp->elevation);

}