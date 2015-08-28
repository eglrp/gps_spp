#include "RinexObsRead.h"



GMO ReadRinexObsFile(string fp)
{
	 GMO obsData;
	 obsData.hdr.MeasureInterval=0;//��Ϊ�������п��п��ޣ���ֹ�����ȸ�ֵ
     int obstypesum=0;
	 cout<<"ReadingRinexOavFile..."<<endl;
	 int satsum=0;
		ifstream inf;//C++��ʽ�������ļ�
		inf.open(fp.c_str());//ע��·��

		if(!inf)
		 {
		 cout<<"�ļ��򿪴���";
		 exit(0);
		 }

		string obs_header;
		getline(inf,obs_header);//���ļ�ͷ

		while(obs_header.substr(60,72).compare(0,13,"END OF HEADER")!=0)
			
		{
			
			if(obs_header.substr(60,79).compare(0,20,"RINEX VERSION / TYPE")==0)
			{
				obsData.hdr.FormatVersion   =  obs_header.substr(0,9);
				obsData.hdr.FileTypeObsStr  =  obs_header.substr(20,1);
				obsData.hdr.PositionSystem  =  obs_header.substr(40,1);
				//cout<<obsData.hdr.PositionSystem<<endl;
				
			}
			if(obs_header.substr(60,78).compare(0,19,"APPROX POSITION XYZ")==0)
			{
				obsData.hdr.AppX =  atof(obs_header.substr(0,14).c_str());
				obsData.hdr.AppY =  atof(obs_header.substr(14,14).c_str());
				obsData.hdr.AppZ =  atof(obs_header.substr(28,14).c_str());
				//cout<<obsData.hdr.AppX<<"   "<<obsData.hdr.AppY<<"   "<<obsData.hdr.AppZ<<endl;
			
			}
			if(obs_header.substr(60,79).compare(0,20,"ANTENNA: DELTA H/E/N")==0)
			{
				obsData.hdr.AntHeight  =  atof(obs_header.substr(0,14).c_str());
				obsData.hdr.AntEast    =  atof(obs_header.substr(14,14).c_str());
				obsData.hdr.AntWest    =  atof(obs_header.substr(28,14).c_str());
				//cout<<obsData.hdr.AntHeight<<"   "<<obsData.hdr.AntEast<<"   "<<obsData.hdr.AntWest<<endl;
			
			}
			if(obs_header.substr(60,79).compare(0,20,"WAVELENGTH FACT L1/2")==0)
			{
				obsData.hdr.WaveFract  =   atoi(obs_header.substr(0,6).c_str());
				obsData.hdr.SDFreq     =   atoi(obs_header.substr(6,6).c_str());
				//cout<<obsData.hdr.WaveFract<<"   "<<obsData.hdr.SDFreq<<endl;
			
			}
			if(obs_header.substr(60,78).compare(0,19,"# / TYPES OF OBSERV")==0)
			{
				obsData.hdr.MeasureTypeNum  =  atoi(obs_header.substr(0,6).c_str());

                obstypesum=obsData.hdr.MeasureTypeNum;
				for(int i=0;i<obstypesum;i++)
				{
					obsData.hdr.ObsType[i]  =  obs_header.substr(10+i*6,2);
				
				}
				//cout<<obsData.hdr.MeasureTypeNum<<endl;
				//cout<<obsData.hdr.ObsType[0]<<"   "<<obsData.hdr.ObsType[1]<<endl;
			}
			if(obs_header.substr(60,67).compare(0,8,"INTERVAL")==0)
			{
				obsData.hdr.MeasureInterval  =  atof(obs_header.substr(0,10).c_str());
				//cout<<obsData.hdr.MeasureInterval<<endl;
			}
			if(obs_header.substr(60,76).compare(0,17,"TIME OF FIRST OBS")==0)
			{
				obsData.hdr.m_startTime.year   =  atoi(obs_header.substr(6,6).c_str());
				obsData.hdr.m_startTime.month  =  atoi(obs_header.substr(6,6).c_str());
				obsData.hdr.m_startTime.day    =  atoi(obs_header.substr(12,6).c_str());
				obsData.hdr.m_startTime.hour   =  atoi(obs_header.substr(18,6).c_str());
				obsData.hdr.m_startTime.minute =  atoi(obs_header.substr(24,6).c_str());
				obsData.hdr.m_startTime.second =  atof(obs_header.substr(30,13).c_str());

				obsData.hdr.timeSystem         =  obs_header.substr(48,3);
				//cout<<obsData.hdr.m_startTime.year<<"  "<<obsData.hdr.m_startTime.month<<"  ";
				//cout<<obsData.hdr.m_startTime.day<<"   "<<obsData.hdr.m_startTime.hour<<"   ";
				//cout<<obsData.hdr.m_startTime.minute<<"   "<<obsData.hdr.m_startTime.second<<endl;
			}

			getline(inf,obs_header);
			
		}
		
		GMOREC obs;
		int sum;//���ǵ���Ŀ
		string obsvalue_header,obs_str;

		while(inf.peek() != EOF)
		{
			/*obs.epochtime.year   =0;
			obs.epochtime.month  =0;
			obs.epochtime.day    =0;
			obs.epochtime.hour   =0;
			obs.epochtime.minute =0;
			obs.epochtime.second =0;*/
			memset(&obs.epochtime,0,sizeof(obs.epochtime));
			memset(obs.PRN_list,0,sizeof(string)*12);
			memset(obs.obsValue,0,sizeof(double)*108);
			obs.flag          = 0;
			obs.sat_time_bias = 0;
			obs.satsum        = 0;//��ʼ�������⸳ֵ����
			

			//��һ��
           getline(inf,obsvalue_header);
		   obs.epochtime.year=atoi(obsvalue_header.substr(1,2).c_str());
		   
		   if(obs.epochtime.year>90)
			   obs.epochtime.year+=1900;
		   else
			   obs.epochtime.year+=2000;

		   obs.epochtime.month   =  atoi(obsvalue_header.substr(4,2).c_str());
		   obs.epochtime.day     =  atoi(obsvalue_header.substr(7,2).c_str());
		   obs.epochtime.hour    =  atoi(obsvalue_header.substr( 10,2).c_str());
		   obs.epochtime.minute  =  atoi(obsvalue_header.substr(13,2).c_str());
		   obs.epochtime.second  =  atof(obsvalue_header.substr(15,11).c_str());

		   obs.flag              =  atoi(obsvalue_header.substr(28, 1).c_str());
		   obs.satsum            =  atoi(obsvalue_header.substr(29,3).c_str());

		   //cout<<obs.epochtime.year<<"  "<<obs.epochtime.month<<"   "<<obs.epochtime.day<<"   ";
		  // cout<<obs.epochtime.hour<<"   "<<obs.epochtime.minute<<"  "<<obs.epochtime.second<<endl;

		  // cout<<obs.flag<<endl;
		  // cout<<obs.satsum<<endl;
           sum  =  obs.satsum;

		   for(int sat=0;sat<sum;sat++)//�������б�
		   {
			   obs.PRN_list[sat] =  obsvalue_header.substr(32+sat*3,3);
			   
		   }
		  // cout<<obs.PRN_list[0]<<"   "<<obs.PRN_list[sum-1]<<endl;

		   if(obsvalue_header.size()>68)//�ж��Ƿ���ʱ��ƫ�������
		        obs.sat_time_bias = atof(obsvalue_header.substr(68,12).c_str());

		  // cout<<obs.sat_time_bias<<endl;

          if(obstypesum<=5)//С��5���۲�ֵ���͵����
		  {
			  for(int satID=0;satID<sum;satID++)
			  {
				  getline(inf,obs_str);
					  for(int satValue=0;satValue<obstypesum;satValue++)
					  {
						  if((int)obs_str.size()>16*satValue)
						       obs.obsValue[satID][satValue] = atof(obs_str.substr(16*satValue,14).c_str());
						 // cout<<obs.obsValue[satID][satValue]<<endl;

					  }
			  }
		  
		  }

		  else//��������۲�ֵ���͵����
		  {
			  for(int satID2=0;satID2<sum;satID2++)
			  {
				  getline(inf,obs_str);
				  for(int satValue=0;satValue<5;satValue++)
					  {
						  if((int)obs_str.size()>16*satValue)
						    obs.obsValue[satID2][satValue] = atof(obs_str.substr(16*satValue,14).c_str());

					  }
				  getline(inf,obs_str);
				  for(int satValue2=5;satValue2<obstypesum;satValue2++)
					  {
						  if((int)obs_str.size()>16*(satValue2-5))
						     obs.obsValue[satID2][satValue2] = atof(obs_str.substr(16*(satValue2-5),14).c_str());

					  }
 
			  }
		  
		  }


	 obsData.obs.push_back(obs);

		}
		
		return obsData;

}