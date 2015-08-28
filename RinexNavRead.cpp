#include "RinexNavRead.h"


vector<GMNREC> ReadRinexNavFile(string fp)
{
	    GMNREC navdata;
		vector<GMNREC> navRecord;

		COMMONTIME timetran;
        double toed; 
		cout<<"ReadingRinexNavFile..."<<endl;
		ifstream inf;//C++方式打开星历文件
		inf.open(fp.c_str());//注意路径

		if(!inf)
		 {
		 cout<<"文件打开错误";
		 exit(0);
		 }

		string str_header;
		getline(inf,str_header);//读文件头

		while(str_header.substr(60,72).compare(0,13,"END OF HEADER")!=0)//忽略文件头
			
		{
			getline(inf,str_header);
			
		}
		string strData0,strData1,strData2,strData3,strData4,strData5,strData6,strData7;

		

		while(inf.peek() != EOF)
		{   
					getline(inf,strData0);//读第0行
					navdata.PRN    =    atoi(strData0.substr(0,2).c_str());
					timetran.year   =    atoi(strData0.substr(3,2).c_str());
					     if(timetran.year>80)
								timetran.year+=1900;
							else
								timetran.year+=2000;
					timetran.month  =    atoi(strData0.substr(6,2).c_str());
					timetran.day    =    atoi(strData0.substr(9,2).c_str());
					timetran.hour   =    atoi(strData0.substr(12,2).c_str());
					timetran.minute =    atoi(strData0.substr(15,2).c_str());
					timetran.second =    atof(strData0.substr(17,5).c_str());
                    
					CommonTimeToGPSTime (&timetran, &navdata.TOC);

					navdata.a0     =    atof(strData0.substr(22,19).c_str());
					navdata.a1     =    atof(strData0.substr(41,19).c_str());
					navdata.a2     =    atof(strData0.substr(60,19).c_str());

					getline(inf,strData1);//读第一行
					navdata.IODE   =    atof(strData1.substr(3,19).c_str());
					navdata.Crs    =    atof(strData1.substr(22,19).c_str());
					navdata.deltn  =    atof(strData1.substr(41,19).c_str());
					navdata.M0     =    atof(strData1.substr(60,19).c_str());

					getline(inf,strData2);//读第二行
					navdata.Cuc    =    atof(strData2.substr(3,19).c_str());
					navdata.e      =    atof(strData2.substr(22,19).c_str());
					navdata.Cus    =    atof(strData2.substr(41,19).c_str());
					navdata.SqrtA  =    atof(strData2.substr(60,19).c_str());

					getline(inf,strData3);//读第三行
					//navdata.toe    =    atof(strData3.substr(3,19).c_str());
					toed           =    atof(strData3.substr(3,19).c_str());
					navdata.TOE.tow.sn=static_cast<long>(toed);
					navdata.TOE.tow.tos=toed-static_cast<long>(toed);

					navdata.Cic    =    atof(strData3.substr(22,19).c_str());
					navdata.omiga0 =    atof(strData3.substr(41,19).c_str());
					navdata.Cis    =    atof(strData3.substr(60,19).c_str());

					getline(inf,strData4);//读第四行
					navdata.i0       =  atof(strData4.substr(3,19).c_str());
					navdata.Crc      =  atof(strData4.substr(22,19).c_str());
					navdata.omiga    =  atof(strData4.substr(41,19).c_str());
					navdata.omigaDot =  atof(strData4.substr(60,19).c_str());

					getline(inf,strData5);//读第五行
					navdata.iDot              =  atof(strData5.substr(3,19).c_str());
					navdata.CodesOnL2Chanel   =  atof(strData5.substr(22,19).c_str());
					//navdata.weekno            =  atof(strData5.substr(41,19).c_str());
					navdata.TOE.wn            =  static_cast<int>(atof(strData5.substr(41,19).c_str()));
					navdata.L2PdataFlag       =  atof(strData5.substr(60,19).c_str());

					getline(inf,strData6);//读第六行
					navdata.SVAccuracy  =  atof(strData6.substr(3,19).c_str());
					navdata.SVHealth    =  atof(strData6.substr(22,19).c_str());
					navdata.tgd         =  atof(strData6.substr(41,19).c_str());
					navdata.todc        =  atof(strData6.substr(60,19).c_str());

					getline(inf,strData7);//读第七行
					navdata.TransTimeofMsg  =  atof(strData7.substr(3,19).c_str());
					navdata.spare1          =  atof(strData7.substr(22,19).c_str());
		            navdata.spare2          =  0;
		            navdata.spare3          =  0;
//
					navRecord.push_back(navdata);

				
       
	   
		  }

         inf.close();//关闭文件


  return navRecord;

}