/************************************************************************/
/* OfflineMapping_basic                                                 */
/*                                                                      */
/* Authors: Yujie Yan, yan.yuj@husky.neu.edu                            */
/*          Jerome Hajjar, jf.Hajjar@northeastern.edu                   */
/* Date: 2/11/2019                                                      */
/************************************************************************/

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "OfflineMapping_basic.h"

/************************************************************************/
/* GPS Data                                                             */
/************************************************************************/
DataGPS::DataGPS():time(50000),timestamp(50000),lat(50000), lon(50000), alt(50000), enu(50000,3)
{
	num_row = 0;
}

void DataGPS::CSVReader() 
{
	std::string line;
	std::ifstream lineStream(filename);
	
	std::getline(lineStream, line, '\n');

	while(std::getline(lineStream, line, '\n'))
	{
		std::istringstream ss(line);
		std::vector<std::string> s_list(50);
		int i = 0;
		while( std::getline(ss, s_list[i], ',') )
		{
			i++;
		}
		//std::string a1 = s_list[0].substr(5,9);
		time[num_row] = std::stod(s_list[0])/1000000000 - 946687000;
		timestamp[num_row] = std::stod(s_list[0])/1000000000 - 946687000;
		lat[num_row] = std::stod(s_list[6]);
		lon[num_row] = std::stod(s_list[7]);
		alt[num_row] = std::stod(s_list[8]);
		num_row++;

//		if (num_row > 1000)
//			break;
	}

	time.conservativeResize(num_row);
	timestamp.conservativeResize(num_row);
	lat.conservativeResize(num_row);
	lon.conservativeResize(num_row);
	alt.conservativeResize(num_row);

	Geodetic2ENU();
}

void DataGPS::Geodetic2ENU()
{
        const double a = 6378137.0;          // WGS-84 Earth semimajor axis (m)
        const double b = 6356752.314245;     // Derived Earth semiminor axis (m)
        const double a_sq = a * a;
        const double b_sq = b * b;
        const double e_sq = 1 - b_sq/a_sq;    // Square of Eccentricity

        //Converts geodetic to Earth-Centered Earth-Fixed (ECEF)
        Eigen::ArrayXd N = sin(lat*M_PI/180);
        N = 1 - N.square() * e_sq;
		Eigen::ArrayXd N2 = N.inverse();
        N = a * N.sqrt().inverse();

        Eigen::ArrayXd sin_lat = (lat*M_PI/180).sin();
        Eigen::ArrayXd cos_lat = (lat*M_PI/180).cos();
        Eigen::ArrayXd sin_lon = (lon*M_PI/180).sin();
        Eigen::ArrayXd cos_lon = (lon*M_PI/180).cos();

        Eigen::ArrayXd x1 = (N + alt)*cos_lat*cos_lon;
        Eigen::ArrayXd y1 = (N + alt)*cos_lat*sin_lon;
        Eigen::ArrayXd z1 = (N*b_sq/a_sq + alt)*sin_lat;
        //Converts ECEF to Local ENU
        x1 = x1-x1(0);
        y1 = y1-y1(0);
        z1 = z1-z1(0);

        Eigen::ArrayXd x = -x1*sin_lon + y1*cos_lon;
        
        Eigen::ArrayXd y = -x1*sin_lat*cos_lon - y1*sin_lat*sin_lon + z1*cos_lat;
        
        Eigen::ArrayXd z =  x1*cos_lat*cos_lon + y1*cos_lat*sin_lon + z1*sin_lat;

        enu.block(0,0,num_row,1) = x;
		enu.block(0,1,num_row,1) = y;
		enu.block(0,2,num_row,1) = z;

		enu.conservativeResize(num_row, 3);
}

void DataGPS::WriteENU()
{
	std::ofstream file("ENU.txt");
	if (file.is_open())
	{
		file << enu << '\n';
	}
}

/************************************************************************/
/* IMU Pose Data                                                        */
/************************************************************************/
DataIMU::DataIMU():time(150000),timestamp(150000), quat(150000,4) 
{
	num_row = 0;
}

void DataIMU::CSVReader() 
{
	std::string line;
	std::ifstream lineStream(filename);
	
	std::getline(lineStream, line, '\n');

	while(std::getline(lineStream, line, '\n'))
	{
		std::istringstream ss(line);
		std::vector<std::string> s_list(300);
		int i = 0;
		while( std::getline(ss, s_list[i], ',') )
		{
			i++;
		}
		//std::string a1 = s_list[0].substr(5,9);
		time[num_row] = std::stod(s_list[0])/1000000000 - 946687000;
		timestamp[num_row] = std::stod(s_list[0])/1000000000 - 946687000;
		quat(num_row,0) = std::stod(s_list[4]);
		quat(num_row,1) = std::stod(s_list[5]);
		quat(num_row,2) = std::stod(s_list[6]);
		quat(num_row,3) = std::stod(s_list[7]);
		num_row++;

//		if (num_row > 10000)
//			break;
	}

	time.conservativeResize(num_row);
	timestamp.conservativeResize(num_row);
	quat.conservativeResize(num_row, 4);
}

void DataIMU::WriteQuat()
{
	std::ofstream file("Quat.txt");
	if (file.is_open())
	{
		file << quat << '\n';
	}
}

/************************************************************************/
/* Process the dji_sdk Data                                             */
/************************************************************************/
double inteploate(Eigen::ArrayXd xData, Eigen::ArrayXd yData, double x, bool extrapolate)
{
	int size = xData.size();

	int i = 0;
	if ( x >= xData[size-2])
	{
		i = size - 2;
	}
	else
	{
		while ( x > xData[i+1]) i++;
	}

	double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];
	if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
   {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
   }
	double dydx = ( yR - yL ) / ( xR - xL );
	
	return yL + dydx * ( x - xL );
}

Eigen::Matrix3d Quat2Rotm(double x, double y, double z, double w)
{
	Eigen::Matrix3d R(3,3);
	R(0,0) = 1 - 2*y*y - 2*z*z;
	R(0,1) = 2*x*y - 2*z*w;
	R(0,2) = 2*x*z + 2*y*w;
	R(1,0) = 2*x*y + 2*z*w;
	R(1,1) = 1 - 2*x*x - 2*z*z;
	R(1,2) = 2*y*z - 2*x*w;
	R(2,0) = 2*x*z - 2*y*w;
	R(2,1) = 2*y*z + 2*x*w;
	R(2,2) = 1 - 2*x*x - 2*y*y;

	return R;
}

Eigen::Matrix3d Axang2Rotm(double x, double y, double z, double theta)
{
	Eigen::Matrix3d R(3,3);
	R(0,0) = cos(theta) + x*x*(1-cos(theta));
	R(0,1) = x*y*(1-cos(theta)) - z*sin(theta);
	R(0,2) = x*z*(1-cos(theta)) + y*sin(theta);
	R(1,0) = y*x*(1-cos(theta)) + z*sin(theta);
	R(1,1) = cos(theta) + y*y*(1-cos(theta));
	R(1,2) = y*z*(1-cos(theta)) - x*sin(theta);
	R(2,0) = z*x*(1-cos(theta)) - y*sin(theta);
	R(2,1) = z*y*(1-cos(theta)) + x*sin(theta);
	R(2,2) = cos(theta) + z*z*(1-cos(theta));
	return R;
}