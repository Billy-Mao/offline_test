/************************************************************************/
/*                                                                      */
/* Authors: Yujie Yan, yan.yuj@husky.neu.edu                            */
/*          Jerome Hajjar, jf.Hajjar@northeastern.edu                   */
/* Date: 2/11/2019                                                      */
/************************************************************************/

#include <Eigen/StdVector>
#include <Eigen/Geometry>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

/************************************************************************/
/* GPS Data                                                             */
/************************************************************************/
class DataGPS{
public:
	Eigen::ArrayXd time;
	Eigen::ArrayXd timestamp;
	Eigen::ArrayXd lat;
	Eigen::ArrayXd lon;
	Eigen::ArrayXd alt;
	Eigen::ArrayX3d enu;
	std::string filename;
	int num_row;

	DataGPS();
	void CSVReader();
	void Geodetic2ENU();
	void WriteENU();
};

/************************************************************************/
/* IMU Pose Data                                                        */
/************************************************************************/
class DataIMU{
	public:
	Eigen::ArrayXd time;
	Eigen::ArrayXd timestamp;
	Eigen::ArrayXXd quat;
	int num_row;
	std::string filename;

	DataIMU();
	void CSVReader();
	void WriteQuat();
};

/************************************************************************/
/* Process the dji_sdk Data                                             */
/*                                                                      */
/************************************************************************/
double inteploate(Eigen::ArrayXd xData, Eigen::ArrayXd yData, double x, bool extrapolate);

Eigen::Matrix3d Quat2Rotm(double x, double y, double z, double w);

Eigen::Matrix3d Axang2Rotm(double x, double y, double z, double theta);