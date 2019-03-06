/************************************************************************/
/* OfflineMapping_basic                                                 */
/* Local map constructed using pre-calibrated transformation            */
/* Local maps are registered using ICP naively                          */
/* To provide a prototype for integrating data streaming                */
/*                                                                      */
/* Authors: Yujie Yan, yan.yuj@husky.neu.edu                            */
/*          Jerome Hajjar, jf.Hajjar@northeastern.edu                   */
/* Date: 2/11/2019                                                      */
/************************************************************************/

#include "OfflineMapping_basic.h"

int main (int argc, char** argv)
{   //Initialize Inputs
//	std::string DataName = "180831_UAS_Lab";
	int RPM_velo = 1200, RPM_gimbal = 6;
	int SPR_gimbal = 60/RPM_gimbal;
	std::string TestName = "test2_" + std::to_string(static_cast<long long>(SPR_gimbal)) + "spr_" \
		                    + std::to_string(static_cast<long long>(RPM_velo)) + "RPM";
	//std::string Path2pcd = "E:/Data/DHS_Project/" + DataName + "/" + TestName + "/pcds/";
	//std::string Path2sdk = "E:/Data/DHS_Project/" + DataName + "/" + TestName + "/dji_sdk/";
	//std::string Path2pcd = "C:/Users/yanyujie0824/Documents/PHD_workspace/cpp_mapping/cpp_test2_10spr_1200RPM/pcds/";
	//std::string Path2sdk = "C:/Users/yanyujie0824/Documents/PHD_workspace/cpp_mapping/cpp_test2_10spr_1200RPM/dji_sdk/";
        // std::string Path2pcd = "/home/billymao/Desktop/Zhong_ws/OfflineMapping_basic/Data/pcds/";
	//std::string Path2sdk = "/home/billymao/Desktop/Zhong_ws/OfflineMapping_basic/Data/dji_sdk/";
         std::string Path2pcd = "/home/ubuntu/offline_mapping/Data/pcds/";
	 std::string Path2sdk = "/home/ubuntu/offline_mapping/Data/dji_sdk/";




	//Read user-inputs

	//Compute Key Parameters
	double A = -(360*RPM_gimbal/(double)RPM_velo)*M_PI/180;
	Eigen::Matrix<double, 3, 1> z_axis(0.9748076, -0.2188575, 0.0102304);
	Eigen::Matrix4d T_sensor2FLU;
	T_sensor2FLU << 0.1879196,  0.9739312,  -0.1270599,  0,
	0.0243581, -0.1339460,  -0.9906892,  0,
	-0.9818823,  0.1830751,  -0.0488942, -0.3065,
	0,          0,           0,          1;


	//Read GPS information
	DataGPS gps1;
	gps1.filename = Path2sdk + "gps.csv";
	gps1.CSVReader();
	gps1.WriteENU();

	DataIMU imu1;
	imu1.filename = Path2sdk + "imu.csv";
	imu1.CSVReader();
	//imu1.WriteQuat();

	//start mapping
	pcl::PointCloud<pcl::PointXYZ>::Ptr map (new pcl::PointCloud<pcl::PointXYZ> ());
	pcl::PointCloud<pcl::PointXYZ>::Ptr map_down (new pcl::PointCloud<pcl::PointXYZ> ());
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ> ());
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_trans (new pcl::PointCloud<pcl::PointXYZ> ());
	pcl::PointCloud<pcl::PointXYZ>::Ptr map_local (new pcl::PointCloud<pcl::PointXYZ> ());
	pcl::PointCloud<pcl::PointXYZ>::Ptr map_local_down (new pcl::PointCloud<pcl::PointXYZ> ());
	pcl::PointCloud<pcl::PointXYZ>::Ptr map_local_down_trans (new pcl::PointCloud<pcl::PointXYZ> ());
	pcl::PointCloud<pcl::PointXYZ>::Ptr map_local_down_icp (new pcl::PointCloud<pcl::PointXYZ> ());
	std::string name_pcd, name_out, path2cloud;
	std::ifstream lineStream(Path2pcd + "indices.txt");
	std::ofstream file1;
	int num_scan = 0, num_scan_start = 101, nn = 0, i;
	double tt_pcd, xi_gps, yi_gps, zi_gps, xi_imu, yi_imu, zi_imu, wi_imu;
	Eigen::Matrix4d RR_init, rotm1, rotm2, RRi;
	pcl::VoxelGrid<pcl::PointXYZ> sor;

	std::string path2map = "/home/ubuntu/offline_mapping/Data/map_v2p.pcd";
    pcl::io::loadPCDFile (path2map, *map);
	sor.setInputCloud(map);
	sor.setLeafSize(0.15f, 0.15f, 0.15f);
    sor.filter (*map_down);

	//name_out = "map_down.pcd";
    //name_out = "C:/Users/yanyujie0824/Documents/PHD_workspace/cpp_mapping/cpp_test2_10spr_1200RPM/cpp_map_v2p_point2point/" + name_out;
    //pcl::io::savePCDFile (name_out, *map_down, true);

	pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;
	icp.setTransformationEpsilon (1e-8);
    icp.setMaxCorrespondenceDistance (5);
    icp.setMaximumIterations(10000);


	RR_init.setIdentity();
	while(std::getline(lineStream, name_pcd, '\n'))
	{
		++num_scan;
		if (num_scan < num_scan_start)
			continue;

		std::istringstream ss(name_pcd);
		std::vector<std::string> s_list(10);
		i = 0;
		while( std::getline(ss, s_list[i], '.') )
		{
			++i;
		}
		tt_pcd = std::stod(s_list[0]) - 946687000 + std::stod(s_list[1])/1000000000;

		xi_gps = inteploate(gps1.timestamp, gps1.enu.block(0,0,gps1.num_row,1), tt_pcd, true);
		yi_gps = inteploate(gps1.timestamp, gps1.enu.block(0,1,gps1.num_row,1), tt_pcd, true);
		zi_gps = inteploate(gps1.timestamp, gps1.enu.block(0,2,gps1.num_row,1), tt_pcd, true);

		xi_imu = inteploate(imu1.timestamp, imu1.quat.block(0,0,imu1.num_row,1), tt_pcd, true);
		yi_imu = inteploate(imu1.timestamp, imu1.quat.block(0,1,imu1.num_row,1), tt_pcd, true);
		zi_imu = inteploate(imu1.timestamp, imu1.quat.block(0,2,imu1.num_row,1), tt_pcd, true);
		wi_imu = inteploate(imu1.timestamp, imu1.quat.block(0,3,imu1.num_row,1), tt_pcd, true);

		rotm1.setIdentity();
		rotm1.block(0,0,3,3) = Quat2Rotm(xi_imu, yi_imu, zi_imu, wi_imu);

		rotm2.setIdentity();
		rotm2.block(0,0,3,3) = Axang2Rotm(z_axis(0), z_axis(1), z_axis(2), A*(num_scan-1));

		RRi = rotm1 * T_sensor2FLU * rotm2;
		RRi(0,3) = RRi(0,3) + xi_gps;
		RRi(1,3) = RRi(1,3) + yi_gps;
		RRi(2,3) = RRi(2,3) + zi_gps;

		path2cloud = Path2pcd + name_pcd;
		pcl::io::loadPCDFile (path2cloud, *cloud);
		pcl::transformPointCloud (*cloud, *cloud_trans, RRi);

		//==================Outputs the point cloud for debugging======================
		//name_out = "rotm1_" + std::to_string(static_cast<long long>(num_scan)) + ".txt";
		//name_out = "C:/Users/yanyujie0824/Documents/PHD_workspace/cpp_mapping/cpp_test2_10spr_1200RPM/cpp_map_v2p_point2point/rotm1/" + name_out;
		//file1.open(name_out);
		//if (file1.is_open())
		//    file1 << rotm1 << '\n';
		//file1.close();

		//name_out = "rotm2_" + std::to_string(static_cast<long long>(num_scan)) + ".txt";
		//name_out = "C:/Users/yanyujie0824/Documents/PHD_workspace/cpp_mapping/cpp_test2_10spr_1200RPM/cpp_map_v2p_point2point/rotm2/" + name_out;
		//file1.open(name_out);
		//if (file1.is_open())
		//    file1 << rotm2 << '\n';
		//file1.close();

		//name_out = "RRi_" + std::to_string(static_cast<long long>(num_scan)) + ".txt";
		//name_out = "C:/Users/yanyujie0824/Documents/PHD_workspace/cpp_mapping/cpp_test2_10spr_1200RPM/cpp_map_v2p_point2point/RRi/" + name_out;
		//file1.open(name_out);
		//if (file1.is_open())
		//    file1 << RRi << '\n';
		//file1.close();
		//=============================End===========================================

		if (nn==0)
		    *map_local = *cloud_trans;
		else
			*map_local += *cloud_trans;

		++nn;
		if (nn != 10)
			continue;

        sor.setInputCloud(map_local);
        sor.filter (*map_local_down);
		pcl::transformPointCloud (*map_local_down, *map_local_down_trans, RR_init);

        icp.setInputSource(map_local_down_trans);
        icp.setInputTarget(map_down);
        icp.align(*map_local_down_icp);

		*map = *map_down + *map_local_down_icp;
		RR_init = icp.getFinalTransformation().cast<double>() * RR_init;

		sor.setInputCloud(map);
		sor.filter (*map_down);
		nn = 0;

		//==================Outputs the point cloud for debugging======================
		//name_out = "RR_init_" + std::to_string(static_cast<long long>(num_scan)) + ".txt";
		//name_out = "C:/Users/yanyujie0824/Documents/PHD_workspace/cpp_mapping/cpp_test2_10spr_1200RPM/cpp_map_v2p_point2point/RR_init/" + name_out;
		//file1.open(name_out);
		//if (file1.is_open())
		//    file1 << RR_init << '\n';
		//file1.close();

		//name_out = "local_" + std::to_string(static_cast<long long>(num_scan)) + ".pcd";
		//name_out = "C:/Users/yanyujie0824/Documents/PHD_workspace/cpp_mapping/cpp_test2_10spr_1200RPM/cpp_map_v2p_point2point/0run_RR/" + name_out;
		//pcl::io::savePCDFile (name_out, *map_local_down_trans, true);

		//name_out = "local_" + std::to_string(static_cast<long long>(num_scan)) + ".pcd";
		//name_out = "C:/Users/yanyujie0824/Documents/PHD_workspace/cpp_mapping/cpp_test2_10spr_1200RPM/cpp_map_v2p_point2point/0run/" + name_out;
		//pcl::io::savePCDFile (name_out, *map_local_down, true);

		//name_out = "map_" + std::to_string(static_cast<long long>(num_scan)) + ".pcd";
		//name_out = "C:/Users/yanyujie0824/Documents/PHD_workspace/cpp_mapping/cpp_test2_10spr_1200RPM/cpp_map_v2p_point2point/map/" + name_out;
                //name_out = "C:/Users/yanyujie0824/Documents/PHD_workspace/cpp_mapping/cpp_test2_10spr_1200RPM/cpp_map_v2p_point2point/map/" + name_out;
		//pcl::io::savePCDFile (name_out, *map_down, true);
		//=============================End===========================================
	}
	name_out = "OfflineMap.pcd";
	name_out = "/home/ubuntu/offline_mapping/output_map/" + name_out;
	pcl::io::savePCDFile (name_out, *map, true);
	return 0;
}
