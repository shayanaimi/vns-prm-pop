/*
 * VNSVNSGOPPath.cpp
 *
 *  Created on: 22. 3. 2016
 *      Author: Robert Penicka
 */

#include "vns_obst_op_solver.h"
#include <typeinfo>

using crl::logger;

using namespace op;
using namespace crl;
using namespace crl::gui;
using namespace opendubins;

namespace op {

template<>
void VnsPrmOPSolver<HeapPoint2D>::save_prm_points(std::string filename) {
	INFO("save_prm_point_path 2d");

	std::string dir = output;
	std::string file = getOutputIterPath(filename, dir);
	assert_io(createDirectory(dir), "Can not create file in path'" + file + "'");
	std::ofstream out(file.c_str());
	assert_io(out.good(), "Cannot create path '" + file + "'");
	std::string delimiter = ", ";

	if (out.is_open()) {
		std::vector<HeapNode<HeapPoint2D>*> points = prm->get_points();
		for (int var = 0; var < points.size(); ++var) {
			out << points[var]->data.x << delimiter << points[var]->data.y << std::endl;
		}
		out.close();
	} else {
		std::cerr << "Cannot open " << filename << std::endl;
	}
}

template<>
void VnsPrmOPSolver<HeapPoint2DHeading>::save_prm_points(std::string filename) {
	INFO("save_prm_point_path 2d heading");

	std::string dir = output;
	std::string file = getOutputIterPath(filename, dir);
	assert_io(createDirectory(dir), "Can not create file in path'" + file + "'");
	std::ofstream out(file.c_str());
	assert_io(out.good(), "Cannot create path '" + file + "'");
	std::string delimiter = ", ";

	if (out.is_open()) {
		std::vector<HeapNode<HeapPoint2DHeading>*> points = prm->get_points();
		for (int var = 0; var < points.size(); ++var) {
			out << points[var]->data.x << delimiter << points[var]->data.y << delimiter << points[var]->data.phi
					<< std::endl;
		}
		out.close();
	} else {
		std::cerr << "Cannot open " << filename << std::endl;
	}
}

template<>
void VnsPrmOPSolver<HeapPoint3D>::save_prm_points(std::string filename) {
	INFO("save_prm_point_path 2d");

	std::string dir = output;
	std::string file = getOutputIterPath(filename, dir);
	assert_io(createDirectory(dir), "Can not create file in path'" + file + "'");
	std::ofstream out(file.c_str());
	assert_io(out.good(), "Cannot create path '" + file + "'");
	std::string delimiter = ", ";

	if (out.is_open()) {
		std::vector<HeapNode<HeapPoint3D>*> points = prm->get_points();
		for (int var = 0; var < points.size(); ++var) {
			out << points[var]->data.x << delimiter << points[var]->data.y << delimiter << points[var]->data.z
					<< std::endl;
		}
		out.close();
	} else {
		std::cerr << "Cannot open " << filename << std::endl;
	}
}

//not sure why we commented this out, but it works this way and doesn't compile otherwise.
template<>
void VnsPrmOPSolver<HeapPoint2D>::save_prm_point_path(std::string filename, VnsSopPath<HeapPoint2D> * toShow) {
	// //INFO_VAR(returnedPath.size());
	// VnsSopPath<HeapPoint2D> path_to_save = this->tourVNSGOPPath;
	// if (toShow != NULL) {
	// 	path_to_save = *toShow;
	// }
	// std::vector<IndexSOP> returnedPath = path_to_save.getPath();

	// if (returnedPath.size() >= 2) {
	// 	std::string dir = output;
	// 	std::string file = getOutputIterPath(filename, dir);
	// 	assert_io(createDirectory(dir), "Can not create file in path'" + file + "'");
	// 	std::ofstream out(file.c_str());
	// 	assert_io(out.good(), "Cannot create path '" + file + "'");
	// 	std::string delimiter = ", ";
	// 	if (out.is_open()) {
	// 		std::vector<HeapNode<HeapPoint2D>> plan_all;
	// 		for (int var = 1; var < returnedPath.size(); ++var) {
	// 			int node_id_nodesAllClustersfrom =
	// 					[returnedPath[var - 1].clusterIndex][returnedPath[var - 1].nodeIndex].id;
	// 			int node_id_to = nodesAllClusters[returnedPath[var].clusterIndex][returnedPath[var].nodeIndex].id;
	// 			plan_with_length<HeapPoint2D> found_path = prm->plan(node_id_from, node_id_to);
	// 			if (var == 1) {
	// 				plan_all.insert(std::end(plan_all), std::begin(found_path.plan), std::end(found_path.plan));
	// 			} else {
	// 				plan_all.insert(std::end(plan_all), std::begin(found_path.plan) + 1, std::end(found_path.plan));
	// 			}
	// 		}
	// 		for (int var = 0; var < plan_all.size(); ++var) {
	// 			out << plan_all[var].data.x << delimiter << plan_all[var].data.y << std::endl;
	// 		}
	// 		out.close();
	// 	} else {
	// 		std::cerr << "Cannot open " << filename << std::endl;
	// 	}

	// }
}

template<>
void VnsPrmOPSolver<HeapPoint2DHeading>::save_prm_point_path(std::string filename,
		VnsSopPath<HeapPoint2DHeading> * toShow) {
	VnsSopPath<HeapPoint2DHeading> path_to_save = this->tourVNSGOPPath;
	if (toShow != NULL) {
		path_to_save = *toShow;
	}
	std::vector<IndexSOP> returnedPath = path_to_save.getPath();
	//INFO_VAR(returnedPath.size());
	if (returnedPath.size() >= 2) {
		std::string dir = output;
		std::string file = getOutputIterPath(filename, dir);
		assert_io(createDirectory(dir), "Can not create file in path'" + file + "'");
		std::ofstream out(file.c_str());
		assert_io(out.good(), "Cannot create path '" + file + "'");
		std::string delimiter = ", ";
		if (out.is_open()) {
			std::vector<HeapNode<HeapPoint2DHeading>> plan_all;
			for (int var = 1; var < returnedPath.size(); ++var) {
				int node_id_from =
						nodesAllClusters[returnedPath[var - 1].clusterIndex][returnedPath[var - 1].nodeIndex].id;
				int node_id_to = nodesAllClusters[returnedPath[var].clusterIndex][returnedPath[var].nodeIndex].id;
				plan_with_length<HeapPoint2DHeading> found_path = prm->plan(node_id_from, node_id_to);
				if (var == 1) {
					plan_all.insert(std::end(plan_all), std::begin(found_path.plan), std::end(found_path.plan));
				} else {
					plan_all.insert(std::end(plan_all), std::begin(found_path.plan) + 1, std::end(found_path.plan));
				}
			}
			for (int var = 0; var < plan_all.size(); ++var) {
				out << plan_all[var].data.x << delimiter << plan_all[var].data.y << delimiter << plan_all[var].data.phi
						<< std::endl;
			}
			out.close();
		} else {
			std::cerr << "Cannot open " << filename << std::endl;
		}

	}
}

template<>
void VnsPrmOPSolver<HeapPoint3D>::save_prm_point_path(std::string filename, VnsSopPath<HeapPoint3D> * toShow) {
	//INFO_VAR(returnedPath.size());
	VnsSopPath<HeapPoint3D> path_to_save = this->tourVNSGOPPath;
	if (toShow != NULL) {
		path_to_save = *toShow;
	}
	std::vector<IndexSOP> returnedPath = path_to_save.getPath();

	if (returnedPath.size() >= 2) {
		std::string dir = output;
		std::string file = getOutputIterPath(filename, dir);
		assert_io(createDirectory(dir), "Can not create file in path'" + file + "'");
		std::ofstream out(file.c_str());
		assert_io(out.good(), "Cannot create path '" + file + "'");
		std::string delimiter = ", ";
		if (out.is_open()) {
			std::vector<HeapNode<HeapPoint3D>> plan_all;
			for (int var = 1; var < returnedPath.size(); ++var) {
				int node_id_from =
						nodesAllClusters[returnedPath[var - 1].clusterIndex][returnedPath[var - 1].nodeIndex].id;
				int node_id_to = nodesAllClusters[returnedPath[var].clusterIndex][returnedPath[var].nodeIndex].id;
				plan_with_length<HeapPoint3D> found_path = prm->plan(node_id_from, node_id_to);
				if (var == 1) {
					plan_all.insert(std::end(plan_all), std::begin(found_path.plan), std::end(found_path.plan));
				} else {
					plan_all.insert(std::end(plan_all), std::begin(found_path.plan) + 1, std::end(found_path.plan));
				}
			}
			for (int var = 0; var < plan_all.size(); ++var) {
				out << plan_all[var].data.x << delimiter << plan_all[var].data.y << delimiter << plan_all[var].data.z << std::endl;
			}
			out.close();
		} else {
			std::cerr << "Cannot open " << filename << std::endl;
		}

	}
}

template<>
void VnsPrmOPSolver<HeapPoint2D>::fillCityNodes(OP_Prolem<HeapPoint3D> &problem) {
	INFO("fillCityNodes begin");
	nodesAllClusters.resize(problem.samples.size());
	//this->nodesAllClusters = problem.samples;
	for (int clid = 0; clid < nodesAllClusters.size(); ++clid) {
		nodesAllClusters[clid].resize(problem.samples[clid].size());
		for (int nodeid = 0; nodeid < nodesAllClusters[clid].size(); ++nodeid) {
			nodesAllClusters[clid][nodeid].cluster_id = problem.samples[clid][0].cluster_id;
			nodesAllClusters[clid][nodeid].id = problem.samples[clid][0].id;
			nodesAllClusters[clid][nodeid].reward = problem.samples[clid][0].reward;
			nodesAllClusters[clid][nodeid].data.x = problem.samples[clid][0].data.x;
			nodesAllClusters[clid][nodeid].data.y = problem.samples[clid][0].data.y;
		}
		maximalRewardAll += nodesAllClusters[clid][0].reward;
		cluster_rewards[clid] = nodesAllClusters[clid][0].reward;
	}
	INFO("nodesAllClusters filled")

	cities_nodes.reserve(nodesAllClusters.size());
	for (int var = 0; var < nodesAllClusters.size(); ++var) {
		HeapNode<HeapPoint2D> *city_node = new HeapNode<HeapPoint2D>();
		city_node->data.x = nodesAllClusters[var][0].data.x;
		city_node->data.y = nodesAllClusters[var][0].data.y;
		city_node->node_id = nodesAllClusters[var][0].id;
		city_node->cluster_id = nodesAllClusters[var][0].cluster_id;
		city_node->city_node = true;
		cities_nodes.push_back(city_node);
	}
	INFO("city nodes set");
	INFO("fillCityNodes end");
}



template<>
	bool VnsPrmOPSolver<HeapPoint2DHeading>::testCollision(std::vector<MeshObject*> obstacles, MeshObject* object, HeapNode<HeapPoint2DHeading>* node) {
		Position3D object_position;
		object_position.x = node->data.x;
		object_position.y = node->data.y;
		object_position.z = 0;
		object_position.yaw = 0;
		object_position.pitch = 0;
		object_position.roll = 0;
		return MeshObject::collide(&obstacles, object, object_position);
	}

template<>
void VnsPrmOPSolver<HeapPoint2DHeading>::fillCityNodes(OP_Prolem<HeapPoint3D> &problem) {
	//this function is where most of shaya's edits were made
	INFO("fillCityNodes begin");
	INFO_VAR(dubins_resolution);
	double temp_radius = 30; //radius to be used when generating neighborhood samples
	int act_node_id = 0;
	nodesAllClusters.resize(num_clusters);

	for (int clid = 0; clid < nodesAllClusters.size(); ++clid) {
		//nodesAllClusters
		nodesAllClusters[clid].resize(dubins_resolution*neighborhood_places);
		cities_nodes.reserve(nodesAllClusters.size() * dubins_resolution * neighborhood_places);
		int id_within_cluster = 0;
		maximalRewardAll += problem.samples[clid][0].reward;
		cluster_rewards[clid] = problem.samples[clid][0].reward;
		double angle_increment = M_2PI / neighborhood_places; //the angle by which we'll increment to discretize points around the hood
		for (int circle_spot = 0; circle_spot < neighborhood_places; ++circle_spot){ //iterate through different positions around the circle	
			double spot_angle = angle_increment * circle_spot; //the angle associated with this "neighbor"
			double spot_x = (problem.samples[clid][0].data.x) + neighborhood_radius * cos(spot_angle); //basic trig to determine x and y coordinates based on angle
			double spot_y = (problem.samples[clid][0].data.y) + neighborhood_radius * sin(spot_angle);

			// Before adding the new point to the system, check for collisions
			HeapNode<HeapPoint2DHeading>* city_node = new HeapNode<HeapPoint2DHeading>(); //temporary city node just for testing purposes
			city_node->data.x = spot_x;
			city_node->data.y = spot_y;
			if (!testCollision(mesh_obstacles, mesh_robot, city_node)) {
				//if the node is okay, it needs to be added to two data structures; nodesallclusters and cities_nodes
				delete city_node;
				for (int nodeid = 0; nodeid < dubins_resolution; ++nodeid) { //iterate through different headings in this circle spot
					//first add the node to nodesallclusters
					nodesAllClusters[clid][id_within_cluster].cluster_id = problem.samples[clid][0].cluster_id;
					nodesAllClusters[clid][id_within_cluster].id = act_node_id;
					act_node_id += 1;
					nodesAllClusters[clid][id_within_cluster].data.x = spot_x;
					nodesAllClusters[clid][id_within_cluster].data.y = spot_y;
					nodesAllClusters[clid][id_within_cluster].reward = problem.samples[clid][0].reward;
					nodesAllClusters[clid][id_within_cluster].data.phi = nodeid * M_2PI / dubins_resolution;

					//now setting city node data
					HeapNode<HeapPoint2DHeading>* city_node = new HeapNode<HeapPoint2DHeading>();
					city_node->data.x = spot_x;
					city_node->data.y = spot_y;
					city_node->data.phi = nodesAllClusters[clid][id_within_cluster].data.phi;
					city_node->data.radius = dubins_radius;
					city_node->node_id = nodesAllClusters[clid][id_within_cluster].id;
					city_node->cluster_id = nodesAllClusters[clid][id_within_cluster].cluster_id;
					city_node->city_node = true;
					cities_nodes.push_back(city_node);
					id_within_cluster++;

				}
			}
			else {
				//if there is a collision, nodesallclusters needs to be resized because memory was already reserved
				delete city_node;
				nodesAllClusters[clid].resize(nodesAllClusters[clid].size() - dubins_resolution);
			}				
		}
	}
	INFO("nodesAllClusters filled");

	num_clusters = nodesAllClusters.size();
	
	
	INFO("city nodes set");
	this->prm->set_num_headings(dubins_resolution * neighborhood_places); 
	INFO("fillCityNodes end");
}

template<>
void VnsPrmOPSolver<HeapPoint3D>::fillCityNodes(OP_Prolem<HeapPoint3D> &problem) {
	INFO("fillCityNodes begin");
	this->nodesAllClusters = problem.samples;
	for (int clid = 0; clid < nodesAllClusters.size(); ++clid) {
		maximalRewardAll += nodesAllClusters[clid][0].reward;
		cluster_rewards[clid] = nodesAllClusters[clid][0].reward;
	}
	INFO("nodesAllClusters filled")

	cities_nodes.reserve(nodesAllClusters.size());
	for (int var = 0; var < nodesAllClusters.size(); ++var) {
		HeapNode<HeapPoint3D> *city_node = new HeapNode<HeapPoint3D>();
		city_node->data.x = nodesAllClusters[var][0].data.x;
		city_node->data.y = nodesAllClusters[var][0].data.y;
		city_node->data.z = nodesAllClusters[var][0].data.z;
		city_node->node_id = nodesAllClusters[var][0].id;
		city_node->cluster_id = nodesAllClusters[var][0].cluster_id;
		city_node->city_node = true;
		cities_nodes.push_back(city_node);
	}
	INFO("city nodes set");
	INFO("fillCityNodes end");
}

template<>
void VnsPrmOPSolver<HeapPoint2D>::drawVisibilityPath(int usleepTime, std::vector<HeapNode<HeapPoint2D>> * toShow) {
	if (canvas) {
		//INFO("clear path");
		*canvas << canvas::CLEAR << "path" << "path";
		std::vector<HeapNode<HeapPoint2D>> returnedPath = *toShow;
		CShape greenLine("green", "green", 3, 3);
		for (int var = 1; var < returnedPath.size(); ++var) {
			*canvas << canvas::LINE << greenLine;
			*canvas << returnedPath[var - 1].data.x << returnedPath[var - 1].data.y;
			*canvas << returnedPath[var].data.x << returnedPath[var].data.y;

		}
		*canvas << canvas::END;
		canvas->redraw();
	}
}

template<>
void VnsPrmOPSolver<HeapPoint2DHeading>::drawVisibilityPath(int usleepTime,
		std::vector<HeapNode<HeapPoint2DHeading>> * toShow) {
	if (canvas) {
		//INFO("drawVisibilityPath dubins begin");
		*canvas << canvas::CLEAR << "path" << "path";
		std::vector<HeapNode<HeapPoint2DHeading>> returnedPath = *toShow;

		CShape greenLine("green", "green", 3, 3);
		State previous_state;
		for (int var = 1; var < returnedPath.size(); ++var) {

			State start_state(returnedPath[var - 1].data.x, returnedPath[var - 1].data.y,
					returnedPath[var - 1].data.phi);
			State stop_state(returnedPath[var].data.x, returnedPath[var].data.y, returnedPath[var].data.phi);
			//INFO_VAR(start_state)
			//INFO_VAR(stop_state)
			if (var > 1) {
				if (previous_state.point.x != start_state.point.x || previous_state.point.y != start_state.point.y
						|| previous_state.ang != start_state.ang) {
					ERROR("not continuous points");
					exit(1);
				}
			}
			previous_state = stop_state;
			//INFO("dubins from "<<start_state<<" radius "<<returnedPath[var - 1].data.radius<<" to "<<stop_state<<" radius "<<returnedPath[var].data.radius)
			double radius = returnedPath[var - 1].data.radius;
			Dubins dub(start_state, stop_state, radius);
			*canvas << greenLine << dub;

		}
		*canvas << canvas::END;
		canvas->redraw();
		//INFO("drawVisibilityPath dubins end");
	}
}

template<>
void VnsPrmOPSolver<HeapPoint3D>::drawVisibilityPath(int usleepTime, std::vector<HeapNode<HeapPoint3D>> * toShow) {
	if (canvas) {
		//INFO("clear path");
		*canvas << canvas::CLEAR << "path" << "path";
		std::vector<HeapNode<HeapPoint3D>> returnedPath = *toShow;
		CShape greenLine("green", "green", 3, 3);
		for (int var = 1; var < returnedPath.size(); ++var) {
			*canvas << canvas::LINE << greenLine;
			*canvas << returnedPath[var - 1].data.x << returnedPath[var - 1].data.y;
			*canvas << returnedPath[var].data.x << returnedPath[var].data.y;

		}
		*canvas << canvas::END;
		canvas->redraw();
	}
}

template<>
std::vector<std::vector<double>> VnsPrmOPSolver<HeapPoint2D>::getTrajectoryPointVectors(
		std::vector<HeapPoint2D> samples_traj) {
	std::vector<std::vector<double>> vec;
	for (int var1 = 0; var1 < samples_traj.size(); ++var1) {
		std::vector<double> single_vec = { samples_traj[var1].x, samples_traj[var1].y, this->fly_altitude, 0 };
		vec.push_back(single_vec);
	}
	return vec;
}

template<>
std::vector<std::vector<double>> VnsPrmOPSolver<HeapPoint2DHeading>::getTrajectoryPointVectors(
		std::vector<HeapPoint2DHeading> samples_traj) {
	std::vector<std::vector<double>> vec;
	for (int var1 = 0; var1 < samples_traj.size(); ++var1) {
		std::vector<double> single_vec = { samples_traj[var1].x, samples_traj[var1].y, this->fly_altitude,
				samples_traj[var1].phi };
		vec.push_back(single_vec);
	}
	return vec;
}

template<>
std::vector<std::vector<double>> VnsPrmOPSolver<HeapPoint3D>::getTrajectoryPointVectors(
		std::vector<HeapPoint3D> samples_traj) {
	std::vector<std::vector<double>> vec;
	for (int var1 = 0; var1 < samples_traj.size(); ++var1) {
		std::vector<double> single_vec = { samples_traj[var1].x, samples_traj[var1].y, this->fly_altitude, 0 };
		vec.push_back(single_vec);
	}
	return vec;
}

}
