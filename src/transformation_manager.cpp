#include <iostream>
#include <queue>

#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>

#include "transformation_manager.h"
#include "input_output_manager.h"
#include "common_data_structures_and_functions.h"
#include "partial_geometry.h"

using namespace std;
using namespace ROOT::Math;

transformation_manager::transformation_manager(partial_geometry* a_geometry)
{
  the_geometry = a_geometry;
}

bool transformation_manager::transform_to_global_tracking_frame(vector<double> clusters[], vector<line_3d> cluster_lines[])
{
  XYZVector temp_vector;
  double dummy_double;
  
  line_3d temp_strip;
  
  RotationX rx;
  RotationY ry;
  RotationZ rz;
  
  Rotation3D total_rotation;
  
  //const double to_radians = Pi() / 180.;
  
  for (int i = 0; i < the_geometry->get_number_of_plaques(); i++)
  {
    rx.SetAngle(the_geometry->get_plaque(i)->rotation.x() * to_radians);
    ry.SetAngle(the_geometry->get_plaque(i)->rotation.y() * to_radians);
    rz.SetAngle(the_geometry->get_plaque(i)->rotation.z() * to_radians);
    
    total_rotation = rz * ry * rx;
    
    for (int j = 0; j < clusters[i].size(); j++)
    {
      // dummy_double = (clusters[i].at(j) - (((double)(the_geometry->get_plate(i)->number_of_strips)) - 1.) / 2.) * 60 * um;
      dummy_double = (clusters[i].at(j) - (((double)(the_geometry->get_plaque(i)->number_of_strips)) - 1.) / 2.) * the_geometry->get_plaque(i)->pitch;
      
      temp_vector = dummy_double * the_geometry->get_primary_direction_1();

      temp_vector = total_rotation * temp_vector;
      temp_vector = temp_vector + the_geometry->get_plaque(i)->position;
      temp_strip.point = temp_vector;
      
      // this is already done in partial_geometry; use it
      temp_vector = the_geometry->get_primary_direction_2();
      temp_vector = total_rotation * temp_vector;
      temp_strip.direction = temp_vector;
      cluster_lines[i].push_back(temp_strip);
    }
  }
  
  return true;
}

bool transformation_manager::merge_plaques_in_same_plate(vector<line_3d> plaque_cluster_lines[], vector<line_3d> cluster_lines[])
{
  int current_element, next_element;
  int current_plate, next_plate;
  queue<int> cluster_queue;
  
  int number_of_plates = 0;
  
  for (int i = 0; i < max_number_of_plates; i++)
  {
    cluster_lines[i].clear();
  }
  
  for (int i = 0; i < the_geometry->get_number_of_plaques(); i++)
  {
    cluster_queue.push(i);
  }
  
  while (!cluster_queue.empty())
  {
    current_element = cluster_queue.front();
//cout << "---- " << current_element << endl;
    current_plate = the_geometry->get_plaque(current_element)->plate;
    cluster_queue.pop();
    
    cluster_lines[number_of_plates] = plaque_cluster_lines[current_element];
    
    next_element = cluster_queue.front();
    next_plate = the_geometry->get_plaque(next_element)->plate;
    
    while (next_plate == current_plate)
    {
      cluster_queue.pop();
      
      for (int i = 0; i < plaque_cluster_lines[next_element].size(); i++)
      {
        cluster_lines[number_of_plates].push_back(plaque_cluster_lines[next_element].at(i));
      }
      
      if (!cluster_queue.empty())
      {
        next_element = cluster_queue.front();
        next_plate = the_geometry->get_plaque(next_element)->plate;
      }
      else
      {
        break;
      }
    }
    
    number_of_plates++;
  }
    
  return true;
}
