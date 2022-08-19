#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>

#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>

#include "common_data_structures_and_functions.h"
#include "global_data_dispatcher.h"
#include "partial_geometry.h"

using namespace std;
using namespace ROOT::Math;

partial_geometry::partial_geometry()
:actual_number_of_plates(0), number_of_tracking_regions(0), number_of_groups(0), primary_direction_1_name(""),
    primary_direction_2_name(""), secondary_direction_1_name(""),
    secondary_direction_2_name("")
{
  cout << "Created an empty geometry." << endl;
}

int partial_geometry::get_number_of_plaques()
{
  return actual_number_of_plaques;
}

int partial_geometry::get_number_of_plates()
{
  return actual_number_of_plates;
}

int partial_geometry::get_number_of_plates_in_tracking_region(int a_region)
{
  if (a_region < number_of_tracking_regions)
  {
    return the_tracking_regions[a_region].number_of_plates;
  }
  else
  {
    return -1;
  }
}

int partial_geometry::get_number_of_tracking_regions()
{
  return number_of_tracking_regions;
}

int partial_geometry::get_number_of_groups()
{
  return number_of_groups;
}

int partial_geometry::get_number_of_groups_in_tracking_region(int a_region)
{
  if (a_region < number_of_tracking_regions)
  {
    return the_tracking_regions[a_region].number_of_groups;
  }
  else
  {
    return -1;
  }
}

int partial_geometry::get_number_of_plates_in_group(int a_group)
{
  if (a_group < number_of_groups)
  {
    return the_plate_groups[a_group].number_of_plates_in_group;
  }
  else
  {
    return -1;
  }
}

int partial_geometry::get_number_of_upstream_plates()
{
  return the_tracking_regions[0].number_of_plates;
}

int partial_geometry::get_number_of_midstream_plates()
{
  return the_tracking_regions[1].number_of_plates;
}

int partial_geometry::get_number_of_downstream_plates()
{
  return the_tracking_regions[2].number_of_plates;
}

int partial_geometry::get_number_of_midstream_plus_downstream_plates()
{
  return (the_tracking_regions[1].number_of_plates + the_tracking_regions[2].number_of_plates);
}

int partial_geometry::get_number_of_upstream_groups()
{
  return the_tracking_regions[0].number_of_groups;
}

int partial_geometry::get_number_of_midstream_groups()
{
  return the_tracking_regions[1].number_of_groups;
}

int partial_geometry::get_number_of_downstream_groups()
{
  return the_tracking_regions[2].number_of_groups;
}

int partial_geometry::get_number_of_midstream_plus_downstream_groups()
{
  return the_tracking_regions[1].number_of_groups + the_tracking_regions[2].number_of_groups;
}

target* partial_geometry::get_target()
{
  return &the_target;
}

XYZVector partial_geometry::get_primary_direction_1()
{
  return primary_direction_1;
}

XYZVector partial_geometry::get_primary_direction_2()
{
  return primary_direction_2;
}

XYZVector partial_geometry::get_secondary_direction_1()
{
  return primary_direction_2;
}

XYZVector partial_geometry::get_secondary_direction_2()
{
  return secondary_direction_2;
}

string partial_geometry::get_primary_direction_1_name()
{
  return primary_direction_1_name;
}

string partial_geometry::get_primary_direction_2_name()
{
  return primary_direction_2_name;
}

string partial_geometry::get_secondary_direction_1_name()
{
  return secondary_direction_1_name;
}

string partial_geometry::get_secondary_direction_2_name()
{
  return secondary_direction_2_name;
}

/*
tracking_region* partial_geometry::get_tracking_region(int a_region)
{
  if (a_region < number_of_tracking_regions)
  {
    return &the_tracking_regions[a_region];
  }
  else
  {
    return nullptr;
  }
}
*/

plate_group* partial_geometry::get_plate_group(int a_group)
{
  if (a_group < number_of_groups)
  {
    return &the_plate_groups[a_group];
  }
  else
  {
    return nullptr;
  }
}

ssd_plate* partial_geometry::get_plate(int a_plate_index)
{
  if (a_plate_index < actual_number_of_plates)
  {
    return &the_plates[a_plate_index];
  }
  else
  {
    return nullptr;
  }
}

ssd_plate* partial_geometry::get_plate_by_id(int a_plate_id)
{
  // implement later
  
  return nullptr;
}

ssd_plaque* partial_geometry::get_plaque(int a_plaque_index)
{
  if (a_plaque_index < actual_number_of_plaques)
  {
    return &the_plaques[a_plaque_index];
  }
  else
  {
    return nullptr;
  }
}


bool partial_geometry::add_plaque(ssd_plaque a_plaque)
{
  if (actual_number_of_plaques == max_number_of_plaques)
  {
    cout << "Requested plaque to add above maxx allowed number. Ignoring..." << endl;
    return false;
  }
  else
  {
    the_plaques[actual_number_of_plaques] = a_plaque;
    actual_number_of_plaques++;
    return true;
  }
}

bool partial_geometry::add_target(target a_target)
{
  the_target = a_target;
  
  return true;
}

bool partial_geometry::set_primary_direction_1(XYZVector a_direction)
{
  primary_direction_1 = a_direction;
  primary_direction_1 = primary_direction_1.Unit();
  
  return true;
}

bool partial_geometry::set_primary_direction_2(XYZVector a_direction)
{
  primary_direction_2 = a_direction;
  primary_direction_2 = primary_direction_2.Unit();
  
  return true;
}

bool partial_geometry::set_secondary_direction_1(XYZVector a_direction)
{
  secondary_direction_1 = a_direction;
  secondary_direction_1 = secondary_direction_1.Unit();
  
  return true;
}

bool partial_geometry::set_secondary_direction_2(XYZVector a_direction)
{
  secondary_direction_2 = a_direction;
  secondary_direction_2 = secondary_direction_2.Unit();
  
  return true;
}

bool partial_geometry::set_primary_direction_1_name(string a_name)
{
  primary_direction_1_name = a_name;
  
  return true;
}

bool partial_geometry::set_primary_direction_2_name(string a_name)
{
  primary_direction_2_name = a_name;
  
  return true;
}

bool partial_geometry::set_secondary_direction_1_name(string a_name)
{
  secondary_direction_1_name = a_name;
  
  return true;
}

bool partial_geometry::set_secondary_direction_2_name(string a_name)
{
  secondary_direction_2_name = a_name;
  
  return true;
}

bool partial_geometry::set_tracking_region_names(int a_region, std::string a_name)
{
  if (a_region < max_number_of_tracking_regions)
  {
    tracking_region_names[a_region] = a_name;
    
    return true;
  }
  else
  {
    return false;
  }
}

bool partial_geometry::classify_plates()
{
  bool classification_succeeded = true;
  
  bool matching_plaques_in_plate;
  
  bool good_cosangle;
  
  int plaque_index;
  
  vector<int> plaques_in_plates[max_number_of_plates];
  vector<int> plates_in_groups[max_number_of_groups];
  
  string dummy_string;
  
  map<int, int> group_index_map;
  group_index_map.clear();
  
  vector<int> group_index_vector;
  group_index_vector.clear();
  
  int scales_counter;
  
  XYZVector sum_coordinates;
  
  for (int i = 0; i < max_number_of_plates; i++)
  {
    plaques_in_plates[i].clear();
  }
  
  for (int i = 0; i < max_number_of_groups; i++)
  {
    plates_in_groups[i].clear();
  }
  
  for (int i = 0; i < actual_number_of_plaques; i++)
  {
    plaques_in_plates[the_plaques[i].plate].push_back(i);
    
    the_plaques[i].pitch = the_plaques[i].size.x() / ((double) the_plaques[i].number_of_strips);
  }
  
  for (int i = 0; i < max_number_of_tracking_regions; i++)
  {
    the_tracking_regions[i].number_of_groups = 0;
    the_tracking_regions[i].number_of_plates = 0;
  }
  
  for (int i = 0; i < max_number_of_tracking_regions; i++)
  {
    if (tracking_region_names[i] != "NA")
    {
      number_of_tracking_regions++;
    }
  }
  
  for (int i = 0; i < number_of_tracking_regions; i++)
  {
    the_tracking_regions[i].name = tracking_region_names[i];
  }
  
  actual_number_of_plates = 0;
  
  for (int i = 0; i < max_number_of_plates; i++)
  {
    if (plaques_in_plates[i].size() != 0)
    {
      actual_number_of_plates++;
      
      the_plates[i].number_of_member_plaques = plaques_in_plates[i].size();
      
      for (int j = 0; j < plaques_in_plates[i].size(); j++)
      {
        the_plates[i].member_plaques[j] = plaques_in_plates[i].at(j);
      }
    }
  }
  
  for (int i = 0; i < max_number_of_plates; i++)
  {
    sum_coordinates.SetCoordinates(0., 0., 0.);

    for (int j = 0; j < plaques_in_plates[i].size(); j++)
    {
      sum_coordinates += the_plaques[the_plates[i].member_plaques[j]].position;
    }
    
    the_plates[i].position = sum_coordinates / ((double) plaques_in_plates[i].size());
  }
  
  matching_plaques_in_plate = true;
  
cout << actual_number_of_plates << endl;
  
  for (int i = 0; i < actual_number_of_plates; i++)
  {
    if (the_plates[i].number_of_member_plaques == 1)
    {
      plaque_index = the_plates[i].member_plaques[0];
      
      the_plates[i].intended_to_measure = the_plaques[plaque_index].type;
      the_plates[i].group = the_plaques[plaque_index].group;
      the_plates[i].tracking_region = the_plaques[plaque_index].tracking_region;

      good_cosangle = check_cosangle(plaque_index);
      
      if (!good_cosangle) matching_plaques_in_plate = false;
    }
    else if (the_plates[i].number_of_member_plaques == 2)
    {
      for (int j = 0; j < 2; j++)
      {
        plaque_index = the_plates[i].member_plaques[j];
        
        good_cosangle = check_cosangle(plaque_index);
        
        if (!good_cosangle) matching_plaques_in_plate = false;
      }
      
      if (the_plaques[the_plates[i].member_plaques[0]].type == the_plaques[the_plates[i].member_plaques[1]].type)
      {
        the_plates[i].intended_to_measure = the_plaques[the_plates[i].member_plaques[0]].type;
      }
      else
      {
        matching_plaques_in_plate = false;
      }
      
      if (the_plaques[the_plates[i].member_plaques[0]].group == the_plaques[the_plates[i].member_plaques[1]].group)
      {
        the_plates[i].group = the_plaques[the_plates[i].member_plaques[0]].group;
      }
      else
      {
        matching_plaques_in_plate = false;
      }
      
      if (the_plaques[the_plates[i].member_plaques[0]].tracking_region == the_plaques[the_plates[i].member_plaques[1]].tracking_region)
      {
        the_plates[i].tracking_region = the_plaques[the_plates[i].member_plaques[0]].tracking_region;
      }
      else
      {
        matching_plaques_in_plate = false;
      }
    }
    else if (the_plates[i].number_of_member_plaques == 3)
    {
      for (int j = 0; j < 3; j++)
      {
        plaque_index = the_plates[i].member_plaques[j];
        
        good_cosangle = check_cosangle(plaque_index);
        
        if (!good_cosangle) matching_plaques_in_plate = false;
      }
      
      if (the_plaques[the_plates[i].member_plaques[0]].type == the_plaques[the_plates[i].member_plaques[1]].type &&
          the_plaques[the_plates[i].member_plaques[1]].type == the_plaques[the_plates[i].member_plaques[2]].type)
      {
        the_plates[i].intended_to_measure = the_plaques[the_plates[i].member_plaques[0]].type;
      }
      else
      {
        matching_plaques_in_plate = false;
      }
      
      if (the_plaques[the_plates[i].member_plaques[0]].group == the_plaques[the_plates[i].member_plaques[1]].group &&
          the_plaques[the_plates[i].member_plaques[1]].group == the_plaques[the_plates[i].member_plaques[2]].group)
      {
        the_plates[i].group = the_plaques[the_plates[i].member_plaques[0]].group;
      }
      else
      {
        matching_plaques_in_plate = false;
      }
      
      if (the_plaques[the_plates[i].member_plaques[0]].tracking_region == the_plaques[the_plates[i].member_plaques[1]].tracking_region &&
          the_plaques[the_plates[i].member_plaques[1]].tracking_region == the_plaques[the_plates[i].member_plaques[2]].tracking_region)
      {
        the_plates[i].tracking_region = the_plaques[the_plates[i].member_plaques[0]].tracking_region;
      }
      else
      {
        matching_plaques_in_plate = false;
      }
    }
    else
    {
      // only sets of 1, 2 or 3 plates are allowed in this phase of EMPHATIC
      cout << "Number of plaques in a plate not allowed: " << the_plates[i].number_of_member_plaques << " plaques in plate " << i << endl;
    }
  }
  
cout << "!!!! " << the_plates[12].tracking_region << endl;
  
  
  if (matching_plaques_in_plate)
  {
    for (int i = 0; i < actual_number_of_plates; i++)
    {
      if (group_index_map[the_plates[i].group] == 0)
      {
        group_index_map[the_plates[i].group]++;
        group_index_vector.push_back(the_plates[i].group);
      }
    }
    
    for (int i = 0; i < group_index_vector.size(); i++)
    {
      if (group_index_vector[i] != i)
      {
        cout << "Wrong order of plate groups (must be consequtive integers starting at 0)..." << endl;
        classification_succeeded = false;
      }
    }
  }
  else
  {
    classification_succeeded = false;
  }
  
  
  
  if (classification_succeeded)
  {
    number_of_groups = group_index_vector.size();
  
    for (int i = 0; i < max_number_of_groups; i++)
    {
      the_plate_groups[i].number_of_plates_in_group = 0;
    }
  
    for (int i = 0; i < actual_number_of_plates; i++)
    {
      if (the_plates[i].tracking_region < number_of_tracking_regions)
      {
        the_tracking_regions[the_plates[i].tracking_region].plate_index_list[the_tracking_regions[the_plates[i].tracking_region].number_of_plates] = i;
        the_tracking_regions[the_plates[i].tracking_region].number_of_plates++;
      }
      else
      {
        cout << "Invalid tracking region: " << the_plates[i].tracking_region << " for plate: " << i << endl;
        classification_succeeded = false;
      }
    
      if (the_plates[i].group < number_of_groups)
      {
        the_plate_groups[the_plates[i].group].plate_index_list[the_plate_groups[the_plates[i].group].number_of_plates_in_group] = i;
        the_plate_groups[the_plates[i].group].number_of_plates_in_group++;
      }
      else
      {
        cout << "Invalid group: " << the_plates[i].group << endl;
        classification_succeeded = false;
      }
    }
  }
  
  for (int i = 0; i < number_of_groups; i++)
  {
    if (the_plate_groups[i].number_of_plates_in_group == 3)
    {
      the_triple_groups[number_of_triple_groups] = i;
      number_of_triple_groups++;
      
      dummy_string  = the_plates[the_plate_groups[i].plate_index_list[0]].intended_to_measure;
      dummy_string += "_";
      dummy_string += the_plates[the_plate_groups[i].plate_index_list[1]].intended_to_measure;
      dummy_string += "_";
      dummy_string += the_plates[the_plate_groups[i].plate_index_list[2]].intended_to_measure;
      the_plate_groups[i].group_type = dummy_string;
      
      if (the_plates[the_plate_groups[i].plate_index_list[0]].tracking_region == the_plates[the_plate_groups[i].plate_index_list[1]].tracking_region &&
          the_plates[the_plate_groups[i].plate_index_list[0]].tracking_region == the_plates[the_plate_groups[i].plate_index_list[2]].tracking_region)
      {
        the_plate_groups[i].tracking_region = the_plates[the_plate_groups[i].plate_index_list[0]].tracking_region;
      }
      else
      {
        cout << "Plates in the same group are in different tracking regions..." << endl;
        classification_succeeded = false;
      }
    }
    else if (the_plate_groups[i].number_of_plates_in_group == 2)
    {
      dummy_string  = the_plates[the_plate_groups[i].plate_index_list[0]].intended_to_measure;
      dummy_string += "_";
      dummy_string += the_plates[the_plate_groups[i].plate_index_list[1]].intended_to_measure;
      the_plate_groups[i].group_type = dummy_string;
      
      if (the_plates[the_plate_groups[i].plate_index_list[0]].tracking_region == the_plates[the_plate_groups[i].plate_index_list[1]].tracking_region)
      {
        the_plate_groups[i].tracking_region = the_plates[the_plate_groups[i].plate_index_list[0]].tracking_region;
      }
      else
      {
        cout << "Plates in the same group are in different tracking regions..." << endl;
        classification_succeeded = false;
      }
    }
    else if (the_plate_groups[i].number_of_plates_in_group == 1)
    {
      dummy_string = the_plates[the_plate_groups[i].plate_index_list[0]].intended_to_measure;
      the_plate_groups[i].group_type = dummy_string;
      the_plate_groups[i].tracking_region = the_plates[the_plate_groups[i].plate_index_list[0]].tracking_region;
    }
    else
    {
      cout << "Number of plates in group different than 1, 2 or 3: group " << i << " with " << the_plate_groups[i].number_of_plates_in_group << " groups." << endl;
      
      classification_succeeded = false;
    }
  }
  
  for (int i = 0; i < number_of_groups; i++)
  {
    the_tracking_regions[the_plate_groups[i].tracking_region].plate_index_list[the_tracking_regions[the_plate_groups[i].tracking_region].number_of_groups] = i;
    the_tracking_regions[the_plate_groups[i].tracking_region].number_of_groups++;
  }
  
  for (int i = 0; i < number_of_groups; i++)
  {
    if (the_plate_groups[i].number_of_plates_in_group == 2)
    {
      // there can be more allowed combinations if arbitrary directions are allowed or a group can be
      // X_D1 or something like that
      if (!(the_plate_groups[i].group_type == "X_Y" || the_plate_groups[i].group_type == "Y_X" ||
            the_plate_groups[i].group_type == "D1_D2" || the_plate_groups[i].group_type == "D2_D1"))
      {
        cout << "A group of plates is found that are not compatible: group " << i << " is type " << the_plate_groups[i].group_type << endl;
        classification_succeeded = false;
      }
    }
    else if (the_plate_groups[i].number_of_plates_in_group == 3)
    {
      // ditto; so a smarter way can be found...
      if (!(the_plate_groups[i].group_type == "X_Y_D1" || the_plate_groups[i].group_type == "Y_X_D1" ||
            the_plate_groups[i].group_type == "X_D1_Y" || the_plate_groups[i].group_type == "Y_D1_X" ||
            the_plate_groups[i].group_type == "D1_X_Y" || the_plate_groups[i].group_type == "D1_Y_X" ||
            the_plate_groups[i].group_type == "X_Y_D2" || the_plate_groups[i].group_type == "Y_X_D2" ||
            the_plate_groups[i].group_type == "X_D2_Y" || the_plate_groups[i].group_type == "Y_D2_X" ||
            the_plate_groups[i].group_type == "D2_X_Y" || the_plate_groups[i].group_type == "D2_Y_X" ))
          /*
          ||
            the_plate_groups[i].group_type == "D1_D2_X" || the_plate_groups[i].group_type == "D2_D1_X" ||
            the_plate_groups[i].group_type == "D1_X_D2" || the_plate_groups[i].group_type == "D2_X_D1" ||
            the_plate_groups[i].group_type == "X_D1_D2" || the_plate_groups[i].group_type == "X_D2_D1" ||
            the_plate_groups[i].group_type == "D1_D2_Y" || the_plate_groups[i].group_type == "D2_D1_Y" ||
            the_plate_groups[i].group_type == "D1_Y_D2" || the_plate_groups[i].group_type == "D2_Y_D1" ||
            the_plate_groups[i].group_type == "Y_D1_D2" || the_plate_groups[i].group_type == "Y_D2_D1"))
           */
      {
        cout << "A group of plates is found that are not compatible: group " << i << " is type " << the_plate_groups[i].group_type << endl;
        classification_succeeded = false;
      }
    }
  }
  
  // now the other direction
  if (magnetic_field_general_direction_name != primary_direction_1_name)
  {
    cout << "Forbidden magnetic field direction. Code needs modification..." << endl;
    classification_succeeded = false;
  }

  number_of_midstream_X_plates = 0;
  number_of_midstream_Y_plates = 0;
  number_of_diagonal_plates_midstream = 0;
  number_of_d1_plates_midstream = 0;
  number_of_d2_plates_midstream = 0;
  
  if (classification_succeeded)
  {
    for (int i = the_tracking_regions[0].number_of_plates; i < the_tracking_regions[0].number_of_plates + the_tracking_regions[1].number_of_plates; i++)
    {
      if (the_plates[i].intended_to_measure == "X")
      {
        midstream_X_plates[number_of_midstream_X_plates] = i;
        number_of_midstream_X_plates++;
      }
      else if (the_plates[i].intended_to_measure == "Y")
      {
        midstream_Y_plates[number_of_midstream_Y_plates] = i;
        number_of_midstream_Y_plates++;
      }
      else if (the_plates[i].intended_to_measure == "D1")
      {
        d1_plates_midstream[number_of_d1_plates_midstream] = i;
        number_of_d1_plates_midstream++;
        
        diagonal_plates_midstream[number_of_diagonal_plates_midstream] = i;
        number_of_diagonal_plates_midstream++;
      }
      else if (the_plates[i].intended_to_measure == "D2")
      {
        d2_plates_midstream[number_of_d2_plates_midstream] = i;
        number_of_d2_plates_midstream++;
        
        diagonal_plates_midstream[number_of_diagonal_plates_midstream] = i;
        number_of_diagonal_plates_midstream++;
      }
    }
    
    // this is assuming single-plaque plates midstream; otherwise need more scales
    
    scales_counter = 0;
    
    for (int i = 0; i < number_of_midstream_X_plates - 1; i++)
    {
      for (int j = i + 1; j < number_of_midstream_X_plates; j++)
      {
        midstream_X_scales_min[scales_counter] = (the_plaques[midstream_X_plates[j]].position.Z() - (the_target.position.Z() - the_target.size.Z() / 2.)) / (the_plaques[midstream_X_plates[i]].position.Z() - (the_target.position.Z() - the_target.size.Z() / 2.));
        midstream_X_scales_max[scales_counter] = (the_plaques[midstream_X_plates[j]].position.Z() - (the_target.position.Z() + the_target.size.Z() / 2.)) / (the_plaques[midstream_X_plates[i]].position.Z() - (the_target.position.Z() + the_target.size.Z() / 2.));
        scales_counter++;
      }
    }

    scales_counter = 0;

    for (int i = 0; i < number_of_midstream_Y_plates - 1; i++)
    {
      for (int j = i + 1; j < number_of_midstream_Y_plates; j++)
      {
        midstream_Y_scales_min[scales_counter] = (the_plaques[midstream_Y_plates[j]].position.Z() - (the_target.position.Z() - the_target.size.Z() / 2.)) / (the_plaques[midstream_Y_plates[i]].position.Z() - (the_target.position.Z() - the_target.size.Z() / 2.));
        midstream_Y_scales_max[scales_counter] = (the_plaques[midstream_Y_plates[j]].position.Z() - (the_target.position.Z() + the_target.size.Z() / 2.)) / (the_plaques[midstream_Y_plates[i]].position.Z() - (the_target.position.Z() + the_target.size.Z() / 2.));
        scales_counter++;
      }
    }
  }
  
  number_of_upstream_XY_groups = 0;
  number_of_midstream_XY_groups = 0;
  number_of_downstream_XY_groups = 0;
  
  number_of_upstream_DD_groups = 0;
  number_of_midstream_DD_groups = 0;
  number_of_downstream_DD_groups = 0;

  for (int i = 0; i < number_of_groups; i++)
  {
    if (the_plate_groups[i].tracking_region == 0)
    {
      if (the_plate_groups[i].group_type == "X_Y" || the_plate_groups[i].group_type == "Y_X" ||
          the_plate_groups[i].group_type == "X_Y_D1" || the_plate_groups[i].group_type == "Y_X_D1" ||
          the_plate_groups[i].group_type == "X_D1_Y" || the_plate_groups[i].group_type == "Y_D1_X" ||
          the_plate_groups[i].group_type == "D1_X_Y" || the_plate_groups[i].group_type == "D1_Y_X" ||
          the_plate_groups[i].group_type == "X_Y_D2" || the_plate_groups[i].group_type == "Y_X_D2" ||
          the_plate_groups[i].group_type == "X_D2_Y" || the_plate_groups[i].group_type == "Y_D2_X" ||
          the_plate_groups[i].group_type == "D2_X_Y" || the_plate_groups[i].group_type == "D2_Y_X")
      {
        upstream_XY_groups[number_of_upstream_XY_groups] = i;
        number_of_upstream_XY_groups++;
      }
      else if (the_plate_groups[i].group_type == "D1_D2" || the_plate_groups[i].group_type == "D2_D1")
      {
        if (number_of_upstream_DD_groups == 0)
        {
          upstream_DD_groups[number_of_upstream_DD_groups] = i;
          number_of_upstream_DD_groups++;
        }
        else
        {
          cout << "Too many DD groups in tracking region 0. Correct the maximum value in common_data_structures_and_functions.h." << endl;
          cout << "Quitting..." << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    else if (the_plate_groups[i].tracking_region == 1)
    {
      if (the_plate_groups[i].group_type == "X_Y" || the_plate_groups[i].group_type == "Y_X" ||
          the_plate_groups[i].group_type == "X_Y_D1" || the_plate_groups[i].group_type == "Y_X_D1" ||
          the_plate_groups[i].group_type == "X_D1_Y" || the_plate_groups[i].group_type == "Y_D1_X" ||
          the_plate_groups[i].group_type == "D1_X_Y" || the_plate_groups[i].group_type == "D1_Y_X" ||
          the_plate_groups[i].group_type == "X_Y_D2" || the_plate_groups[i].group_type == "Y_X_D2" ||
          the_plate_groups[i].group_type == "X_D2_Y" || the_plate_groups[i].group_type == "Y_D2_X" ||
          the_plate_groups[i].group_type == "D2_X_Y" || the_plate_groups[i].group_type == "D2_Y_X")
      {
        midstream_XY_groups[number_of_midstream_XY_groups] = i;
        number_of_midstream_XY_groups++;
      }
      else if (the_plate_groups[i].group_type == "D1_D2" || the_plate_groups[i].group_type == "D2_D1")
      {
        if (number_of_midstream_DD_groups == 0)
        {
          midstream_DD_groups[number_of_midstream_DD_groups] = i;
          number_of_midstream_DD_groups++;
        }
        else
        {
          cout << "Too many DD groups in tracking region 1. Correct the maximum value in common_data_structures_and_functions.h." << endl;
          cout << "Quitting..." << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    else if (the_plate_groups[i].tracking_region == 2)
    {
      if (the_plate_groups[i].group_type == "X_Y" || the_plate_groups[i].group_type == "Y_X" ||
          the_plate_groups[i].group_type == "X_Y_D1" || the_plate_groups[i].group_type == "Y_X_D1" ||
          the_plate_groups[i].group_type == "X_D1_Y" || the_plate_groups[i].group_type == "Y_D1_X" ||
          the_plate_groups[i].group_type == "D1_X_Y" || the_plate_groups[i].group_type == "D1_Y_X" ||
          the_plate_groups[i].group_type == "X_Y_D2" || the_plate_groups[i].group_type == "Y_X_D2" ||
          the_plate_groups[i].group_type == "X_D2_Y" || the_plate_groups[i].group_type == "Y_D2_X" ||
          the_plate_groups[i].group_type == "D2_X_Y" || the_plate_groups[i].group_type == "D2_Y_X")
      {
        downstream_XY_groups[number_of_downstream_XY_groups] = i;
        number_of_downstream_XY_groups++;
      }
      else if (the_plate_groups[i].group_type == "D1_D2" || the_plate_groups[i].group_type == "D2_D1")
      {
        if (number_of_downstream_DD_groups == 0)
        {
          downstream_DD_groups[number_of_downstream_DD_groups] = i;
          number_of_downstream_DD_groups++;
        }
        else
        {
          cout << "Too many DD groups in tracking region 2. Correct the maximum value in common_data_structures_and_functions.h." << endl;
          cout << "Quitting..." << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    else
    {
      cout << "Unknown tracking region: " << the_plate_groups[i].tracking_region << endl;
      classification_succeeded = false;
    }
  }
  
  if (DEBUG_PRINT_LEVEL == 3)
  {
    cout << "Geometry directions:" << endl;
    cout << primary_direction_1 << endl;
    cout << primary_direction_1_name << endl;
    cout << primary_direction_2 << endl;
    cout << primary_direction_2_name << endl;
    cout << secondary_direction_1 << endl;
    cout << secondary_direction_1_name << endl;
    cout << secondary_direction_2 << endl;
    cout << secondary_direction_2_name << endl;
  
    cout << "-----------" << endl;
    cout << "-----------" << endl;

    cout << "Tracking regions: " << endl;
    
    for (int i = 0; i < number_of_tracking_regions; i++)
    {
      cout << tracking_region_names[i] << endl;
    }
  
    cout << "-----------" << endl;
    cout << "-----------" << endl;
  
    cout << "Groups: " << number_of_groups << endl;
  
    for (int i = 0; i < number_of_groups; i++)
    {
      cout << the_plate_groups[i].number_of_plates_in_group << " " << the_plate_groups[i].group_type << " " << the_plate_groups[i].tracking_region << endl;
    
      for (int j = 0; j < the_plate_groups[i].number_of_plates_in_group; j++)
      {
        cout << the_plate_groups[i].plate_index_list[j] << endl;
      }
    }
    
    cout << "XY groups: " << number_of_upstream_XY_groups + number_of_midstream_XY_groups + number_of_downstream_XY_groups << endl;
    cout << "upstream: ";
    for (int i = 0; i < number_of_upstream_XY_groups; i++)
    {
      cout << upstream_XY_groups[i] << " ";
    }
    cout << endl;
    cout << "midstream: ";
    for (int i = 0; i < number_of_midstream_XY_groups; i++)
    {
      cout << midstream_XY_groups[i] << " ";
    }
    cout << endl;
    cout << "downstream: ";
    for (int i = 0; i < number_of_downstream_XY_groups; i++)
    {
      cout << downstream_XY_groups[i] << " ";
    }
    cout << endl;
    
    cout << "DD groups: " << number_of_upstream_DD_groups + number_of_midstream_DD_groups + number_of_downstream_DD_groups << endl;
    cout << "upstream: ";
    for (int i = 0; i < number_of_upstream_DD_groups; i++)
    {
      cout << upstream_DD_groups[i] << " ";
    }
    cout << endl;
    cout << "midstream: ";
    for (int i = 0; i < number_of_midstream_DD_groups; i++)
    {
      cout << midstream_DD_groups[i] << " ";
    }
    cout << endl;
    cout << "downstream: ";
    for (int i = 0; i < number_of_downstream_DD_groups; i++)
    {
      cout << downstream_DD_groups[i] << " ";
    }
    cout << endl;
    
    cout << "Triple groups: " << number_of_triple_groups << endl;
    for (int i = 0; i < number_of_triple_groups; i++)
    {
      cout << the_triple_groups[i] << " ";
    }
    cout << endl;
  
    cout << "-----------" << endl;
    cout << "-----------" << endl;
  
    /*
    
    cout << "Plaques: " << actual_number_of_plaques << endl;
  
    for (int i = 0; i < actual_number_of_plaques; i++)
    {
      cout << "-----------" << endl;
      cout << "plaque " << i << ":" << endl;
      cout << the_plaques[i].plaque_id << endl;
      cout << the_plaques[i].tracking_region << endl;
      cout << the_plaques[i].plate << endl;
      cout << the_plaques[i].group << endl;
      cout << the_plaques[i].type << endl;
      cout << the_plaques[i].efficiency << endl;
      cout << the_plaques[i].position << endl;
      cout << the_plaques[i].rotation << endl;
      cout << the_plaques[i].size << endl;
      cout << the_plaques[i].number_of_strips << endl;
      cout << the_plaques[i].strip_direction << endl;
      cout << the_plaques[i].pitch << endl;
    }
    
    cout << "-----------" << endl;
    cout << "-----------" << endl;
    */
    cout << "Plates: " << actual_number_of_plates << endl;
    
    for (int i = 0; i < actual_number_of_plates; i++)
    {
      cout << "-----------" << endl;
      cout << "plate " << i << ":" << endl;
      cout << the_plates[i].number_of_member_plaques << endl;
      
      cout << "plaques: ";
      for (int j = 0; j < the_plates[i].number_of_member_plaques; j++)
      {
        cout << the_plates[i].member_plaques[j] << " ";
      }
      cout << endl;
      cout << the_plates[i].group << endl;
      cout << the_plates[i].tracking_region << endl;
      cout << the_plates[i].intended_to_measure << endl;
      cout << the_plates[i].position << endl;
    }
    
    cout << "-----------" << endl;
    cout << "-----------" << endl;

     
    cout << "Midstream Scales:" << endl;

    cout << "-----------" << endl;

    cout << "X min: " << endl;
    for (int i = 0; i < number_of_midstream_X_plates * (number_of_midstream_X_plates - 1) / 2; i++)
    {
      cout << midstream_X_scales_min[i] << " ";
    }
    cout << endl;
    
    cout << "X max: " << endl;
    for (int i = 0; i < number_of_midstream_X_plates * (number_of_midstream_X_plates - 1) / 2; i++)
    {
      cout << midstream_X_scales_max[i] << " ";
    }
    cout << endl;

    cout << "Y min: " << endl;
    for (int i = 0; i < number_of_midstream_Y_plates * (number_of_midstream_Y_plates - 1) / 2; i++)
    {
      cout << midstream_Y_scales_min[i] << " ";
    }
    cout << endl;

    cout << "Y max: " << endl;
    for (int i = 0; i < number_of_midstream_Y_plates * (number_of_midstream_Y_plates - 1) / 2; i++)
    {
      cout << midstream_Y_scales_max[i] << " ";
    }
    cout << endl;

    cout << "-----------" << endl;
    cout << "-----------" << endl;

    cout << "Magnetic field general direction:" << endl;
    
    cout << "-----------" << endl;

    cout << magnetic_field_general_direction << " " << magnetic_field_general_direction_name << endl;
    
    cout << "-----------" << endl;
    cout << "-----------" << endl;

    cout << "Target: " << endl;
    cout << the_target.position << endl;
    cout << the_target.rotation << endl;
    cout << the_target.size << endl;
    cout << the_target.material << endl;
    
    cout << "-----------" << endl;
    cout << "-----------" << endl;
  }
  
  
  return classification_succeeded;
}

ROOT::Math::XYZVector partial_geometry::get_magnetic_field_general_direction()
{
  return magnetic_field_general_direction;
}

bool partial_geometry::set_magnetic_field_general_direction(ROOT::Math::XYZVector a_direction)
{
  magnetic_field_general_direction = a_direction;
}

bool partial_geometry::set_magnetic_field_general_direction_name(std::string a_magnetic_field_general_direction_name)
{
  magnetic_field_general_direction_name = a_magnetic_field_general_direction_name;
}

double partial_geometry::get_target_min_z()
{
  // assuming target perfectly aligned
  return (the_target.position.Z() - the_target.size.Z() / 2.);
}

double partial_geometry::get_target_max_z()
{
  // assuming target perfectly aligned
  return (the_target.position.Z() + the_target.size.Z() / 2.);
}

bool partial_geometry::check_cosangle(int a_plaque_index)
{
  bool check_success = true;
  
  RotationX rx;
  RotationY ry;
  RotationZ rz;
  
  Rotation3D total_rotation;
  
  XYZVector dummy_vector_1, dummy_vector_2;
  
  double max_cosangle = the_global_data_dispatcher->get_reconstruction_options()->maximum_cosangle_between_intended_and_actual_strip_direction;

  rx.SetAngle(the_plaques[a_plaque_index].rotation.x() * to_radians);
  ry.SetAngle(the_plaques[a_plaque_index].rotation.y() * to_radians);
  rz.SetAngle(the_plaques[a_plaque_index].rotation.z() * to_radians);
  
  total_rotation = rz * ry * rx;
  
  // the_plates[i].strip_direction.SetCoordinates(0., 1., 0.);
  the_plaques[a_plaque_index].strip_direction = primary_direction_2;
  the_plaques[a_plaque_index].strip_direction = total_rotation * the_plaques[a_plaque_index].strip_direction;
  
  if (the_plaques[a_plaque_index].type == primary_direction_1_name)
  {
    dummy_vector_1 = primary_direction_1.Unit();
    dummy_vector_2 = the_plaques[a_plaque_index].strip_direction.Unit();
    
    if (fabs(dummy_vector_1.Dot(dummy_vector_2)) < max_cosangle)
    {
      cout << "Plate intended_to_measure and rotation strip directions do not coincide: " << a_plaque_index << endl;
      check_success = false;
    }
  }
  else if (the_plaques[a_plaque_index].type == primary_direction_2_name)
  {
    dummy_vector_1 = primary_direction_2.Unit();
    dummy_vector_2 = the_plaques[a_plaque_index].strip_direction.Unit();
    
    if (fabs(dummy_vector_1.Dot(dummy_vector_2)) < max_cosangle)
    {
      cout << "Plate intended_to_measure and rotation strip directions do not coincide: " << a_plaque_index << endl;
      check_success = false;
    }
  }
  else if (the_plaques[a_plaque_index].type == secondary_direction_1_name)
  {
    dummy_vector_1 = secondary_direction_1.Unit();
    dummy_vector_2 = the_plaques[a_plaque_index].strip_direction.Unit();
    
    if (fabs(dummy_vector_1.Dot(dummy_vector_2)) < max_cosangle)
    {
      cout << "Plate intended_to_measure and rotation strip directions do not coincide: " << a_plaque_index << endl;
      check_success = false;
    }
  }
  else if (the_plaques[a_plaque_index].type == secondary_direction_2_name)
  {
    dummy_vector_1 = secondary_direction_2.Unit();
    dummy_vector_2 = the_plaques[a_plaque_index].strip_direction.Unit();
    
    if (fabs(dummy_vector_1.Dot(dummy_vector_2)) < max_cosangle)
    {
      cout << "Plate intended_to_measure and rotation strip directions do not coincide: " << a_plaque_index << endl;
      check_success = false;
    }
  }
  else
  {
    cout << "Unknown intended_to_measure direction..." << endl;
    check_success = false;
  }
  
  the_plaques[a_plaque_index].normal_direction.SetCoordinates(0., 0., 1.);
  the_plaques[a_plaque_index].normal_direction = total_rotation * the_plaques[a_plaque_index].normal_direction;
  
  the_plaques[a_plaque_index].plaque_plane_equation.a = the_plaques[a_plaque_index].normal_direction.x();
  the_plaques[a_plaque_index].plaque_plane_equation.b = the_plaques[a_plaque_index].normal_direction.y();
  the_plaques[a_plaque_index].plaque_plane_equation.c = the_plaques[a_plaque_index].normal_direction.z();
  the_plaques[a_plaque_index].plaque_plane_equation.d = -(the_plaques[a_plaque_index].normal_direction.x() * the_plaques[a_plaque_index].position.x() +
                                           the_plaques[a_plaque_index].normal_direction.y() * the_plaques[a_plaque_index].position.y() +
                                           the_plaques[a_plaque_index].normal_direction.z() * the_plaques[a_plaque_index].position.z());
  
  return check_success;
}

bool partial_geometry::get_midstream_X_plates(int& a_number_of_midstream_X_plates, int a_midstream_X_plates[])
{
  if (number_of_midstream_X_plates > 1)
  {
    a_number_of_midstream_X_plates = number_of_midstream_X_plates;
  
    for (int i = 0; i < number_of_midstream_X_plates; i++)
    {
      a_midstream_X_plates[i] = midstream_X_plates[i];
    }
  
    return true;
  }
  else
  {
    return false;
  }
}

bool partial_geometry::get_midstream_Y_plates(int& a_number_of_midstream_Y_plates, int a_midstream_Y_plates[])
{
  if (number_of_midstream_Y_plates > 1)
  {
    a_number_of_midstream_Y_plates = number_of_midstream_Y_plates;
    
    for (int i = 0; i < number_of_midstream_Y_plates; i++)
    {
      a_midstream_Y_plates[i] = midstream_Y_plates[i];
    }
    
    return true;
  }
  else
  {
    return false;
  }
}

bool partial_geometry::get_midstream_X_scales(double a_midstream_X_scales_max[], double a_midstream_X_scales_min[])
{
  if (number_of_midstream_X_plates > 1)
  {
    for (int i = 0; i < number_of_midstream_X_plates * (number_of_midstream_X_plates - 1) / 2; i++)
    {
      a_midstream_X_scales_max[i] = midstream_X_scales_max[i];
      a_midstream_X_scales_min[i] = midstream_X_scales_min[i];
    }
    
    return true;
  }
  else
  {
    return false;
  }
}

bool partial_geometry::get_midstream_Y_scales(double a_midstream_Y_scales_max[], double a_midstream_Y_scales_min[])
{
  if (number_of_midstream_Y_plates > 1)
  {
    for (int i = 0; i < number_of_midstream_Y_plates * (number_of_midstream_Y_plates - 1) / 2; i++)
    {
      a_midstream_Y_scales_max[i] = midstream_Y_scales_max[i];
      a_midstream_Y_scales_min[i] = midstream_Y_scales_min[i];
    }
    
    return true;
  }
  else
  {
    return false;
  }
}

bool partial_geometry::get_upstream_XY_groups(int& a_number_of_upstream_XY_groups, int a_upstream_XY_groups[])
{
  a_number_of_upstream_XY_groups = number_of_upstream_XY_groups;
  
  for (int i = 0; i < number_of_upstream_XY_groups; i++)
  {
    a_upstream_XY_groups[i] = upstream_XY_groups[i];
  }
  
  return true;
}

bool partial_geometry::get_midstream_XY_groups(int& a_number_of_midstream_XY_groups, int a_midstream_XY_groups[])
// this includes triple groups
{
  a_number_of_midstream_XY_groups = number_of_midstream_XY_groups;
  
  for (int i = 0; i < number_of_midstream_XY_groups; i++)
  {
    a_midstream_XY_groups[i] = midstream_XY_groups[i];
  }
  
  return true;
}

bool partial_geometry::get_downstream_XY_groups(int& a_number_of_downstream_XY_groups, int a_downstream_XY_groups[])
// this also includes triple groups
{
  a_number_of_downstream_XY_groups = number_of_downstream_XY_groups;
  
  for (int i = 0; i < number_of_downstream_XY_groups; i++)
  {
    a_downstream_XY_groups[i] = downstream_XY_groups[i];
  }
  
  return true;
}

bool partial_geometry::get_types_of_midstream_plates(int& number_of_x_plates, ssd_plate x_plates[], int x_plate_indices[], int& number_of_y_plates, ssd_plate y_plates[], int y_plate_indices[], int& number_of_d1_plates, ssd_plate d1_plates[], int d1_plate_indices[], int& number_of_d2_plates, ssd_plate d2_plates[], int d2_plate_indices[], int& number_of_d_plates, ssd_plate d_plates[], int d_plate_indices[])
{
  number_of_d1_plates = 0;
  number_of_d2_plates = 0;
  
  number_of_x_plates = number_of_midstream_X_plates;
  for (int i = 0; i < number_of_midstream_X_plates; i++)
  {
    x_plates[i] = the_plates[midstream_X_plates[i]];
    x_plate_indices[i] = midstream_X_plates[i];
  }
  
  number_of_y_plates = number_of_midstream_Y_plates;
  for (int i = 0; i < number_of_midstream_Y_plates; i++)
  {
    y_plates[i] = the_plates[midstream_Y_plates[i]];
    y_plate_indices[i] = midstream_Y_plates[i];
  }
  
  number_of_d1_plates = number_of_d1_plates_midstream;
  for (int i = 0; i < number_of_d1_plates_midstream; i++)
  {
    d1_plates[i] = the_plates[d1_plates_midstream[i]];
    d1_plate_indices[i] = d1_plates_midstream[i];
  }
  
  number_of_d2_plates = number_of_d2_plates_midstream;
  for (int i = 0; i < number_of_d2_plates_midstream; i++)
  {
    d2_plates[i] = the_plates[d2_plates_midstream[i]];
    d2_plate_indices[i] = d2_plates_midstream[i];
  }
  
  number_of_d_plates = number_of_diagonal_plates_midstream;
  for (int i = 0; i < number_of_diagonal_plates_midstream; i++)
  {
    d_plates[i] = the_plates[diagonal_plates_midstream[i]];
    d_plate_indices[i] = diagonal_plates_midstream[i];
  }
  
  return true;
}

// rename this to find_line_plane_intersection_point
bool find_line_plate_intersection_point(line_3d a_track_line, plane_3d plate_plane, XYZVector& intersection_point)
{
  double t;
  double dummy;
  
  dummy = plate_plane.a * a_track_line.direction.x() + plate_plane.b * a_track_line.direction.y() + plate_plane.c * a_track_line.direction.z();
  
  if (dummy == 0) return false;
  
  t = - (plate_plane.d + plate_plane.a * a_track_line.point.x() + plate_plane.b * a_track_line.point.y() + plate_plane.c * a_track_line.point.z()) / dummy;
  
  intersection_point = a_track_line.point + t * a_track_line.direction;
  
  return true;
}

double find_point_to_line_3d_distance(XYZVector a_point, line_3d a_line)
{
  return sqrt(a_line.direction.Cross(a_point - a_line.point).Mag2());
}

bool find_line_to_line_distance(line_3d line_1, line_3d line_2, double& distance)
{
  if (line_1.direction.Cross(line_2.direction).Mag2() == 0)
  {
    // arbitrarty small number to mean zero
    if (fabs(sqrt( (line_1.point - line_2.point).Cross(line_1.direction).Mag2())) < 0.000000001)
    {
      distance = 0;
      return true;
    }
    else
    {
      return false;
    }
  }
  
  distance = fabs((line_1.point - line_2.point).Dot(line_1.direction.Cross(line_2.direction))) / sqrt(line_1.direction.Cross(line_2.direction).Mag2());
  
  return true;
}

bool find_closest_point(line_3d line_1, line_3d line_2, XYZVector& closest_point)
{
  XYZVector d;
  double t1, t2;
  
  d = line_1.point - line_2.point;
  
  t1 = (line_1.direction.Cross(line_2.direction)).Dot(d.Cross(line_2.direction)) / sqrt((line_1.direction.Cross(line_2.direction)).Mag2());
  
  t2 = (line_1.direction.Cross(line_2.direction)).Dot(d.Cross(line_1.direction)) / sqrt((line_1.direction.Cross(line_2.direction)).Mag2());
  
  closest_point = 0.5 * (line_1.point + t1 * line_1.direction + line_2.point + t2 * line_2.direction);
  
  return true;
}

void track::sort_cluster_lines_by_plate_index()
{
  pair<line_3d, int> a_pair;
  vector<pair<line_3d, int>> the_pairs;
  
  for (int i = 0; i < cluster_lines.size(); i++)
  {
    a_pair.first = cluster_lines.at(i);
    a_pair.second = plate_index_of_cluster_line.at(i);
    the_pairs.push_back(a_pair);
  }
    
  sort(the_pairs.begin(), the_pairs.end(), compare_cluster_line_plate);
  
  cluster_lines.clear();
  plate_index_of_cluster_line.clear();
  
  for (int i = 0; i < the_pairs.size(); i++)
  {
    cluster_lines.push_back(the_pairs.at(i).first);
    plate_index_of_cluster_line.push_back(the_pairs.at(i).second);
  }
}

bool track::find_intersection_points_no_field()
{
  // used for no field and in field for upstream and maybe midstream regions; uses the assumed start and end groups
  // assumes in the regions to be used or previous regions there are no multi-plaque plates, so plaque and plate numbers are the same
  
  int first_plate = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(assumed_start_group)->plate_index_list[0];
  int last_plate = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(assumed_end_group)->plate_index_list[the_global_data_dispatcher->get_partial_geometry()->get_plate_group(assumed_end_group)->number_of_plates_in_group - 1];
  
  for (int i = first_plate; i <= last_plate; i++)
  {
    track_plate_intersection_points[i].SetCoordinates(0., 0., 0.);
  }
  
  for (int i = first_plate; i <= last_plate; i++)
  {
    find_line_plate_intersection_point(track_line, the_global_data_dispatcher->get_partial_geometry()->get_plaque(the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->member_plaques[0])->plaque_plane_equation, track_plate_intersection_points[i]); 
  }
  
  return true;
}

bool track::find_intersection_points_rough()
{
  for (int i = 0; i < max_number_of_plates; i++)
  {
    track_plate_intersection_points[i].SetCoordinates(0., 0., 0.);
  }
  
  for (int i = 0; i < max_number_of_plates; i++)
  {
    if (i >= plate_index_of_cluster_line.at(0) && i <= plate_index_of_cluster_line.at(plate_index_of_cluster_line.size() - 1)) // assumes cluster lines are oredered by plate index
    {
      find_line_plate_intersection_point(track_line, the_global_data_dispatcher->get_partial_geometry()->get_plaque(the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->member_plaques[0])->plaque_plane_equation, track_plate_intersection_points[i]); // this is rough; the plate plane is assumed to be the same as the first plaque plane
    }
  }
 
  return true;
}

void pattern::reset()
{
  number_of_x_lines = 0;
  number_of_y_lines = 0;
  equations_x.clear();
  equations_y.clear();
  pattern_x_coordinates.clear();
  pattern_y_coordinates.clear();
  x_tracks.clear();
  y_tracks.clear();
  combined_tracks.clear();
}

bool find_plane_of_two_parallel_lines(line_3d line_1, line_3d line_2, plane_3d& a_plane)
{
  bool plane_found = true;
  
  XYZVector the_normal;
  
  if (line_1.direction.Dot(line_2.direction) > 0.9999) // redo with cut
  {
    the_normal = line_1.direction.Cross(line_1.point - line_2.point);
    the_normal = the_normal.Unit();
    a_plane.a = the_normal.X();
    a_plane.b = the_normal.Y();
    a_plane.c = the_normal.Z();
    a_plane.d = - the_normal.Dot(line_1.point); // redo with the point that has the lower uncertainty
  }
  else
  {
    plane_found = false;
  }
  
  return plane_found;
}

void find_3d_line_from_projections(double x[], double z_x[], double y[], double z_y[], XYZVector& point, XYZVector& direction)
{
  XYZVector t_x, t_y;
  
  t_x.SetCoordinates(x[1] - x[0], 0., z_x[1] - z_x[0]);
  t_x = t_x / (t_x.Z());
  
  t_y.SetCoordinates(0., y[1] - y[0], z_y[1] - z_y[0]);
  t_y = t_y / (t_y.Z());
  
  direction.SetCoordinates(t_x.X(), t_y.Y(), t_x.Z());
  
  point.SetCoordinates(x[0] - z_x[0] * (x[1] - x[0]) / (z_x[1] - z_x[0]), y[0] - z_y[0] * (y[1] - y[0]) / (z_y[1] - z_y[0]), 0.);
}

XYZVector find_point_on_line_by_z(line_3d a_line, double a_z)
{
  XYZVector the_point;
  
  the_point = a_line.point + a_line.direction * ((a_z - a_line.point.Z()) / a_line.direction.Z());
  
  return the_point;
}

void combine_3d_line_from_projections(track x_track, track y_track, line_3d& a_line)
{
  XYZVector t_x, t_y;
  
  double z_middle, z_first, z_last;
  
  if (x_track.cluster_lines[0].point.Z() < y_track.cluster_lines[0].point.Z())
  {
    z_first = x_track.cluster_lines[0].point.Z();
  }
  else
  {
    z_first = y_track.cluster_lines[0].point.Z();
  }
  
  if (x_track.cluster_lines[1].point.Z() < y_track.cluster_lines[1].point.Z())
  {
    z_last = x_track.cluster_lines[1].point.Z();
  }
  else
  {
    z_last = y_track.cluster_lines[1].point.Z();
  }
  
  z_middle = (z_first + z_last) / 2.;
  
  t_x.SetCoordinates(x_track.cluster_lines[1].point.X() - x_track.cluster_lines[0].point.X(), 0., x_track.cluster_lines[1].point.Z() - x_track.cluster_lines[0].point.Z());
  t_x = t_x / (t_x.Z());
  
  t_y.SetCoordinates(0., y_track.cluster_lines[1].point.Y() - y_track.cluster_lines[0].point.Y(), y_track.cluster_lines[1].point.Z() - y_track.cluster_lines[0].point.Z());
  t_y = t_y / (t_y.Z());
  
  a_line.direction.SetCoordinates(t_x.X(), t_y.Y(), t_x.Z());
  a_line.direction = a_line.direction.Unit();
    
  a_line.point.SetCoordinates(x_track.cluster_lines[0].point.X() - x_track.cluster_lines[0].point.Z() * (x_track.cluster_lines[1].point.X() - x_track.cluster_lines[0].point.X()) / (x_track.cluster_lines[1].point.Z() - x_track.cluster_lines[0].point.Z()), y_track.cluster_lines[0].point.Y() - y_track.cluster_lines[0].point.Z() * (y_track.cluster_lines[1].point.Y() - y_track.cluster_lines[0].point.Y()) / (y_track.cluster_lines[1].point.Z() - y_track.cluster_lines[0].point.Z()), 0.);
  
  a_line.point.SetX(a_line.point.X() + z_middle * a_line.direction.X() / a_line.direction.Z());
  a_line.point.SetY(a_line.point.Y() + z_middle * a_line.direction.Y() / a_line.direction.Z());
  a_line.point.SetZ(z_middle);
}

bool compare_multi_track_selections(const multi_track_selections& first, const multi_track_selections& second)
{
  return (first.selection_quality < second.selection_quality);
}

bool is_intersection_point_inside_plaque(XYZVector intersection_point, ssd_plaque* p_dummy_plaque)
{
  bool operation_success;
  
  RotationZ rz;
  XYZVector dummy_point;
  
  // assuming reasonably good alignment; can be done with the actual plane
  rz.SetAngle(p_dummy_plaque->rotation.z() * to_radians);
  dummy_point = rz * intersection_point;
  
  if (dummy_point.X() > p_dummy_plaque->position.X() - p_dummy_plaque->size.Y() / 2. &&
      dummy_point.X() < p_dummy_plaque->position.X() + p_dummy_plaque->size.Y() / 2. &&
      dummy_point.Y() > p_dummy_plaque->position.Y() - p_dummy_plaque->size.X() / 2. &&
      dummy_point.Y() < p_dummy_plaque->position.Y() + p_dummy_plaque->size.X() / 2.)
  {
    operation_success = true;
  }
  else
  {
    operation_success = false;
  }
  
  return operation_success;
}

void find_2d_plane_equation(string subspace, XYZVector point_1, XYZVector point_2, line_3d the_equation_line)
{
  // !! needs implementation
}
