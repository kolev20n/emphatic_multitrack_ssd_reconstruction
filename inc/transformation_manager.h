// Class for ssd to lab reference frame transformation.
// Includes rotations and translations due to misalignment.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef transformation_manager_h
#define transformation_manager_h

#include <vector>

#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>

#include "common_data_structures_and_functions.h"
#include "partial_geometry.h"

class input_output_manager;
class partial_geometry;

class transformation_manager
{
public:
  transformation_manager(partial_geometry* a_geometry);
  bool transform_to_global_tracking_frame(std::vector<double> clusters[], std::vector<line_3d> cluster_lines[]);
  bool merge_plaques_in_same_plate(std::vector<line_3d> plaque_cluster_lines[], std::vector<line_3d> cluster_lines[]);
private:
  partial_geometry* the_geometry;
};

#endif
