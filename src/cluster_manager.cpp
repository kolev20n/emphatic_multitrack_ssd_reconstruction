#include <vector>
#include <iostream>
#include <string>

#include "common_data_structures_and_functions.h"
#include "cluster_manager.h"


using namespace std;

cluster_manager::cluster_manager()
{
  cout << "Created the_cluster_manager." << endl;
}

bool cluster_manager::cluster(vector<hit_pair> hit_pair_list[], vector<double> clusters[], int actual_number_of_plaques, string cluster_interpretation_method)
{
  int cluster_multiplicity[max_number_of_plaques][max_number_of_strips_on_plaque];
  // theoretically there can be as many clusters in the plane as hits
  int clusters_per_plane[max_number_of_plaques][max_number_of_strips_on_plaque];
  bool single_cluster_in_each_plane[max_number_of_plaques];
  vector<hit_pair> dummy_hit_pair_vector;
  vector<vector<hit_pair>> temp_clusters[max_number_of_plaques];

  int index1, index2;
  int cluster_multiplicity_counter;
  int hit_multiplicity_counter;
  
  int hit_counter = 0;
  double sum_energy_value;
  double sum_strip_number;
  
  for (int i = 0; i < max_number_of_plaques; i++)
  {
    single_cluster_in_each_plane[i] = true;
    
    for (int j = 0; j < max_number_of_strips_on_plaque; j++)
    {
      cluster_multiplicity[i][j] = 0;
      clusters_per_plane[i][j] = 0;
    }
  }
  
  for (int i = 0; i < actual_number_of_plaques; i++)
  {
    // hit pairs seems to come sorted by strip number from the simulation,
    // this is not necessarily true from experiment, is it?
    // so just in case
    sort(hit_pair_list[i].begin(), hit_pair_list[i].end(), compare_hit_pairs);
    
    temp_clusters[i].clear();
    dummy_hit_pair_vector.clear();
    
    if (hit_pair_list[i].size() >= max_number_of_strips_on_plaque)
    {
      cout << "More than maximum assumed hits on a plaque: " << hit_pair_list[i].size() << " is more than " << max_number_of_strips_on_plaque << endl;
      cout << "Adjust value in cluster_manager.cpp!" << endl;
      cout << "Quitting..." << endl;
      exit(EXIT_FAILURE);
    }
    
    if (hit_pair_list[i].empty())
    {
      cluster_multiplicity[i][0]++;
      clusters_per_plane[i][0]++;
      single_cluster_in_each_plane[i] = false;
    }
    else
    {
      if (hit_pair_list[i].size() == 1)
      {
        cluster_multiplicity[i][1]++;
        clusters_per_plane[i][1]++;
        dummy_hit_pair_vector.push_back(hit_pair_list[i].at(0));
        temp_clusters[i].push_back(dummy_hit_pair_vector);
      }
      else
      {
        index1 = 0;
        index2 = 0;
        
        cluster_multiplicity_counter = 0;
        
        // dummy_hit_pair_vector.push_back(hit_pair_list[j].at(index1));
        // clusters[j].push_back(dummy_hit_pair_vector);
        
        dummy_hit_pair_vector.push_back(hit_pair_list[i].at(index1));
        
        do
        {
          hit_multiplicity_counter = 0;
          
          do
          {
            hit_multiplicity_counter++;
            index2++;
            
            if (index2 == hit_pair_list[i].size())
            {
              if (dummy_hit_pair_vector.size() != 0) temp_clusters[i].push_back(dummy_hit_pair_vector);
              
              break;
            }
            
            if (hit_pair_list[i].at(index2).strip_number == hit_pair_list[i].at(index2 - 1).strip_number + 1)
            {
              dummy_hit_pair_vector.push_back(hit_pair_list[i].at(index2));
            }
            else
            {
              temp_clusters[i].push_back(dummy_hit_pair_vector);
              dummy_hit_pair_vector.clear();
              dummy_hit_pair_vector.push_back(hit_pair_list[i].at(index2));
            }
            
            
          }
          while (hit_pair_list[i].at(index1).strip_number + index2 - index1 == hit_pair_list[i].at(index2).strip_number);
          
          cluster_multiplicity[i][hit_multiplicity_counter]++;
          cluster_multiplicity_counter++;
          
          index1 = index2;
        }
        while (index1 < hit_pair_list[i].size());
        
        clusters_per_plane[i][cluster_multiplicity_counter]++;
        
        if (cluster_multiplicity_counter != 1)
        {
          single_cluster_in_each_plane[i] = false;
        }
      }
    }
    
    /*
    for (int j = 0; j < temp_clusters[i].size(); j++)
    {
      cout << "(";
      for (int k = 0; k < temp_clusters[i].at(j).size(); k++)
      {
        cout << "(" << temp_clusters[i].at(j).at(k).strip_number << " " << temp_clusters[i].at(j).at(k).strip_number << ") ";
      }
      cout << "\b) ";
    }
    cout << endl;
    */
     
    if (cluster_interpretation_method == "weighted")
    {
      for (int j = 0; j < temp_clusters[i].size(); j++)
      {
        sum_energy_value = 0.;
        sum_strip_number = 0.;
          
        for (int k = 0; k < temp_clusters[i].at(j).size(); k++)
        {
          hit_counter++;
          sum_energy_value += temp_clusters[i].at(j).at(k).energy_value;
          sum_strip_number += temp_clusters[i].at(j).at(k).energy_value * temp_clusters[i].at(j).at(k).strip_number;
        }

        clusters[i].push_back(sum_strip_number / sum_energy_value);
      }
    }
    
  } // end loop on plates
  
  return true;
  
}

