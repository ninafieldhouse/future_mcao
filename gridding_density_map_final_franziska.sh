#!/bin/bash

module load dyn_tools

cd /net/thermo/atmosdyn/fschnyde/data/lagranto/cao_traj_density_final

# file_events="/net/thermo/atmosdyn/fschnyde/data/lagranto/case_study/all_times_993727"

# nlines=$(wc -l $file_events | awk '{print $1}')
# echo $nlines

# for ((i=1; i<=$nlines; i++)); do
# 	# read date from datelist
# 	date=$(awk -v line=$i 'NR==line{print $1}' $file_events)
#     echo $date

#     gridding -npa -a gridding_input_BCAR_warming_along_trajectory_DTHweighted_case_study_${date} density_map_BCAR_warming_along_trajectory_DTHweighted_case_study_${date} 100 -180 180 -90 90 0.5 0.5
    
# done

# NA
# gridding -npa gridding_input_BCAR_trajectory_origin_NA_12h_warming_rate density_map_BCAR_trajectory_origin_NA_12h_warming_rate_nearest_neighbour 100 -180 180 -90 90 0.5 0.5
# gridding -npa -a gridding_input_BCAR_warming_along_trajectory_DTHweighted_NA_12h_warming_rate density_map_BCAR_warming_along_trajectory_DTHweighted_NA_12h_warming_rate_nearest_neighbour 0 -180 180 -90 90 0.5 0.5

# gridding -npa gridding_input_BCAR_trajectory_origin_NA_12h_cooling_rate density_map_BCAR_trajectory_origin_NA_12h_cooling_rate_nearest_neighbour 100 -180 180 -90 90 0.5 0.5
# gridding -npa -a gridding_input_BCAR_cooling_along_trajectory_DTHweighted_NA_12h_cooling_rate density_map_BCAR_cooling_along_trajectory_DTHweighted_NA_12h_cooling_rate_nearest_neighbour 0 -180 180 -90 90 0.5 0.5

gridding -npa gridding_input_BCAR_trajectory_origin_NA_12h_cooling_rate_group2_2days density_map_BCAR_trajectory_origin_NA_12h_cooling_rate_group2_2days_nearest_neighbor 100 -180 180 -90 90 0.5 0.5
gridding -npa gridding_input_BCAR_trajectory_origin_NA_12h_cooling_rate_group2_4days density_map_BCAR_trajectory_origin_NA_12h_cooling_rate_group2_4days_nearest_neighbor 100 -180 180 -90 90 0.5 0.5
gridding -npa gridding_input_BCAR_trajectory_origin_NA_12h_cooling_rate_group2_6days density_map_BCAR_trajectory_origin_NA_12h_cooling_rate_group2_6days_nearest_neighbor 100 -180 180 -90 90 0.5 0.5

# Siberia
# gridding -npa gridding_input_BCAR_trajectory_origin_Siberia_12h_warming_rate density_map_BCAR_trajectory_origin_Siberia_12h_warming_rate_nearest_neighbour 100 -180 180 -90 90 0.5 0.5
# gridding -npa -a gridding_input_BCAR_warming_along_trajectory_DTHweighted_Siberia_12h_warming_rate density_map_BCAR_warming_along_trajectory_DTHweighted_Siberia_12h_warming_rate_nearest_neighbour 0 -180 180 -90 90 0.5 0.5

# gridding -npa gridding_input_BCAR_trajectory_origin_Siberia_12h_cooling_rate density_map_BCAR_trajectory_origin_Siberia_12h_cooling_rate_nearest_neighbour 100 -180 180 -90 90 0.5 0.5
# gridding -npa -a gridding_input_BCAR_cooling_along_trajectory_DTHweighted_Siberia_12h_cooling_rate density_map_BCAR_cooling_along_trajectory_DTHweighted_Siberia_12h_cooling_rate_nearest_neighbour 0 -180 180 -90 90 0.5 0.5

# gridding -npa gridding_input_BCAR_trajectory_origin_Siberia_12h_cooling_rate_group2_2days density_map_BCAR_trajectory_origin_Siberia_12h_cooling_rate_group2_2days_nearest_neighbor 100 -180 180 -90 90 0.5 0.5
# gridding -npa gridding_input_BCAR_trajectory_origin_Siberia_12h_cooling_rate_group2_4days density_map_BCAR_trajectory_origin_Siberia_12h_cooling_rate_group2_4days_nearest_neighbor 100 -180 180 -90 90 0.5 0.5
# gridding -npa gridding_input_BCAR_trajectory_origin_Siberia_12h_cooling_rate_group2_6days density_map_BCAR_trajectory_origin_Siberia_12h_cooling_rate_group2_6days_nearest_neighbor 100 -180 180 -90 90 0.5 0.5

# gridding -npa gridding_input_BCAR_last_time_step_NA density_map_BCAR_last_time_step_NA 100 -180 180 -90 90 0.5 0.5
gridding -npa gridding_input_BCAR_last_time_step_Siberia density_map_BCAR_last_time_step_Siberia 100 -180 180 -90 90 0.5 0.5


