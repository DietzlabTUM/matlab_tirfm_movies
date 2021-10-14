# matlab_tirfm_movies

These and other related scripts were written by Dr. Matthias Schickinger, Dr. Philip Ketterer, Dr. Jonas Funke, Anna-Katharina Pumm, Dr. Wouter Engelen and Eva Bertosin

The main scripts are in matlab_tirfm_movies folder, named from 1 to 5.

1-process_movies_manual:
load movies, compute STD DEV of whole movie, pick and export all particles (spots)
To change: run('/Users/username/Documents/MOVIE_TIRFM/functions/startup.m')

2-get_std_dev_of_spots:
calculate STD DEV of all spots
to change: folders to load spots

3-classify_sum_tifs:
classification of all spots according to their STD DEV (classes: stationary, diffusive, moving, unclear)
only "moving" spots will be further processed

4-analyse_sorted_particles:
find center of spots, fit PSF to spots, find spots in each frame with VWCM
to change: folders to load spots

5-classify_spots:
classify moving spots into switching between 2 or 3 positions, stationary, rotating, diffusive (no full circle) and unclear
calculare RMSD, cumulative angle, angular traces, angular velocity, distance from center for desired type of spots (switching, rotating, diffusive)

The folders "functions" and "TOOLBOX_MOVIES" contain support scripts.

DOI: 10.5281/zenodo.5568757
