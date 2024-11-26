# PIEZO cluster track analysis

This set of Matlab scripts are used to analyze output tables produced by the ImageJ plugins Trackmate and TraJClassifier, used initially on live TIRF microscopy imaging of PIEZO-tag clusters movement (see Schaefer et al., JBC 2023: https://doi.org/10.1016/j.jbc.2023.104782 , or Verkest et al., Nat.Comm 2022 https://doi.org/10.1038/s41467-022-28974-6 for examples). 

It allows visualization of the different classes of track identified by the two plugins (static image with tracks overlay, as well as overlay video of the tracks and the original microscopy movie) ,  produces MSD graphs and summary tables about various parameters (proportion for each subtype of track, duration, intensity, diameter etc.). 

![](https://github.com/clementverkest/PIEZO_Track_Analysis/tracks_movie.gif)

# Requirements
Scripts generated and tested with Matlab (v2023b). Required input data tables produced with Trackmate (v7.13.2) and TraJClassifier.

 Required main scripts (TIRF_Tracks_MSD.m) and associated subroutines (ExtractAndAssignTracks.m, PlotAllTracks.m, calculate_msd.m, average_msd.m, calculate_diffusion_coefficient.m, GaussFit2D.m, ExtractMSD_AVG_Results.m, ExtractTracksClustersResults.m)



# How to use
- Trackmate data tables required and needed to be exported from the ImageJ plugin: allspots .csv table, spots and tracks (optional) .csv tables, Track in xml format to be read by TraJclassifier
- TraJClassifier data tables: all four track classes and Parents .csv tables
- All the tables should be stored in a single folder, together with the corresponding .tif stack used in Trackmate and TraJClassifier
- Make sure all the subroutines required are accessible by Matlab
- Select the desire options at the beginning of the script TIRF_Tracks_MSD.m, by typing either "true" or "false". Make sure to indicate the proper pixel size before running TIRF_Tracks_MSD.m (should be the same value as the one used during Trackmate analysis)
- Run TIRF_Tracks_MSD.m, select the folder where all the tables and tif stack are
- Generates csv result tables in a new folder, as well as graphs.
- An example is provided as a zip file: download and unzip it. It contains the required trackmate and trajclassifier tables, as well as the tif file of the first frame of the TIRF acquisition sequence.




# Citation
Upcoming work Verkest et al., 2024, in preparation.

Do not forget to cite Trackmate v7 related paper (DOI: 10.1038/s41592-022-01507-1) and TraJClassifier (DOI: 10.1371/journal.pone.0170165)
