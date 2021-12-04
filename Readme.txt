Step 1.
  calculated ray path and takeoff angle as well from ray tracing method in 1D layered Vs
  1. use [make_1DMOD.m] to construct velocity grid model (ex. MOD_CWB1D)

  2. use fortran code [pbr] to calculate ray path for S wave (code is modified from Tsai, 2015)
     format of input is (pbr  EventLongitude EventLatitude EventDepth StaLongitude StaLatitude StaDepth), StaDepth det on -200 meter 
       depth to prevent out of velocity model
     output is 'out_route_P' and 'out_route_S' for P and S wave (they are almost the same if the velocity structure is consistent
     output format including 6 lines as line 1:n_route (how many route points) ,s_treval_time (in total), 2:source [Lat] [Lon] [Dep], 
       3:sta [Lat] [Lon] [Dep], 4:dep_from_source (n_route points) to station, 5:lat(n_route points), 6:lon(n_route points)

  3. use matlab function [cal_route_takeoff_angle.m] to get takeoff angle
     format of input is EventLongitude EventLatitude EventDepth StaLongitude StaLatitude
       [takeoff_angle_P,treval_time_P,takeoff_angle_S,treval_time_S]=cal_route_takeoff_angle(evlo,evla,evdp,stlo,stla)
     
Step 2. 
  calculated radiation pattern using fault informations such as (Strike,Dip,Rake,Takeoff Angle,AZIMUTH from source to station)
  1. use matlab function [rpgen.m] to get radiation pattern (code is from reference Grzegorz Kwiatek, 2020)
     Due to the code traditionally is for calculate radiation pattern of tensile rupture on the fault, the gamma should be set as 0 to get radiation pattern for pure shear rupture

  2. the output GS is the RMS of two Horizontal, GSH and GSV are the radiation pattern for individual directions
     [get_radiation.m] is the code I used to combine SSHAC flatfile and calculated radiation pattern
     take (tmp1_database_TW_ShaCru.mat) as input, it is previous part of FAS flatfile in SSHAC database
     the output result is (database_TW_ShaCru.csv) for example


References
Tsai, C.Y. (2015), The rapid travel-time sequence method for earthquake locating and application of the local earthquake early warning arrays in Taiwan, master thesis of National Cheng-Kung University, pp182. (in Chinese with English abstract)
Grzegorz Kwiatek (2020). Radiation pattern from shear-tensile seismic source (https://www.mathworks.com/matlabcentral/fileexchange/43524-radiation-pattern-from-shear-tensile-seismic-source), MATLAB Central File Exchange. Retrieved August 7, 2020.