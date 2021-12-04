layered_Vs_1D_CWB=[1.96 2.62 3.03 3.35 3.61 3.71 3.95 4.21 4.49 4.68 4.72 4.79 4.8 4.74 4.86 4.92 5.49];
layered_Vp_1D_CWB=[3.48 4.48 5.25 5.83 6.21 6.41 6.83 7.29 7.77 8.05 8.16 8.34 8.35 8.2 8.4 8.51 8.7];
layers=length(layered_Vs_1D_CWB); %how many layers in 1D velocity structure
bld3=0.1; %degree of grid block size in Horizontal?
bld4=1; %km of grid block size in vertical?
lon_c=118:bld3:125;lat_c=19.5:bld3:26.5;
layered_D_lower=[2 4 9 13 17 25 30 35 50 70 90 110 140 170 200 240 280];
layered_D_upper(2:length(layered_D_lower))=layered_D_lower(1:length(layered_D_lower)-1);layered_D_upper(1)=0;
dep_c=0:1:279;
for kk=1:length(dep_c);
  index_start_tmp=find((layered_D_lower-dep_c(kk))>0);tmp1=find(min(layered_D_lower(index_start_tmp)-dep_c(kk)));index_start(kk)=index_start_tmp(tmp1);
  if(index_start(kk)>1&layered_D_lower(index_start(kk)-1)==dep_c(kk))
    index_start(kk)=index_start(kk)-1;
  end
end

dep_c(length(dep_c)+1)=dep_c(length(dep_c))+bld4;
for kk=1:length(dep_c);
  if(kk<length(dep_c))
    vel_p(1:length(lon_c),1:length(lat_c),kk)=layered_Vp_1D_CWB(index_start(kk));
  elseif(kk==length(dep_c))
    vel_p(1:length(lon_c),1:length(lat_c),kk)=8.7;
  end
end
for kk=1:length(dep_c);
  if(kk<length(dep_c))
    vel_s(1:length(lon_c),1:length(lat_c),kk)=layered_Vs_1D_CWB(index_start(kk));
  elseif(kk==length(dep_c))
    vel_s(1:length(lon_c),1:length(lat_c),kk)=5.49;
  end
end
nlon_c=length(lon_c);nlat_c=length(lat_c);ndep_c=length(dep_c);%ndep_c=length(layered_Vs_1D_CWB);
fid=fopen('MOD_CWB1D','w');
fprintf(fid,'%3.1f %g %g %g %g\n',bld3,bld4,nlon_c,nlat_c,ndep_c);
fprintf(fid,'%g ',lon_c(1:nlon_c));fprintf(fid,'\n');
fprintf(fid,'%g ',lat_c(1:nlat_c));fprintf(fid,'\n');
fprintf(fid,'%g ',dep_c(1:ndep_c));fprintf(fid,'\n');
for k=1:ndep_c
  for j=1:nlat_c
    fprintf(fid,'%4.2f ',vel_p(1:nlon_c,j,k));fprintf(fid,'\n');
  end
end
for k=1:ndep_c
  for j=1:nlat_c
    fprintf(fid,'%4.2f ',vel_s(1:nlon_c,j,k));fprintf(fid,'\n');
  end
end

