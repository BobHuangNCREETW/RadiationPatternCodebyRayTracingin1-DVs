function [takeoff_angle_P,treval_time_P,takeoff_angle_S,treval_time_S]=cal_route_takeoff_angle(evlo,evla,evdp,stlo,stla)
unix(['pbr ',num2str(evlo),' ',num2str(evla),' ',num2str(evdp),' ',num2str(stlo),' ',num2str(stla),' -200']); %elevation of station fixed to -200m to prevent out of range
fid=fopen('out_route_P','r');
nu=fgetl(fid);
route_point_num_tmp=fgetl(fid);
[route_point_num_P,treval_time_P]=strread(route_point_num_tmp,'%f %f');
source_tmp=fgetl(fid);[evla,evlo,evdp]=strread(source_tmp,'%f %f %f');
sta_tmp=fgetl(fid);[stla,stlo,stdp]=strread(sta_tmp,'%f %f %f');
[ECEF_source_X,ECEF_source_Y,ECEF_source_Z]=geog2ECEF(evla,evlo,-1000*evdp); %depth should be nagetive
[ECEF_sta_X,ECEF_sta_Y,ECEF_sta_Z]=geog2ECEF(stla,stlo,-1000*stdp);
[ECEF_core_X,ECEF_core_Y,ECEF_core_Z]=geog2ECEF(0,0,-6378137); %earth core point
cordi_route=fscanf(fid,'%f',[route_point_num_P 3]);
route_dep_P=cordi_route(:,1);route_lat_P=cordi_route(:,2);route_lon_P=cordi_route(:,3);
point_num=2;
x_start=route_lon_P(1);
y_start=route_lat_P(1);
z_start=route_dep_P(1);
x_start1=route_lon_P(point_num);
y_start1=route_lat_P(point_num);
z_start1=route_dep_P(point_num);
[X0,Y0,Z0]=geog2ECEF(0,0,-6378137); %earth core point
[ECEF_route_P1_X,ECEF_route_P1_Y,ECEF_route_P1_Z]=geog2ECEF(y_start,x_start,-1000*z_start);
[ECEF_route_P2_X,ECEF_route_P2_Y,ECEF_route_P2_Z]=geog2ECEF(y_start1,x_start1,-1000*z_start1);
[proj_route_P1_X,proj_route_P1_Y,proj_route_P1_Z]=proj_3Dto2D(ECEF_source_X,ECEF_source_Y,ECEF_source_Z,ECEF_sta_X,ECEF_sta_Y,ECEF_sta_Z,ECEF_core_X,ECEF_core_Y,ECEF_core_Z,ECEF_route_P1_X,ECEF_route_P1_Y,ECEF_route_P1_Z);
[proj_route_P2_X,proj_route_P2_Y,proj_route_P2_Z]=proj_3Dto2D(ECEF_source_X,ECEF_source_Y,ECEF_source_Z,ECEF_sta_X,ECEF_sta_Y,ECEF_sta_Z,ECEF_core_X,ECEF_core_Y,ECEF_core_Z,ECEF_route_P2_X,ECEF_route_P2_Y,ECEF_route_P2_Z);
P1=[proj_route_P1_X,proj_route_P1_Y,proj_route_P1_Z];
P2=[proj_route_P2_X,proj_route_P2_Y,proj_route_P2_Z];
P0=[X0,Y0,Z0];
[rad,ang]=pvp_sub(P1,P0,P2);
takeoff_angle_P=ang;
fclose(fid);
%%
fid=fopen('out_route_S','r');
nu=fgetl(fid);
route_point_num_tmp=fgetl(fid);
[route_point_num_S,treval_time_S]=strread(route_point_num_tmp,'%f %f');
source_tmp=fgetl(fid);[evla,evlo,evdp]=strread(source_tmp,'%f %f %f');
sta_tmp=fgetl(fid);[stla,stlo,stdp]=strread(sta_tmp,'%f %f %f');
cordi_route=fscanf(fid,'%f',[route_point_num_S 3]);
route_dep_S=cordi_route(:,1);route_lat_S=cordi_route(:,2);route_lon_S=cordi_route(:,3);
x_start=route_lon_S(1);y_start=route_lat_S(1);z_start=route_dep_S(1);
x_start1=route_lon_S(2);y_start1=route_lat_S(2);z_start1=route_dep_S(2);
[X0,Y0,Z0]=geog2ECEF(0,0,-6378137); %earth core point
[ECEF_route_S1_X,ECEF_route_S1_Y,ECEF_route_S1_Z]=geog2ECEF(y_start,x_start,-1000*z_start);
[ECEF_route_S2_X,ECEF_route_S2_Y,ECEF_route_S2_Z]=geog2ECEF(y_start1,x_start1,-1000*z_start1);
[proj_route_S1_X,proj_route_S1_Y,proj_route_S1_Z]=proj_3Dto2D(ECEF_source_X,ECEF_source_Y,ECEF_source_Z,ECEF_sta_X,ECEF_sta_Y,ECEF_sta_Z,ECEF_core_X,ECEF_core_Y,ECEF_core_Z,ECEF_route_S1_X,ECEF_route_S1_Y,ECEF_route_S1_Z);
[proj_route_S2_X,proj_route_S2_Y,proj_route_S2_Z]=proj_3Dto2D(ECEF_source_X,ECEF_source_Y,ECEF_source_Z,ECEF_sta_X,ECEF_sta_Y,ECEF_sta_Z,ECEF_core_X,ECEF_core_Y,ECEF_core_Z,ECEF_route_S2_X,ECEF_route_S2_Y,ECEF_route_S2_Z);
P1=[proj_route_S1_X,proj_route_S1_Y,proj_route_S1_Z];
P2=[proj_route_S2_X,proj_route_S2_Y,proj_route_S2_Z];
P0=[X0,Y0,Z0];
[rad,ang]=pvp_sub(P1,P0,P2);
takeoff_angle_S=ang;
fclose(fid);

function [rad,ang]=pvp_sub(p,pA,pB)
vA=(pA-p)/norm(pA-p);
vB=(pB-p)/norm(pB-p);
[rad,ang]=VvV(vA,vB);

function [rad,ang]=VvV(vA,vB)
vA=vA/norm(vA);
vB=vB/norm(vB);
rad=acos(dot(vA,vB));
ang=rad*180/pi;

function [X,Y,Z]=geog2ECEF(lat,lon,height)
%height in meter, surface=0, core=-6378137;
ellipsoid(1)=6378137; %長軸
ellipsoid(2)=298.257223563; %扁率
phi=D2R(lat);
small_lambda=D2R(lon);
b=ellipsoid(1)-(ellipsoid(1)*(1/ellipsoid(2)));
e2=(ellipsoid(1)^2-b^2)/(ellipsoid(1)^2);
v=ellipsoid(1)/((1-e2*(sin(phi)^2))^0.5);
X=double((v+height)*cos(phi)*cos(small_lambda));
Y=double((v+height)*cos(phi)*sin(small_lambda));
Z=double((v*(1-e2)+height)*sin(phi));


function [rad]=D2R(theta)
rad=theta*pi/180;


function [proj_x,proj_y,proj_z]=proj_3Dto2D(plane_x1,plane_y1,plane_z1,plane_x2,plane_y2,plane_z2,plane_x3,plane_y3,plane_z3,point_x0,point_y0,point_z0)
syms x y z;
A=[plane_x1,plane_y1,plane_z1];B=[plane_x2,plane_y2,plane_z2];C=[plane_x3,plane_y3,plane_z3];
D=[ones(4,1),[[x,y,z];A;B;C]];
detd=det(D);
cc=coeffs(detd);
if(length(cc)==4)
  coef_con=str2num(char(cc(1)));% plane function coef_x*x+coef_y*y+coef_z*z+coef_con=0
  coef_z=str2num(char(cc(2)));
  coef_y=str2num(char(cc(3)));
  coef_x=str2num(char(cc(4)));
elseif(length(cc)==3)
  coef_z=str2num(char(cc(1)));
  coef_y=str2num(char(cc(2)));
  coef_x=str2num(char(cc(3)));
  coef_con=0;
elseif(length(cc)==2)
  coef_y=str2num(char(cc(1)));
  coef_x=str2num(char(cc(2)));
  coef_z=0;coef_con=0;
end
proj_x=((coef_y^2+coef_z^2)*point_x0-coef_x*(coef_y*point_y0+coef_z*point_z0+coef_con))/(coef_x^2+coef_y^2+coef_z^2);
proj_y=((coef_x^2+coef_z^2)*point_y0-coef_y*(coef_x*point_x0+coef_z*point_z0+coef_con))/(coef_x^2+coef_y^2+coef_z^2);
proj_z=((coef_x^2+coef_y^2)*point_z0-coef_z*(coef_x*point_x0+coef_y*point_y0+coef_con))/(coef_x^2+coef_y^2+coef_z^2);
