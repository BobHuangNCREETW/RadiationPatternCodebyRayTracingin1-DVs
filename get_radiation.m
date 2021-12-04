load taiwan_coast84.mat;
num=2;
if(num==1)
  load tmp1_database_TW_DeepCru.mat;
  tmpnm='database_TW_DeepCru.csv';
elseif(num==2)
  load tmp1_database_TW_ShaCru.mat;
  tmpnm='database_TW_ShaCru.csv';
elseif(num==3)
  load tmp1_database_TW_subInter.mat;
  tmpnm='database_TW_subInter.csv';
elseif(num==4)
  load tmp1_database_TW_subIntra.mat;
  tmpnm='database_TW_subIntra.csv';
end
%% find sta & evt cordinate
uni_stanm=unique(sta);
for i=1:length(uni_stanm)
  index_same=find(strcmp(uni_stanm{i},sta)==1);
  uni_stalon(i)=STAlon(index_same(1));uni_stalat(i)=STAlat(index_same(1));
end
uni_EQnm=unique(EQID);
for i=1:length(uni_EQnm)
  index_same=find(strcmp(uni_EQnm{i},EQID)==1);
  uni_EVlon(i)=EVlon(index_same(1));uni_EVlat(i)=EVlat(index_same(1));  
end
%%
layered_Vs_1D_CWB=[1.96 2.62 3.03 3.35 3.61 3.71 3.95 4.21 4.49 4.68 4.72 4.79 4.8 4.74 4.86 4.92 5.49];
layered_Vp_1D_CWB=[3.48 4.48 5.25 5.83 6.21 6.41 6.83 7.29 7.77 8.05 8.16 8.34 8.35 8.2 8.4 8.51 8.7];
PoissonRatio=0.5*((layered_Vp_1D_CWB.^2)-2.*(layered_Vs_1D_CWB.^2))./(layered_Vp_1D_CWB.^2.-(layered_Vs_1D_CWB.^2));
layered_D_lower=[2 4 9 13 17 25 30 35 50 70 90 110 140 170 200 240 280];
AZIMUTH = (0:4:360)*pi/180;
TAKEOFF = (0:2:180)*pi/180;
[AZIMUTH,TAKEOFF] = meshgrid(AZIMUTH,TAKEOFF);
EQID_com='temp';
fid=fopen(tmpnm,'w');
fprintf(fid,'%s\n','smoothed_FAS_ID,fileID,EQID,EVlat,EVlon,EVdep,EVmw,EVml,faulttype,Ztor,sta,STAlat,STAlon,Vs30,Z1_0,FASID,Repi,Rhypo,RJB,Rrup,Comp1_Low-Cut_Freq,Comp1_High-Cut_Freq,Comp2_Low-Cut_Freq,Comp2_High-Cut_Freq,Comp3_Low-Cut_Freq,Comp3_High-Cut_Freq,PGA_comp1,PGA_comp2,PGA_comp3,PGV_comp1,PGV_comp2,PGV_comp3,S-DS-EAS_0.1Hz,S-DS-EAS_0.2Hz,S-DS-EAS_0.33Hz,S-DS-EAS_0.5Hz,S-DS-EAS_1.0Hz,S-DS-EAS_1.33Hz,S-DS-EAS_1.67Hz,S-DS-EAS_2Hz,S-DS-EAS_2.5Hz,S-DS-EAS_3.3Hz,S-DS-EAS_5Hz,S-DS-EAS_10Hz,takeoff_angle,Ray_path_filenm,Strike,Dip,Rake,Rad_Coef_S,Azimuth,Rad_plot_nm');
for i=1:length(smoothed_FAS_ID)
  [Repi_cal(i),az_cal(i)]=distance(STAlat(i),STAlon(i),EVlat(i),EVlon(i));Repi_cal(i)=Repi_cal(i)*(pi/180)*6371;
  az_cal(i)=az_cal(i)+180; %az is azimuth from source to station N=0;E=90;S=180;W=270;
  if(az_cal(i)>360)
    az_cal(i)=az_cal(i)-360;
  end
  index_start_tmp=find((layered_D_lower-EVdep(i))>0);tmp1=find(min(layered_D_lower(index_start_tmp)-EVdep(i)));index_start=index_start_tmp(tmp1);
  sigma_poisson=PoissonRatio(index_start);
  if(strcmp(EQID{i},EQID_com)==0)
    Strike_radiation=Strike(i);Dip_radiation=Dip(i);Rake_radiation=Rake(i);
    [GP, GS, GSH, GSV] = rpgen(Strike_radiation,Dip_radiation,Rake_radiation,0,sigma_poisson,TAKEOFF*180/pi, AZIMUTH*180/pi);
    scale = abs(GS);YP = scale.*cos(AZIMUTH) .* sin(TAKEOFF);
    XP = scale.*sin(AZIMUTH) .* sin(TAKEOFF);ZP = scale.*-cos(TAKEOFF);    
    EQID_com=EQID{i};
  end
  AZIMUTH_tar=az_cal(i)*pi/180;TAKEOFF_tar=takeoff_angle(i)*pi/180;
  [GP_tar, GS_tar, GSH_tar, GSV_tar] = rpgen(Strike(i),Dip(i),Rake(i),0,sigma_poisson, TAKEOFF_tar*180/pi, AZIMUTH_tar*180/pi);
  scale_tar = abs(GS_tar);YP_tar = scale_tar.*cos(AZIMUTH_tar) .* sin(TAKEOFF_tar);
  XP_tar = scale_tar.*sin(AZIMUTH_tar) .* sin(TAKEOFF_tar);ZP_tar = scale_tar.*-cos(TAKEOFF_tar);
%% Read FAS ascii file to plot it
  pwd_EAS='/FAS_smooth_ascii_EAS/';
  fid_EAS=fopen([pwd_EAS,smoothed_FAS_ID{i},'.fas'],'r');
  nu=fgetl(fid_EAS);
  tmp=textscan(fid_EAS,'%f %f %f\n');
  freq=tmp{1};FAS_EAS_DS=tmp{2};sm_FAS_EAS=tmp{3};
  fclose(fid_EAS);
  HP_Hori=max(Comp1_LowCut_Freq(i),Comp2_LowCut_Freq(i));
  ymax_Hori=10^(round(log10(max(sm_FAS_EAS)))+2);
  ymin_Hori=10^(floor(log10(min(sm_FAS_EAS)))-1);
%%  
  rad_outnm=['Radiation_',EQID{i},'_',sta{i}];
  plotYN='Y';
%    plotYN='N';
  if(plotYN=='Y')
%% find path from previous pbr calculation
    Ray_path_file=['Ray_path_all/',Ray_path_filenm{i}];
    fid_ray=fopen(Ray_path_file,'r');
    nu=fgetl(fid_ray);
    route_point_num_tmp=fgetl(fid_ray);
    [route_point_num_S,treval_time_S]=strread(route_point_num_tmp,'%f %f');
    source_tmp=fgetl(fid_ray);[evla,evlo,evdp]=strread(source_tmp,'%f %f %f');
    sta_tmp=fgetl(fid_ray);[stla,stlo,stdp]=strread(sta_tmp,'%f %f %f');
    cordi_route=fscanf(fid_ray,'%f',[route_point_num_S 3]);
    route_dep_S=cordi_route(:,1);route_lat_S=cordi_route(:,2);route_lon_S=cordi_route(:,3);
    fclose(fid_ray);
%%
    figure('Visible','Off');
%    figure(1);    
    H1=subplot(3,2,1);
      plot(uni_stalon,uni_stalat,'v','markerfacecolor',[0.6 0.6 0.6],'markersize',4);hold on;
      plot(uni_EVlon,uni_EVlat,'o','markerfacecolor',[0.6 0.6 0.6],'markersize',4);
      plot(lon_taiwan_coast,lat_taiwan_coast,'k- ');
      plot(STAlon(i),STAlat(i),'v','linewidth',2,'markeredgecolor','k','markerfacecolor','w','markersize',7);hold on;
      plot(EVlon(i),EVlat(i),'o','linewidth',2,'markeredgecolor','k','markerfacecolor','r','markersize',7.5);
      axis([119.4 124 21.4 25.5]);
      axis equal;
      Rrup_title=(ceil(Rrup(i)*100))/100;Vs30_title=(ceil(Vs30(i)*10))/10;
      title(['Mw:',num2str(EVmw(i)),' Rrup:',num2str(Rrup_title),'km'],'fontsize',10,'Fontweight','bold');
    H2=subplot(3,2,[2,4]);
      loglog(freq,sm_FAS_EAS,'- r','linewidth',1.5);hold on;
      xmin=0.01;xmax=100;
      line([HP_Hori HP_Hori],[ymin_Hori ymax_Hori],'linewidth',2,'color','k','linestyle','--');
      axis([xmin,xmax,ymin_Hori,ymax_Hori]);
      Low_useable_cor=HP_Hori;High_useable_cor=50;
      fill_L_x1=[xmin Low_useable_cor Low_useable_cor xmin];fill_L_y1=[ymax_Hori ymax_Hori ymin_Hori ymin_Hori];
      fill(fill_L_x1',fill_L_y1,[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8]);
      fill_H_x1=[High_useable_cor xmax xmax High_useable_cor];fill_H_y1=[ymax_Hori ymax_Hori ymin_Hori ymin_Hori];
      fill(fill_H_x1',fill_H_y1,[0.8 0.8 0.8],'edgecolor',[0.8 0.8 0.8]);
      line([Comp1_LowCut_Freq(i) Comp1_LowCut_Freq(i)],[ymin_Hori ymax_Hori],'linewidth',2,'color','k','linestyle','--');
      line([Comp2_LowCut_Freq(i) Comp2_LowCut_Freq(i)],[ymin_Hori ymax_Hori],'linewidth',2,'color','k','linestyle','--');
      line([High_useable_cor High_useable_cor],[ymin_Hori ymax_Hori],'linewidth',2,'color','k','linestyle','--');
      loglog(freq,sm_FAS_EAS,'- r','linewidth',1.5);
      index_EQnm_dash=find(EQID{i}=='_');EQnm_title=EQID{i};EQnm_title(index_EQnm_dash)=' ';
      title([EQnm_title,' ',sta{i},' Vs30:',num2str(Vs30_title),'m/s'],'fontsize',12,'Fontweight','bold');
      xlabel('freq(Hz)','fontsize',10,'fontweight','bold');ylabel('FAS (g-s)','fontsize',11,'fontweight','bold');
      legend('S-DS EAS','bandpass F_0',3);    
    H3=subplot(3,2,3);
      plot3(route_lon_S,route_lat_S,-1*route_dep_S,'linewidth',2);hold on;
      plot3(route_lon_S(1),route_lat_S(1),-1*route_dep_S(1),'o','linewidth',2,'markeredgecolor','k','markerfacecolor','r','markersize',7.5);
      plot3(route_lon_S(length(route_lon_S)),route_lat_S(length(route_lat_S)),-1*route_dep_S(length(route_dep_S)),'v','linewidth',2,'markeredgecolor','k','markerfacecolor','w','markersize',7);
      grid on;
      az_cal_title=(ceil(az_cal(i)*10))/10;
      title(['Ray path, Takeoff angle:',num2str(takeoff_angle(i)),char(176),sprintf('\n'),'Azimuth:',num2str(az_cal_title),char(176)],'fontsize',10,'Fontweight','bold');
    H4=subplot(3,2,5);
      surf(XP,YP,ZP,GS,'EdgeColor','none','FaceAlpha',0.7,'AmbientStrength',0.7,'EdgeColor','k','EdgeAlpha',0.1);hold on;
      shading interp; %get rid of grid edge lines
      plot3(XP_tar,YP_tar,ZP_tar,'v','linewidth',2,'markeredgecolor','k','markerfacecolor','w','markersize',8);
      title(['S-wave Radiation Pattern',sprintf('\n'),'Strike:',num2str(Strike(i)),',Dip:',num2str(Dip(i)),',Rake:',num2str(Rake(i))],'fontsize',10,'Fontweight','bold');
      xlabel('EW','fontsize',8,'Fontweight','bold');ylabel('NS','fontsize',8,'Fontweight','bold');zlabel('UD','fontsize',8,'Fontweight','bold');
    H5=subplot(3,2,6);
      surf(XP,YP,ZP,GS,'EdgeColor','none','FaceAlpha',0.7,'AmbientStrength',0.7,'EdgeColor','k','EdgeAlpha',0.1);hold on;
      shading interp; %get rid of grid edge lines
      plot3(XP_tar,YP_tar,max(max(ZP)),'v','linewidth',2,'markeredgecolor','k','markerfacecolor','w','markersize',8);
      xlabel('EW','fontsize',8,'Fontweight','bold');ylabel('NS','fontsize',8,'Fontweight','bold');
      GS_tar_title=(ceil(GS_tar*100))/100;
      title(['Coefficient:',num2str(GS_tar_title),sprintf('\n'),'Map View'],'fontsize',10,'Fontweight','bold');
      axis equal;
      view(0,90); %view from top=0 degree & horizontal=90 degree
      colorbar;
      axis([-1,1,-1,1])

      set(H2,'xtick',[0.01,0.1,0.2,0.5,1,2,5,10.0,20.0,50.0,100.0],'xticklabel',[{0.01},{0.1},{0.2},{0.5},{1.0},{2.0},{5.0},{10},{20},{50},{100}]);
      set(H2,'linewidth',2,'Fontweight','bold');set(H2,'layer','top');
      set(H1,'position',[0.07 0.54 0.25 0.4],'xlim',[119.4 124],'ylim',[21.4 25.5],'xtick',[119:1:124],'ytick',[21:1:26]);    
      set(H2,'position',[0.45 0.48 0.45 0.45]);set(H1,'linewidth',2,'Fontweight','bold');set(H1,'layer','top');
      set(H3,'position',[0.05 0.15 0.28 0.28]);set(H3,'linewidth',2,'Fontweight','bold');set(H3,'layer','top');
      set(H4,'position',[0.41 0.09 0.23 0.23]);set(H4,'linewidth',2,'Fontweight','bold');set(H4,'layer','top');
      set(H5,'position',[0.68 0.09 0.23 0.23]);set(H5,'linewidth',2,'Fontweight','bold');set(H5,'layer','top');
    saveas(gcf, rad_outnm , 'png');
    unix(['mv ',rad_outnm,'.png Rad_fig']);
    close all;
  end
  fprintf(fid,'%s,%s,%s,%9.5f,%10.5f,%8.4f,%4.2f,%4.2f,%s,%10.5f,%s,%9.5f,%10.5f,%7.2f,%10.4f,%s,%9.4f,%9.4f,%9.4f,%9.4f,%9.6f,%9.6f,%9.6f,%9.6f,%9.6f,%9.6f,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%6.2f,%s,%6.1f,%6.1f,%6.1f,%6.4f,%6.2f,%s\n',smoothed_FAS_ID{i},fileID{i},EQID{i},EVlat(i),EVlon(i),EVdep(i),EVmw(i),EVml(i),faulttype{i},Ztor(i),sta{i},STAlat(i),STAlon(i),Vs30(i),Z1_0(i),FASID{i},Repi(i),Rhypo(i),RJB(i),Rrup(i),Comp1_LowCut_Freq(i),Comp1_HighCut_Freq(i),Comp2_LowCut_Freq(i),Comp2_HighCut_Freq(i),Comp3_LowCut_Freq(i),Comp3_HighCut_Freq(i),PGA_comp1(i),PGA_comp2(i),PGA_comp3(i),PGV_comp1(i),PGV_comp2(i),PGV_comp3(i),SDSEAS_01Hz(i),SDSEAS_02Hz(i),SDSEAS_033Hz(i),SDSEAS_05Hz(i),SDSEAS_10Hz(i),SDSEAS_133Hz(i),SDSEAS_167Hz(i),SDSEAS_2Hz(i),SDSEAS_25Hz(i),SDSEAS_33Hz(i),SDSEAS_5Hz(i),SDSEAS_10Hz1(i),takeoff_angle(i),Ray_path_filenm{i},Strike(i),Dip(i),Rake(i),GS_tar,az_cal(i),rad_outnm);
end

