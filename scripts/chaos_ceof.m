function [B,D,t,expvar] = chaos_ceof(iopt,window_north,window_south,window_west,window_east,eval_lat,eval_lon,eval_time);
% inputs:
% window_north = mask extent degrees latitude north
% window_south = degrees latitude south (assume negative for south input)
% window_west = degrees longtitude west
% window_east = degrees longtitude east
% eval_lat = latitude where to plot EOF
% eval_lon = longitude where to plot EOF

% outputs:
% B = map of the first mode
% D = map of masked data
% t = time span 
% expvar = the fraction of total variance "explained" by each EOF ie it has the form EXPVAR(N).

close all

% check window inputs, error messages
if window_south >= window_north
    disp("Check your latitudes. North should be greater than South. Negative values for southern hemisphere.")
    return
end
if window_west > 0.0 && window_east > 0.0 && window_west >= window_east
    disp("Check your longtitudes. East should be greater than West. Negative values for western hemisphere.")
    return
end
if window_west < 0.0 && window_east < 0.0 && window_west >= window_east
    disp("Check your longtitudes. East should be greater than West. Negative values for western hemisphere.")
    return
end


% Lat/Lon 
pi = 3.14159265359;
Nth = 100;
Nph = Nth*2;
dth = pi/Nth;
colat_rad = linspace(dth/2, pi-dth/2, Nth);
lon_rad = linspace(dth/2, 2*pi-dth/2, Nph);
longitude = (lon_rad*180/pi)-180;
latitude =  90-(colat_rad*180/pi);
[Lon_grid ,Lat_grid] = meshgrid(longitude,latitude);

index_eval_lat = find(abs(eval_lat - latitude) == min(abs(eval_lat - latitude)));
index_eval_lon = find(abs(eval_lon - longitude) == min(abs(eval_lon - longitude)));
index_window_east = find(abs(window_east - longitude) == min(abs(window_east - longitude)));
index_window_west = find(abs(window_west - longitude) == min(abs(window_west - longitude)));
index_window_north = find(abs(window_north - latitude) == min(abs(window_north - latitude)));
index_window_south = find(abs(window_south - latitude) == min(abs(window_south - latitude)));


% Chaos 6 time span
t=1998.5:0.5:2017.5;
index_eval_time = find(eval_time - t == min(abs(eval_time - t)));
[XX, YY, time_grid] = meshgrid(longitude, latitude, t);

%load Chaos 6 maps of Br, SV, and SA
Br_all = load('Br_Map.mat','-mat');
SV_all = load('SV_Map.mat','-mat');
SA_all = load('SA_Map.mat','-mat');
Br_all = Br_all.Br_Map;
SV_all = SV_all.SV_Map;
SA_all = SA_all.SA_Map;

% Mask out portions of the map
F=zeros(Nth,Nph);

% Logic for different possible windows
% 1) simple window and longitude window across the prime meridian
if  window_south < window_north && window_west < window_east
    F(Lon_grid>window_west & Lon_grid<window_east & Lat_grid<window_north & Lat_grid>window_south) = 1;
end
% 2) wrapped around in longitude window
if window_west > 0.0 && window_east < 0.0
    F(Lon_grid>window_west & Lat_grid<window_north & Lat_grid>window_south) = 1;
    F(Lon_grid<window_east & Lat_grid<window_north & Lat_grid>window_south) = 1;
end


% vectorize map2mat
[D] = map2mat(F,SA_all);

% complex eof
[e,pc,expvar]=calCeof(D,3,4);

% reconstruct first and second mode
Ypred1 = pc(1,:).' * conj(e(1,:));
Ypred2 = pc(2,:).' * conj(e(2,:));
expvar

%go back to maps to plot
[A] = mat2map(F,D);
[B] = mat2map(F,Ypred1);
[C] = mat2map(F,Ypred2);
%%

%define diverging colormap -- I used
%http://colorbrewer2.org/#type=diverging&scheme=PiYG&n=11 to get RGB values
pink_green = [142,1,82
197,27,125
222,119,174
241,182,218
253,224,239
247,247,247
230,245,208
184,225,134
127,188,65
77,146,33
39,100,25] / 255;

colorbar_range = 1500;
load coastlines

% plot
if (iopt ==1)
fig1=figure(1) ;  
set(fig1,'Units','centimeters','Position',[0 30 20 30])
subplot(3,1,1)
contourf(latitude,t,A(:,:,index_eval_lon),'LineColor', 'none');
xlabel('Latitude, Degrees')
ylabel('Time, Year')
title(['Data, Longitude = ' num2str(longitude(index_eval_lon))])
colormap(cool)
caxis([-colorbar_range colorbar_range])
colorbar
subplot(3,1,2)
contourf(latitude,t,real(B(:,:,index_eval_lon)),'LineColor', 'none');
xlabel('Latitude')
ylabel('Time, Year')
title('First Mode')
colormap(cool)
caxis([-colorbar_range colorbar_range])
colorbar
subplot(3,1,3)
contourf(latitude,t,real(C(:,:,index_eval_lon)),'LineColor', 'none');
xlabel('Latitude, Degrees')
ylabel('Time, Year')
title('Second Mode')
colormap(cool)
caxis([-colorbar_range colorbar_range])
colorbar
filename=['EOF_latitude_' num2str(longitude(index_eval_lon)) '.eps'];
%saveas(fig1,filename,'epsc')

fig2=figure(2) ;
set(fig2,'Units','centimeters','Position',[0 30 20 30])
subplot(3,1,1)
contourf(longitude,t,permute(A(:,index_eval_lat,:),[1 3 2]),'LineColor', 'none');
xlabel('Longitude, Degrees')
ylabel('Time, Year')
title(['Data, Latitude = ' num2str(latitude(index_eval_lat))])
colormap(cool)
caxis([-colorbar_range colorbar_range])
colorbar
subplot(3,1,2)
contourf(longitude,t,permute(real(B(:,index_eval_lat,:)),[1 3 2]),'LineColor', 'none');
xlabel('Longitude, Degrees')
ylabel('Time, Year')
title('First Mode')
colormap(cool)
caxis([-colorbar_range colorbar_range])
colorbar
subplot(3,1,3)
contourf(longitude,t,permute(real(C(:,index_eval_lat,:)),[1 3 2]),'LineColor', 'none');
xlabel('Longitude, Degrees')
ylabel('Time, Year')
title('Second Mode')
colormap(cool)
caxis([-colorbar_range colorbar_range])
colorbar
filename=['EOF_longitude_' num2str(latitude(index_eval_lat)) '.eps'];
saveas(fig2,filename,'epsc')

fig3=figure(3) ;
set(fig3,'Units','centimeters','Position',[0 30 20 30])
subplot(3,1,1)
contourf(longitude,latitude,permute(A(index_eval_time,:,:),[2 3 1]),'LineColor', 'none');
hold on
plot(coastlon,coastlat,'k')
xlabel('Longitude, Degrees')
ylabel('Latitude, Degrees')
title(['Data, year = ' num2str(t(index_eval_time))])
colormap(cool)
caxis([-colorbar_range colorbar_range])
colorbar
subplot(3,1,2)
contourf(longitude,latitude,permute(real(B(index_eval_time,:,:)),[2 3 1]),'LineColor', 'none');
hold on
plot(coastlon,coastlat,'k')
xlabel('Longitude, Degrees')
ylabel('Latitude, Degrees')
title('First Mode')
colormap(cool)
caxis([-colorbar_range colorbar_range])
colorbar
subplot(3,1,3)
contourf(longitude,latitude,permute(real(C(index_eval_time,:,:)),[2 3 1]),'LineColor', 'none');
hold on
plot(coastlon,coastlat,'k')
xlabel('Longitude, Degrees')
ylabel('Latitude, Degrees')
title('Second Mode')
colormap(cool)
caxis([-colorbar_range colorbar_range])
colorbar
filename=['EOF_map_' num2str(t(index_eval_time)) '.eps'];
%saveas(fig3,filename,'epsc')

fig5=figure(5);
plot(t(1:end),real(pc(1,1:end)),'LineWidth',2)
hold on
plot(t(1:end),real(pc(2,1:end)),'LineWidth',2)
plot(t(1:end),real(pc(3,1:end)),'LineWidth',2)
legend('mode 1','mode 2', 'mode 3')
plot([t(1) t(end)],[0 0],'k--')
title('Principle Components')
xlabel('Time, Year')
%saveas(fig5,'principle_component_time.eps','epsc')
end


fig6=figure(6);
isosurface(XX, YY, time_grid, permute(real(B),[2 3 1]),800)
axis equal
grid on
xlabel('Longitude')
ylabel('Latitude')
zlabel('time')
title('Mode 1')

% spatial phase
if (iopt==2)
    kx1 = atan2(imag(B(index_eval_time,index_eval_lat,:)),real(B(index_eval_time,index_eval_lat,:)));
    kx2 = atan2(imag(C(index_eval_time,index_eval_lat,:)),real(C(index_eval_time,index_eval_lat,:)));
    fig1=figure(1);
    plot(longitude,squeeze(kx2),'LineWidth',2)
    hold on
    plot([longitude(index_window_west) longitude(index_window_east)],[0 0],'k--')
    plot([longitude(index_window_west) longitude(index_window_east)],[pi/2 pi/2],'k--')
    plot([longitude(index_window_west) longitude(index_window_east)],[pi pi],'k--')
    plot([longitude(index_window_west) longitude(index_window_east)],[-pi/2 -pi/2],'k--')
    plot([longitude(index_window_west) longitude(index_window_east)],[-pi -pi],'k--')
    title(['Spatial Phase, Latitude = ' num2str(latitude(index_eval_lat)) ' Time = ' num2str(t(index_eval_time))])
    xlabel 'Longitude, Degrees'
    ylabel 'Spatial Phase kx'
    filename=['Spatial_Phase_Mode1_year' num2str(t(index_eval_time)) '_lat' num2str(latitude(index_eval_lat)) '.eps'];
    %saveas(fig1,filename,'epsc')
    
    kx1 = atan2(imag(B(:,index_eval_lat,:)),real(B(:,index_eval_lat,:)));
    kx2 = atan2(imag(C(:,index_eval_lat,:)),real(C(:,index_eval_lat,:)));
    figure(2)
    cm = parula(length(t));
    for i=1:length(t);
        plot(longitude,squeeze(kx1(i,:,:)),'o','MarkerEdgeColor',cm(i,:));
        hold on
    end
    colormap(parula)
    cb=colorbar('Limits',[0 1],'Ticks',[0,1],'TickLabels',...
        {num2str(t(1)),num2str(t(end))});
    cb.Label.String = 'Time, year';
    plot([longitude(index_window_west) longitude(index_window_east)],[0 0],'k--')
    plot([longitude(index_window_west) longitude(index_window_east)],[pi/2 pi/2],'k--')
    plot([longitude(index_window_west) longitude(index_window_east)],[pi pi],'k--')
    plot([longitude(index_window_west) longitude(index_window_east)],[-pi/2 -pi/2],'k--')
    plot([longitude(index_window_west) longitude(index_window_east)],[-pi -pi],'k--')
    title(['Spatial Phase, Latitude = ' num2str(latitude(index_eval_lat)) ])
    xlabel 'Longitude, Degrees'
    ylabel 'Spatial Phase kx'
 
end

% Temporal Phase
if (iopt ==3)
    wt1 = atan2(imag(B(:,index_eval_lat,index_eval_lon)),real(B(:,index_eval_lat,index_eval_lon)));
    wt2 = atan2(imag(C(:,index_eval_lat,index_eval_lon)),real(C(:,index_eval_lat,index_eval_lon)));
    figure(1)
    plot(t,wt1,'LineWidth',2)
    hold on
    plot([t(1) t(end)],[0 0],'k--')
    plot([t(1) t(end)],[pi/2 pi/2],'k--')
    plot([t(1) t(end)],[pi pi],'k--')
    plot([t(1) t(end)],[-pi/2 -pi/2],'k--')
    plot([t(1) t(end)],[-pi -pi],'k--')
    title(['Temporal Phase, Latitude = ' num2str(latitude(index_eval_lat)) ' Longitude = ' num2str(longitude(index_eval_lon))])
    xlabel 'Time, years'
    ylabel 'Temporal Phase wt'
    filename=['Temporal_Phase_Mode1_lat' num2str(latitude(index_eval_lat)) '_lon' num2str(longitude(index_eval_lon)) '.eps'];
    %saveas(fig1,filename,'epsc')
    
    
    wt1 = atan2(imag(B(:,index_eval_lat,:)),real(B(:,index_eval_lat,:)));
    wt2 = atan2(imag(C(:,index_eval_lat,:)),real(C(:,index_eval_lat,:)));
    
    figure(2)
    cm = parula(length(longitude(index_window_west:index_window_east)));
    for i=1:length(longitude(index_window_west:index_window_east));
        plot(t,wt1(:,:,i+index_window_west-1),'o','MarkerEdgeColor',cm(i,:));
        hold on
    end
    colormap(parula)
    cb=colorbar('Limits',[0 1],'Ticks',[0,1],'TickLabels',...
        {num2str(longitude(index_window_west)),num2str(longitude(index_window_east))});
    cb.Label.String = 'Longitude, Degree';
    plot([t(1) t(end)],[0 0],'k--')
    plot([t(1) t(end)],[pi/2 pi/2],'k--')
    plot([t(1) t(end)],[pi pi],'k--')
    plot([t(1) t(end)],[-pi/2 -pi/2],'k--')
    plot([t(1) t(end)],[-pi -pi],'k--')
    title(['Temporal Phase, Latitude = ' num2str(latitude(index_eval_lat))])
    xlabel 'Time, years'
    ylabel 'Temporal Phase wt'
end


