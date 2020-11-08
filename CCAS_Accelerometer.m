clc; clear; close all;
format long g
tic 
data=xlsread('Inner_Tower - CX1_2612 - AC - 20201101_015108.csv');
%data=xlsread('Outer_Tower - CX1_2613 - AC - 20201101_022058.csv');
time=(data(:,3)-data(1,3))*86400; %Seconds
Accel_X=data(:,4)/9.81; %m/s^2
Accel_Y=data(:,5)/9.81; %m/s^2
Accel_Z=data(:,6)/9.81; %m/s^2
x=find(isnan(Accel_X)); %NAN entries for X,Y,Z direction (they are the same)
time(x)=[]; %Reomving time locations with NAN entries

%Removing NAN entries
Accel_X=Accel_X(~isnan(Accel_X));
Accel_Y=Accel_Y(~isnan(Accel_Y));
Accel_Z=Accel_Z(~isnan(Accel_Z));

%Normalizing data around 0
Accel_X=Accel_X-Accel_X(1);
Accel_Y=Accel_Y-Accel_Y(1);
Accel_Z=Accel_Z-Accel_Z(1);

%%% Data Analysis %%%
Vel_x_in=0;Vel_y_in=0;Vel_z_in=0;
Disp_x_in=0;Disp_y_in=0;Disp_z_in=0;
for i=1:length(time)-1
   temp_t=[time(i) time(i+1)];
   temp_ax=[Accel_X(i) Accel_X(i+1)];
   temp_ay=[Accel_Y(i) Accel_Y(i+1)];
   temp_az=[Accel_Z(i) Accel_Z(i+1)];
   temp_vx=[Vel_x_in trapz(temp_t, temp_ax)];
   temp_vy=[Vel_y_in trapz(temp_t, temp_ay)];
   temp_vz=[Vel_z_in trapz(temp_t, temp_az)];
   Disp_x(i)=trapz(temp_t, temp_vx);
   Disp_y(i)=trapz(temp_t, temp_vy);
   Disp_z(i)=trapz(temp_t, temp_vz);
   Vel_x_in=temp_vx(2);Vel_y_in=temp_vy(2);Vel_z_in=temp_vz(2);
   Disp_x_in=Disp_x(i);Disp_y_in=Disp_y(i);Disp_z_in=Disp_z(i);
end
Disp_x=Disp_x';Disp_y=Disp_y';Disp_z=Disp_z';

%%% Plotting Section %%%
figure()
hold on; grid on;
plot(time, Accel_X)
xlabel("Time (sec)")
ylabel("Acceleration (m/s^2)")
title("X Accleration vs Time")
figure()
hold on; grid on;
plot(time, Accel_Y)
xlabel("Time (sec)")
ylabel("Acceleration (m/s^2)")
title("Y Accleration vs Time")
figure()
hold on; grid on;
plot(time, Accel_Z)
xlabel("Time (sec)")
ylabel("Acceleration (m/s^2)")
title("Z Accleration vs Time")
figure()
hold on; grid on;
plot(time(2:end), Disp_x)
xlabel("Time (sec)")
ylabel("Position (m)")
title("X Displacement vs Time")
figure()
hold on; grid on;
plot(time(2:end), Disp_y)
xlabel("Time (sec)")
ylabel("Position (m)")
title("Y Displacement vs Time")
figure()
hold on; grid on;
plot(time(2:end), Disp_z)
xlabel("Time (sec)")
ylabel("Position (m)")
title("Z Displacement vs Time")
toc