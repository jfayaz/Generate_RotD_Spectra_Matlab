function [rotd_combine] = RotDxx_Calculations( data1, data2, dt1, dt2, combine_index, T, damping )
% Author:  Sarah Azar
% Revised by: Jawad Fayaz

% INPUT:
% data1: as recorded acceleration 1
% data2: as recorded acceleration 2
% dt1, dt2: time steps of data1 and data2
% combine_index: 0: RotD00 and its corresponding values; 50: RotD50; 100: RotD100

dy    = 100*ones(1,length(T));
alpha = 0; % strain hardening ratio - not used
g     = 9.81; %m/s
RotD180 = zeros(180, length(T));
dt    = 0.005;                 % Dt    = analyis time step

% Time step for analysis is taken here equal to 0.005s
[~,~,H1] = Bilinear_Newmark_v031218( T, damping, dy, alpha, g*data1, dt1, dt );
[~,~,H2] = Bilinear_Newmark_v031218( T, damping, dy, alpha, g*data2, dt2, dt );

% should be the relative displacement time series (not spectrum)
RD1=H1.d'; % matrix [period by data points]
RD2=H2.d'; % matrix [period by data points]
omega = 2*pi./T;

%% Rotation Calculations
for per_index=1:length(T) % loop over periods
    length_min = min(length(RD1(per_index,:)), length(RD2(per_index,:)));
    Rot_rd1 = RD1(per_index,(1:length_min));
    Rot_rd2 = RD2(per_index,(1:length_min));
    for theta=1:180
        RS1 = Rot_rd1*cos(theta/180*pi) + Rot_rd2*sin(theta/180*pi);
        RotD180(theta, per_index) = max(abs(RS1))*omega(per_index)^2/g;
    end
end


RotD00  = min(RotD180);
RotD100 = max(RotD180);
RotD50  = median(RotD180);

switch combine_index
    case 0
         rotd_combine = RotD00;
    case 50
         rotd_combine = RotD50;
    case 100
         rotd_combine = RotD100;
end
        
end

