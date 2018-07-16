function [Voc, S, K, Re] = TEADVinter_module(DB, Temp, A, t)
%TEADVinter Function to find material properties for the given temperature
%   Detailed explanation goes here

%% Read From Material Database

Temp_db  = DB(:,1)';
S_db     = DB(:,2)';
lamda_db = DB(:,3)';
sigma_db = DB(:,4)';

%% Target Value Designation Based on Temperature

Target = Temp(1:length(Temp) - 1) + diff(Temp)/2;

S     = interp1(Temp_db, S_db    , Target, 'linear', 'extrap') * 10^-6;
lamda = interp1(Temp_db, lamda_db, Target, 'linear', 'extrap');
sigma = interp1(Temp_db, sigma_db, Target, 'linear', 'extrap') * 10^2;

K  = (lamda  * A)  /  t; % W/K, Thermal Conductance
Re = t  ./ (sigma  * A); % Ohm, Electrical Resistance

%% Current Calculation (Power Generation Only)

Delta = -diff(Temp);

Voc  = sum(Delta .* S);
%Rtot = sum(Re) + R_L + R_c;

%I = Voc / Rtot;
% 
% Misc(1) = I;
% Misc(2) = Voc;
% Misc(3) = Rtot;

end