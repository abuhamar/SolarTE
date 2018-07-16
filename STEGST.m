% STEGST.m
% TE module sim for Concentrated Solar TE
% - Ameer, Kevin & Je-Hyeong (updated Summer 2018)

function STEGST(infile)
io = rpLib(infile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get input values from Rappture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Material properties database files
rpUtilsProgress(0, 'Reading Material Database');
choice = rpLibGetDouble(io,'input.phase(mat).group(horizontal).group(prop).choice.current');

DB_n(1,1) = 0;
DB_n(2,1) = 100;
DB_p = DB_n;

switch choice
    case 1
        % N-Type
        DB_n(:,2) = rpLibGetDouble(io,'input.phase(mat).string(nSC).current');
        DB_n(:,3) = rpLibGetDouble(io,'input.phase(mat).string(nTC).current');
        DB_n(:,4) = rpLibGetDouble(io,'input.phase(mat).string(nEC).current');
        
        % P-Type
        DB_p(:,2) = rpLibGetDouble(io,'input.phase(mat).string(pSC).current');
        DB_p(:,3) = rpLibGetDouble(io,'input.phase(mat).string(pTC).current');
        DB_p(:,4) = rpLibGetDouble(io,'input.phase(mat).string(pEC).current');        
    case 2
        DB_n = rpLibGetString(io,'input.phase(mat).group(horizontal).group(prop).group(DB).string(DBn).current');
        DB_n = str2num(DB_n);
        DB_p = rpLibGetString(io,'input.phase(mat).group(horizontal).group(prop).group(DB).string(DBp).current');
        DB_p = str2num(DB_p);
end

%% User input parameters ================
rpUtilsProgress(10, 'Preparing Input');

G = rpLibGetDouble(io,'input.phase(mat).string(G).current');            % Solar energy input from one sun [W/m2]
%Conct_array = [1 10:10:100];         % solar concentration [suns]
Conct_array = rpLibGetDouble(io,'input.phase(mat).string(Conct_array).current');;
%Tau_glass = 1;   % transmittance through top glass 
Tau_glass = rpLibGetDouble(io,'input.phase(mat).string(Tau_glass).current');;   % transmittance through top glass 
Alpha_abs = rpLibGetDouble(io,'input.phase(mat).string(Alpha_abs).current');;   % Absorptance of the solar absorber
Emissivity_abs = rpLibGetDouble(io,'input.phase(mat).string(Emissivity_abs).current');;      % Emissivity of the absorber

MatDensity=7.7;         % material density g/cm^3
MatCostPerKg=26;       % $/kg

% T_Top: top temperature is now a variable
T_Bot = rpLibGetDouble(io,'input.phase(mat).string(T_Bot).current');

% TE element dimensions 
A_array=[2*2];     % mm^2f
A_ratio=1;              % ratio A_p/A_n
%L_array= [0.1:0.1:0.9, 1:1:20];     % thickness array to test [mm]
L_array= [3];

Area_abs = rpLibGetDouble(io,'input.phase(mat).string(Total_Area).current');       % Absorber area mm^2     
Total_Area = rpLibGetDouble(io,'input.phase(mat).string(Total_Area).current');     % Total TE module area mm^2 

%FillFactor_array = [0.1:0.1:0.9];
FillFactor_array = rpLibGetDouble(io,'input.phase(mat).string(FillFactor_array).current');

R_c = 0;         % Ohm-cm^2 (Specific contact resistance between TE element and metal electrodes)

% Conductive heat exchange coefficients at 'h'ot side and 'c'old side substartes
Psi_h_array = [1000:1000:10000];       % W/m^2-K
Psi_c_array = [100:100:1000, 1500:500:10000];       % W/m^2-K
%Psi_h_array = [8000];       % W/m^2-K
%Psi_c_array = [3000];       % W/m^2-K

% Convection heat transfer coefficient at the hot side of the receiving slab
%h_conv = 20;       % =0 assuming TE in high vacuum. [W/m^2-K]
h_conv = 20;       % =0 assuming TE in high vacuum. [W/m^2-K]

% Thermal conductivity of filler
%k_filler = 0.04;          % vacuum [W/mK]
k_filler = rpLibGetDouble(io,'input.phase(mat).string(k_filler).current');;          % vacuum [W/mK]
count = 1;
lbar  = 1;
%===================================================================
%%  INDEPENDENT VARIABLES
rpUtilsProgress(20, 'Processing Independent Variable');
choice = rpLibGetDouble(io,'input.phase(mat).group(horizontal).group(ind).choice(choice).current');

switch choice
    case 1
        % ELEMENT THICKNESS
        min = rpLibGetDouble(io,'input.phase(mat).group(horizontal).group(ind).group(thick).number(min).current');
        step = rpLibGetDouble(io,'input.phase(mat).group(horizontal).group(ind).group(thick).number(step).current');
        max = rpLibGetDouble(io,'input.phase(mat).group(horizontal).group(ind).group(thick).number(max).current');
        L_array = min: step: max;
		len = length(L_array);
    case 2
        % ELEMENT AREA
        min = rpLibGetDouble(io,'input.phase(mat).group(horizontal).group(ind).group(area).number(min).current');
        step = rpLibGetDouble(io,'input.phase(mat).group(horizontal).group(ind).group(area).number(step).current');
        max = rpLibGetDouble(io,'input.phase(mat).group(horizontal).group(ind).group(area).number(max).current');
        A_array = min: step: max;
		len = length(A_array);
    case 3
        % LOAD RESSISTANCE RESISTANCE
        min = rpLibGetDouble(io,'input.phase(mat).group(horizontal).group(ind).group(resist).number(min).current');
        step = rpLibGetDouble(io,'input.phase(mat).group(horizontal).group(ind).group(resist).number(step).current');
        max = rpLibGetDouble(io,'input.phase(mat).group(horizontal).group(ind).group(resist).number(max).current');
        R_L_array = min: step: max;
		len = length(R_L_array);
end



%% Run

R_c = 0;         % Ohm-cm^2 (Specific contact resistance) Not a user input! ignore contact resistance!
N = 20;             % number of segments of each leg Not a user input! (see 66, 67)

fprintf('-----------------Starting Simulation-----------------\n\n')
%fprintf(' %1.0f segments for each leg\n\n', N)

DATA=[];
DATAmax=[];

DATA_Thot(:,1)=L_array';
DATA_Tcold(:,1)=L_array';
DATA_Tdiff(:,1)=L_array';
DATA_Pout(:,1)=L_array';
DATA_Qin(:,1)=L_array';
DATA_Eff_TE(:,1)=L_array';
DATA_Eff_system(:,1)=L_array';
DATA_Pout_hcoeff=[];

% for index_Qin=1:length(Q_in_array),
%     Q_in = Q_in_array(index_Qin);

for index_Conct=1:length(Conct_array),
     Conct = Conct_array(index_Conct);

for index_FillFactor=1:length(FillFactor_array)
    FillFactor = FillFactor_array(index_FillFactor);
        
    for index_A=1:length(A_array),
        
        A_n=A_array(index_A);
        A_p=A_ratio*A_array(index_A);
    
        R_c_n = R_c/(A_n*1e-2);     % Contact resistance (Ohm)
        R_c_p = R_c/(A_p*1e-2);     % Contact resistance (Ohm)
        
        N_n = floor(Total_Area / ( A_n / FillFactor + A_p / FillFactor ));
        N_p = floor(Total_Area / ( A_n / FillFactor + A_p / FillFactor ));        

        % Thickness optimization routine
        for index_L=1:length(L_array),
            
            L_n= L_array(index_L);
            L_p= L_array(index_L);
            
            for index_psih=1:length(Psi_h_array)
                for index_psic=1: length(Psi_c_array)
                    Psi_h = Psi_h_array(index_psih);
                    Psi_c = Psi_c_array(index_psic);
                    
                    SolarInput = G*Conct*Tau_glass*Alpha_abs*(Area_abs*1e-6);        % in Watt
                    
                    [Temp_n, Temp_p, Re_int, R_L, I, V_out, P_out, Q_in_te, Q_rad, Q_conv, Eff_TE] = TEADV_mode4_module(SolarInput, Total_Area*1e-6, Emissivity_abs, h_conv, T_Bot, R_c_n, R_c_p, Psi_h, Psi_c, L_n*1e-3, L_p*1e-3, A_n*1e-6, A_p*1e-6, N, N_n, N_p, DB_n, DB_p, LOAD_MATCH, k_filler, FillFactor);

                    Eff_system = P_out/SolarInput*Tau_glass*Alpha_abs*100;
                    
                    T_top = (Temp_n(1)+Temp_p(1))/2;
                    T_hot = (Temp_n(2)+Temp_p(2))/2;
                    T_cold = (Temp_n(N+2)+Temp_p(N+2))/2;
%                     fprintf(' Solar concentration : %f suns\n', Conct);
%                     fprintf(' Leg thickness : %f mm\n', L_n);
%                     fprintf(' Leg area %1.0f mm^2\n', A_n);
%                     fprintf(' Total %1.0f n-type legs and %1.0f p-type legs\n', N_n, N_p);
%                     fprintf(' Top temperature (T_top) is : %f K (%f degC)\n', T_top, T_top-273);
%                     fprintf(' TE hot-side temperature (T_hot) is : %f K (%f degC)\n', T_hot, T_hot-273);
%                     fprintf(' TE cold-side temperature (T_cold) is : %f K (%f degC)\n', T_cold, T_cold-273);
%                     fprintf(' Internal resistance (R_int) is : %f Ohm\n', Re_int);
%                     fprintf(' Resulted current (I) is : %f Amp\n', I);
%                     fprintf(' Voltage output (V_out) is : %f V\n', V_out);
%                     fprintf(' Total solar input : %f Watt\n', SolarInput);
%                     fprintf(' Radiation loss : %f Watt\n', Q_rad);
%                     fprintf(' Total heat input to TE : %f Watt\n', Q_in_te);
%                     fprintf(' Power output (P_out) is : %f Watt\n', P_out);
%                     fprintf(' TE efficiency: %6.3f %% \n', Eff_TE);
%                     fprintf(' System efficiency: %6.3f %% \n\n', Eff_system);

                    %% Final results:
                    % Top side temperature, bottom side temperature, Area of n-type element, area of p-type element, load resistance, internal resistance, current, total voltage_out, Power output, Heat input, Efficieny 
%                     DATA=[DATA; Conct, Psi_h, Psi_c, Total_Area, FillFactor, h_conv, k_filler, T_top, T_hot, T_cold, A_n, L_n, R_L, Re_int, I, V_out, P_out, Q_in_te, Eff_TE, T_hot-T_cold, Eff_system];       % mV, mW 
% 
%                     DATA_Thot(index_L,index_FillFactor+1)=T_hot;
%                     DATA_Tcold(index_L,index_FillFactor+1)=T_cold;
%                     DATA_Tdiff(index_L,index_FillFactor+1)=T_hot-T_cold;
%                     DATA_Qin(index_L,index_FillFactor+1)=Q_in_te;
%                     DATA_Eff_TE(index_L,index_FillFactor+1)=Eff_TE;
%                     DATA_Eff_system(index_L,index_FillFactor+1)=Eff_system;
%                     DATA_PowerCost(index_L,index_FillFactor+1)=(N_n*L_n*A_n+N_p*L_p*A_p)*(MatDensity/1e3)*1e-3*MatCostPerKg./P_out;  % $/Watt
% 
%                     DATA_Pout_hcoeff(index_psic,index_psih)=P_out;
%                    DATA_PowerCost(index_psic,index_FillFactor+1)=(N_n*L_n*A_n+N_p*L_p*A_p)*(MatDensity/1e3)*1e-3*MatCostPerKg./P_out;  % $/Watt
                    
                end
            end
%            DATAmax=[DATAmax; T_Top, T_Bot, A_n, L_n, Re_int, max_Pout_R_L, max_Pout, max_Pout_Qin, max_Pout_Eff, max_Eff_R_L, max_Eff_Pout, max_Eff_Qin, max_Eff]; 
               

    end
    end
end

end

%% Save output values back to Rappture
rpUtilsProgress(80, 'Generating Output Plots');
switch choice
    case 1
        mat = L_array;
        rpLibPutString(io,'output.field(ncurve).about.xaxis.label','Distance from Top Surface (mm)',0);
        rpLibPutString(io,'output.field(ncurve).about.yaxis.label','Thickness (mm)',0);
        rpLibPutString(io,'output.field(pcurve).about.xaxis.label','Distance from Top Surface (mm)',0);
        rpLibPutString(io,'output.field(pcurve).about.yaxis.label','Thickness (mm)',0);
        rpLibPutString(io,'output.curve(P_out).xaxis.label','Thickness (mm)',0);
        rpLibPutString(io,'output.curve(Eff).xaxis.label','Thickness (mm)',0);
        rpLibPutString(io,'output.curve(Cur).xaxis.label','Thickness (mm)',0);
        rpLibPutString(io,'output.curve(V_out).xaxis.label','Thickness (mm)',0);
    case 2
        mat = A_array;
        rpLibPutString(io,'output.field(ncurve).about.xaxis.label','Distance from Top Surface (mm)',0);
        rpLibPutString(io,'output.field(ncurve).about.yaxis.label','Total Area (mm^2)',0);
        rpLibPutString(io,'output.field(pcurve).about.xaxis.label','Distance from Top Surface (mm)',0);
        rpLibPutString(io,'output.field(pcurve).about.yaxis.label','Total Area (mm^2)',0);
        rpLibPutString(io,'output.curve(P_out).xaxis.label','Total Area (mm^2)',0);
        rpLibPutString(io,'output.curve(Eff).xaxis.label','Total Area (mm^2)',0);
        rpLibPutString(io,'output.curve(Cur).xaxis.label','Total Area (mm^2)',0);
        rpLibPutString(io,'output.curve(V_out).xaxis.label','Total Area (mm^2)',0);
    case 3
        mat = R_L_array;
        rpLibPutString(io,'output.field(ncurve).about.xaxis.label','Distance from Top Surface (mm)',0);
        rpLibPutString(io,'output.field(ncurve).about.yaxis.label','Load Resistance (Ohm)',0);
        rpLibPutString(io,'output.field(pcurve).about.xaxis.label','Distance from Top Surface (mm)',0);
        rpLibPutString(io,'output.field(pcurve).about.yaxis.label','Load Resistance (Ohm)',0);
        rpLibPutString(io,'output.curve(P_out).xaxis.label','Load Resistance (Ohm)',0);
        rpLibPutString(io,'output.curve(Eff).xaxis.label','Load Resistance (Ohm)',0);
        rpLibPutString(io,'output.curve(Cur).xaxis.label','Load Resistance (Ohm)',0);
        rpLibPutString(io,'output.curve(V_out).xaxis.label','Load Resistance (Ohm)',0);
end
%% Output Plots
rpUtilsProgress(90, 'Displaying Output Plots');
% fprintf('--------------Printing Temperature Plot--------------\n\n')

xydata = [mat;DATA(:,12)'];
str = sprintf('%12g %12g\n', xydata);
rpLibPutString(io,'output.curve(P_out).component.xy',str,0);

xydata = [mat;DATA(:,14)'];
str = sprintf('%12g %12g\n', xydata);
rpLibPutString(io,'output.curve(Eff).component.xy',str,0);

xydata = [mat;DATA(:,10)'];
str = sprintf('%12g %12g\n', xydata);
rpLibPutString(io,'output.curve(Cur).component.xy',str,0);

xydata = [mat;DATA(:,11)'];
str = sprintf('%12g %12g\n', xydata);
rpLibPutString(io,'output.curve(V_out).component.xy',str,0);

vals = reshape(Mesh_n, (N+3) * length(mat), 1);
str = sprintf('%12g %12g\n', vals);
rpLibPutString(io,'output.field(ncurve).component.values',str,0);
rpLibPutString(io,'output.mesh(ncurve).grid.xaxis.numpoints',sprintf('%12g\n',N+3),0);
rpLibPutString(io,'output.mesh(ncurve).grid.yaxis.numpoints',sprintf('%12g\n',length(mat)),0);

vals = reshape(Mesh_p, (N+3) * length(mat), 1);
str = sprintf('%12g %12g\n', vals);
rpLibPutString(io,'output.field(pcurve).component.values',str,0);
rpLibPutString(io,'output.mesh(pcurve).grid.xaxis.numpoints',sprintf('%12g\n',N+3),0);
rpLibPutString(io,'output.mesh(pcurve).grid.yaxis.numpoints',sprintf('%12g\n',length(mat)),0);

xydata = [DB_n(:,1)';DB_n(:,2)'];
str = sprintf('%12g %12g\n', xydata);
rpLibPutString(io,'output.curve(nSC_Plot).component.xy',str,0);

xydata = [DB_n(:,1)';DB_n(:,3)'];
str = sprintf('%12g %12g\n', xydata);
rpLibPutString(io,'output.curve(nTC_Plot).component.xy',str,0);

xydata = [DB_n(:,1)';DB_n(:,4)'];
str = sprintf('%12g %12g\n', xydata);
rpLibPutString(io,'output.curve(nEC_Plot).component.xy',str,0);

xydata = [DB_p(:,1)';DB_p(:,2)'];
str = sprintf('%12g %12g\n', xydata);
rpLibPutString(io,'output.curve(pSC_Plot).component.xy',str,0);

xydata = [DB_p(:,1)';DB_p(:,3)'];
str = sprintf('%12g %12g\n', xydata);
rpLibPutString(io,'output.curve(pTC_Plot).component.xy',str,0);

xydata = [DB_p(:,1)';DB_p(:,4)'];
str = sprintf('%12g %12g\n', xydata);
rpLibPutString(io,'output.curve(pEC_Plot).component.xy',str,0);

rpUtilsProgress(100, 'Done');

rpLibResult(io);
quit;

fprintf('----------------Simulation Completed!----------------\n\n')

end

