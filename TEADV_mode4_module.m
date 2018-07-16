function [Temp_n, Temp_p, Re_int_tot, R_L, I, V_out, P_out, Q_in2, Q_rad, Q_conv, Eff_TE] = TEADV_mode4_module(SolarInput,  Total_Area, eps1, h_conv, T_Bot, R_c_n, R_c_p, Psi_h, Psi_c, L_n, L_p, A_n, A_p, N, N_n, N_p, DB_n, DB_p, LOAD_MATCH, k_f, FF)
%TEADV Summary of this function goes here
%   Power generation mode 4: Solar Thermoelectric module

%close
%clc
sigma = 5.670e-8;    % Stefan-Boltzmann constant [W m^-2 K^-4]

TempTolerance = 1e-4;       % ratio
TempTol_hotside = 0.001;     % i.e. 0.5 % tolerance

t_n = L_n / N; % m (length of each segment of n-leg)
t_p = L_p / N; % m (length of each segment of p-leg)

%A_fraction_n = A_n / (N_n*A_n+ N_p*A_p);
%A_fraction_p = A_p / (N_n*A_n+ N_p*A_p);

%% Initial conditions

lamda_n = interp1(DB_n(:,1), DB_n(:,3), T_Bot, 'linear', 'extrap');
lamda_p = interp1(DB_p(:,1), DB_p(:,3), T_Bot, 'linear', 'extrap');

T_amb = T_Bot;

% Initial guess of T_Top
go_initial=1;
T_Top= T_Bot;
T_step = 10;
while (go_initial==1) 

    T_Top = T_Top+T_step;
    
    k_avg_n = interp1(DB_n(:,1), DB_n(:,3), (T_Bot+T_Top)/2, 'linear', 'extrap');
    k_avg_p = interp1(DB_p(:,1), DB_p(:,3), (T_Bot+T_Top)/2, 'linear', 'extrap');
    k_avg= (k_avg_n + k_avg_p)/2;

    T_Top2 = T_Bot + (SolarInput-Total_Area*sigma*eps1*(T_Top^4 - T_amb^4)-h_conv*Total_Area*(T_Top-T_amb))*(1/Psi_h + L_n/k_avg/FF + 1/Psi_c)/Total_Area;
    
    if (T_Top2 < T_Top || T_Top >= 800)
        go_initial =0;
    end

end 
%T_Top

% Initial temperature profiles
Temp_n = linspace(T_Top, T_Bot, N+1)';
Temp_p = linspace(T_Top, T_Bot, N+1)';
Temp_n(N+2)=T_Top;
Temp_p(N+2)=T_Top;
    
%% CALCULATION (2D loops)
go_qin_div=1;
iter_qin_div=0;

    
options = optimset('Display','off','TolFun',1e-18);


% Lateral heat exchange in the top heat exchanger to balance n-type and
% p-type hot side temperatures; initially zero
Q_lateral =0;
 
while go_qin_div==1

    iter_qin_div = iter_qin_div + 1;
    
    go   = 1;
    iter(iter_qin_div) = 0;

       
    while go == 1

        iter(iter_qin_div) = iter(iter_qin_div) + 1;
        
        [Voc_n, S_n, K_n, Re_n] = TEADVinter_module(DB_n, Temp_n(1:N+1), A_n, t_n);
        [Voc_p, S_p, K_p, Re_p] = TEADVinter_module(DB_p, Temp_p(1:N+1), A_p, t_p);

        Re_int_tot = N_n*( sum(Re_n) + 2*R_c_n ) + N_p*( sum(Re_p)+ 2*R_c_p);
        if (LOAD_MATCH)
            R_L = Re_int_tot;       
        end

        Re_tot = Re_int_tot + R_L;
        Voc_tot = -N_n*Voc_n + N_p*Voc_p;

        I = Voc_tot / Re_tot; 

        % Total heat input
        
        % n-leg
        f_n = @(x_n) TEADVfun_mode4(x_n, -I, S_n, K_n, Re_n, R_c_n, Total_Area, SolarInput, eps1, h_conv, T_Bot, Psi_h, Psi_c, N, A_n, Q_lateral, k_f, FF, L_n);
        [x_n] = fsolve(f_n, Temp_n, options);

        Q_in_n = A_n/FF*Psi_h*(x_n(N+2)-x_n(1))-(1-FF)/FF*A_n/L_n*k_f*(x_n(1)-x_n(N+1))+ Q_lateral;       

        % p-leg
        f_p = @(x_p) TEADVfun_mode4(x_p, I, S_p, K_p, Re_p, R_c_p, Total_Area, SolarInput, eps1, h_conv, T_Bot, Psi_h, Psi_c, N, A_p, -Q_lateral, k_f, FF, L_p);
        [x_p] = fsolve(f_p, Temp_p, options);

        Q_in_p = A_p/FF*Psi_h*(x_p(N+2)-x_p(1))-(1-FF)/FF*A_p/L_p*k_f*(x_p(1)-x_p(N+1))- Q_lateral;        
               
        if ( ( max(abs((Temp_n - x_n) ./ Temp_n)) < TempTolerance ) && ( max(abs((Temp_p - x_p) ./ Temp_p)) < TempTolerance ) ) 
                go = 0;
        end
        
        %max(abs((Temp_n - x_n) ./ Temp_n))
        %max(abs((Temp_p - x_p) ./ Temp_p))
        
%        x_n(N+2)
%        x_p(N+2)
        
        Temp_n = x_n;
        Temp_p = x_p;        
        
    end     % end of the first while

    
    
    % Re-division of Qin and re-guessing of T_top
    
    T_Top = ( Temp_n(N+2) + Temp_p(N+2) ) / 2;
     
    if abs(Temp_p(N+2) - Temp_n(N+2))/T_Top < TempTol_hotside
        go_qin_div=0;
    elseif Temp_p(N+2) > Temp_n(N+2)
        if Q_in_p > Q_in_n
            Q_lateral = Q_lateral + 0.5*(Q_in_p-Q_in_n);        % Additional heat input into a single n-type element
        else
            Q_lateral = Q_lateral + 0.001*Q_in_p;
        end
    elseif Temp_n(N+2) > Temp_p(N+2)
        if Q_in_n > Q_in_p
            Q_lateral = Q_lateral - 0.5*(Q_in_n-Q_in_p);     % Additional Heat input into a single p-type element
        else
            Q_lateral = Q_lateral - 0.001*Q_in_n;
        end
    end
    
%     Q_lateral
%      Temp_n(N+2)
%      Temp_p(N+2)
%      Q_in_n
%      Q_in_p

end

%     Temp_n(N+2)
%     Temp_p(N+2)
%        Q_in_te_check = ( N_n* Q_in_n  + N_p * Q_in_p )
        
Temp_n(1)     = x_n(N+2);
Temp_n(2:N+2) = x_n(1:N+1);
Temp_n(  N+3) = T_Bot;

Temp_p(1)     = x_p(N+2);
Temp_p(2:N+2) = x_p(1:N+1);
Temp_p(  N+3) = T_Bot;

%% to obtain the final properties
% [Voc_n, S_n, K_n, Re_n] = TEADVinter_module(DB_n, Temp_n, A_n, t_n);
% [Voc_p, S_p, K_p, Re_p] = TEADVinter_module(DB_p, Temp_p, A_p, t_p);

Re_int_tot = N_n*( sum(Re_n) + 2*R_c_n ) + N_p*( sum(Re_p)+ 2*R_c_p);
if (LOAD_MATCH)
    R_L = Re_int_tot;       
end

Re_tot = Re_int_tot + R_L;
Voc_tot = -N_n*Voc_n + N_p*Voc_p;

I = Voc_tot / Re_tot; 

%% Analysis

V_out = I * R_L;
P_out = V_out*I;


% Q_in_total2 = SolarInput - Total_Area*sigma*eps1*(Temp_n(1)^4 - T_amb^4) - h_conv*Total_Area*(T_Top-T_amb)  
% Q_in_total3 = SolarInput - Total_Area*sigma*eps1*(Temp_p(1)^4 - T_amb^4) - h_conv*Total_Area*(T_Top-T_amb)  

Q_in2 = ( N_n* Q_in_n  + N_p * Q_in_p );
Q_rad = Total_Area*sigma*eps1*(T_Top^4 - T_amb^4);
%Q_rad2 = Total_Area*sigma*eps1*(Temp_n(1)^4 - T_amb^4)
%Q_rad3 = Total_Area*sigma*eps1*(Temp_p(1)^4 - T_amb^4)
Q_conv = h_conv*Total_Area*(T_Top-T_amb);
Eff_TE   = P_out / Q_in2 * 100;
%Eff_sys   = P_out / SolarInput * 100;

%% Output Plots

% fprintf(' Resulted current (I) is : %f Amp\n', I)
% fprintf(' Power output (P_out) is : %f Watts\n', P_out)
% fprintf(' with %6.3f%% efficiency\n\n', Eff)
% fprintf('----------------Simulation Completed!----------------\n\n')

end
