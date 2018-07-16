function fun = TEADVfun_mode4(x, I, S, K, Re, R_c, Total_Area, SolarInput, eps1, h_conv, T_Bot, Psi_h, Psi_c, N, A, Q_lateral, k_f, FF, L)
%TEADVfun Function to solve nodal equations of TEADV for Mode 2
%   Detailed explanation goes here

%   Function in which
%   x(1:N+2) = T_1 to T_N+2

%% Variable Initialization

fun = zeros(N+2,1);

sigma = 5.670e-8;    % Stefan-Boltzmann constant [W m^-2 K^-4]

%% Hot-Side Reservoir
Q_in_total = SolarInput - Total_Area*sigma*eps1*(x(N+2)^4 - T_Bot^4) - h_conv*Total_Area*(x(N+2)-T_Bot);        % subtract radiation loss and convection loss (if any) from the absorber surface

fun(N+2,1) = Psi_h * Total_Area * (x(N+2) - x(1)) - Q_in_total;

%% Hot side of TE
Q_in = A/FF*Psi_h*(x(N+2)-x(1))-(1-FF)/FF*A/L*k_f*(x(1)-x(N+1))+ Q_lateral;

Qk1 = (x( 1 ) - x(2)) * K(1);
Qj = 0.5 * (I^2) * Re(1) + (I^2) * R_c;
Qp =  S(1) * x(1) * I;

fun(1,1) = Qj - Qp + Q_in - Qk1;

%% Node-by-Node Temperature

for i=2 : N
    Qj = 0.5 * (I^2) * (Re(i-1) + Re(i));
    Qp = (S(i) - S(i-1)) * x(i) * I;

    fun(i,1) = Qj - Qp + (x(i-1) - x(i)) * K(i-1) - (x(i) - x(i+1)) * K(i);
end

%% Cold-Side Reservoir

Qj = 0.5 * (I^2) * Re(N) + (I^2) * R_c;
Qp = -S(N) * x(N) * I;

fun(N+1,1) = Qj - Qp - Psi_c * A * (x(N+1) - T_Bot) + (x(N) - x(N+1)) * K(N);

end

