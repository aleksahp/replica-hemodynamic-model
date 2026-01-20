clear; clc;

%% Heart Parameters Values [from Comunale et al. 2021]

% Maximum Contraction Force [mmHg/mL]
pars.Emax_LV = 2.8;
pars.Emax_RV = 0.45; 
pars.Emax_LA = 0.13;
pars.Emax_RA = 0.09;

% Ventricular Stiffness [mmHg/mL]
pars.Emin_LV = 0.07;
pars.Emin_RV = 0.035;
pars.Emin_LA = 0.09;
pars.Emin_RA = 0.045;

% Unstressed chamber volume at zero pressure [mL]
pars.Vp0_LV = 20;
pars.Vp0_RV = 30;
pars.Vp0_LA = 3;
pars.Vp0_RA = 7;

% Duration of One Cardiac Cycle [sec]
pars.T = 0.8;

% Time of Heartbeat in Contractile Phase [sec]
pars.tau1_LV = 0.269 * pars.T;
pars.tau1_RV = 0.269 * pars.T;
pars.tau1_LA = 0.110 * pars.T;
pars.tau1_RA = 0.110 * pars.T;

% Time of Heartbeat in Relaxation Phase [sec]
pars.tau2_LV = 0.452 * pars.T;
pars.tau2_RV = 0.452 * pars.T;
pars.tau2_LA = 0.180 * pars.T;
pars.tau2_RA = 0.180 * pars.T;

% Coefficient for Elastance Steepness for Relaxation Phase [-]
pars.m1_LV = 1.32;
pars.m1_RV = 1.32;
pars.m1_LA = 1.99;
pars.m1_RA = 1.99;

% Coefficient for Elastance Steepness for Contractile Phase [-]
pars.m2_LV = 21.9;
pars.m2_RV = 21.9;
pars.m2_LA = 11.2;
pars.m2_RA = 11.2;

% Time-shift for atrial contraction [sec]
pars.tonset_LV = 0;
pars.tonset_RV = 0;
pars.tonset_LA = 0.85 * pars.T;
pars.tonset_RA = 0.85 * pars.T;

% Valve Resistances [mmHg * s./mL]
pars.Rvalve_LV = 0.01; % Aortic.
pars.Rvalve_RV = 0.01; % Pulmonary.
pars.Rvalve_LA = 0.005; % Mitral.
pars.Rvalve_RA = 0.005; % Tricuspid.

%% Model Parameters Values

% Resistances [mmHg * s./mL].
pars.Rsa = 0.0448; % Systemic arteries.
pars.Rsvb = 0.824; % Systemic vascular bed.
pars.Rsv = 0.0269; % Systemic veins.
pars.Rpa = 0.003; % Pulmonary arteries.
pars.Rpvb = 0.0552; % Pulmonary vascular bed.
pars.Rpv = 0.0018; % Pulmonary veins.

% Compliances. [mL./mmHg]
pars.Csa = 0.983; % Systemic arteries.
pars.Csv = 29.499; % Systemic veins.
pars.Cpa = 6.7; % Pulmonary arteries.
pars.Cpv = 15.8; % Pulmonary veins.

%% Initial Values Parameters

% Initial Pressures. [mmHg]
pars.Psa = 100; % Systemic arteries.
pars.Psv = 4; % Systemic veins.
pars.Ppa = 15; % Pulmonary arteries.
pars.Ppv = 8; % Pulmonary veins.

% Initial chamber volume. [mL]
pars.V0_LV = 149.6;
pars.V0_RV = 189.2;
pars.V0_LA = 71;
pars.V0_RA = 67;

%% Calculate scaling factor [k] for each compartment to ensure that max of E(t) = Emax parameter.

% Temporary Cycle Values.
Nt_k = 5000; % 5000 time steps.
t_k = linspace(0, pars.T, Nt_k); % Between 0 sec and one period.

% LV scaling factor.
nm = max(t_k - pars.tonset_LV,0); % Clamp the numerator to be above 0 to avoid complex numbers.
g1_LV = (nm/pars.tau1_LV).^pars.m1_LV;
g2_LV = (nm/pars.tau2_LV).^pars.m2_LV;
g_max_LV = max( (g1_LV./(1+g1_LV)) .* (1./(1+g2_LV)) );
pars.k_LV = (pars.Emax_LV - pars.Emin_LV) / g_max_LV; % Add scaling factor to parameters.

% RV scaling factor.
nm = max(t_k - pars.tonset_RV,0); % Clamp the numerator to be above 0 to avoid complex numbers.
g1_RV = (nm/pars.tau1_RV).^pars.m1_RV;
g2_RV = (nm/pars.tau2_RV).^pars.m2_RV;
g_max_RV = max( (g1_RV./(1+g1_RV)) .* (1./(1+g2_RV)) );
pars.k_RV = (pars.Emax_RV - pars.Emin_RV) / g_max_RV; % Add scaling factor to parameters.

% LA scaling factor.
nm = max(t_k - pars.tonset_LA,0); % Clamp the numerator to be above 0 to avoid complex numbers.
g1_LA = (nm/pars.tau1_LA).^pars.m1_LA;
g2_LA = (nm/pars.tau2_LA).^pars.m2_LA;
g_max_LA = max( (g1_LA./(1+g1_LA)) .* (1./(1+g2_LA)) );
pars.k_LA = (pars.Emax_LA - pars.Emin_LA) / g_max_LA; % Add scaling factor to parameters.

% RA scaling factor.
nm = max(t_k - pars.tonset_RA,0); % Clamp the numerator to be above 0 to avoid complex numbers.
g1_RA = (nm/pars.tau1_RA).^pars.m1_RA;
g2_RA = (nm/pars.tau2_RA).^pars.m2_RA;
g_max_RA = max( (g1_RA./(1+g1_RA)) .* (1./(1+g2_RA)) );
pars.k_RA = (pars.Emax_RA - pars.Emin_RA) / g_max_RA; % Add scaling factor to parameters.

%% Hemodynamic Model

% Cycle Values
Ncycles = 30; % Number of cycles.
dt = 0.001; % Time step.
t = 0:dt:(Ncycles*pars.T); % Time array.
Nt = length(t);

% Solve the ODE system. Function code and system of equations at bottom.
init = [pars.V0_LV; pars.V0_RV; pars.V0_LA; pars.V0_RA; pars.Psa; pars.Psv; pars.Ppa; pars.Ppv];
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
sol = ode15s(@model, [t(1) t(end)], init, opts, pars);

% Evaluate the solution of the ODE system so they can be plotted.
Y = deval(sol,t);

% Extract chamber volumes.
V_LV_arr = Y(1,:);
V_RV_arr = Y(2,:);
V_LA_arr = Y(3,:);
V_RA_arr = Y(4,:);

% Extract vessel pressures.
Psa_arr = Y(5,:);
Psv_arr = Y(6,:);
Ppa_arr = Y(7,:);
Ppv_arr = Y(8,:);

%% Recalculate pressures to plot.

tmod = mod(t, pars.T); % Keep time within the range of a period.

% LV Elastance & Pressure.
nm_LV = max(tmod - pars.tonset_LV,0); % Clamp the numerator to be above 0 to avoid complex numbers.
g1_LV_arr = (nm_LV/pars.tau1_LV).^pars.m1_LV;
g2_LV_arr = (nm_LV/pars.tau2_LV).^pars.m2_LV;

E_LV_arr = pars.k_LV .* ((g1_LV_arr) ./ (1 + g1_LV_arr)) .* (1 ./ (1 + g2_LV_arr)) + pars.Emin_LV; % Elastance of left ventricle.
P_LV_arr = E_LV_arr .* (V_LV_arr - pars.Vp0_LV); % Pressure of left ventricle.

% RV Elastance & Pressure
nm_RV = max(tmod - pars.tonset_RV,0); % Clamp the numerator to be above 0 to avoid complex numbers.
g1_RV_arr = (nm_RV/pars.tau1_RV).^pars.m1_RV;
g2_RV_arr = (nm_RV/pars.tau2_RV).^pars.m2_RV;

E_RV_arr = pars.k_RV .* ((g1_RV_arr) ./ (1 + g1_RV_arr)) .* (1 ./ (1 + g2_RV_arr)) + pars.Emin_RV; % Elastance of right ventricle.
P_RV_arr = E_RV_arr .* (V_RV_arr - pars.Vp0_RV); % Pressure of right ventricle.

%% Plot PV Loops.

cycle_num = Ncycles; % Focus on last cycle.
t_start = (cycle_num-1) * pars.T;
t_end = cycle_num * pars.T;
idx = find(t >= t_start & t <= t_end);

% Extract values.
V_LV_cycle = V_LV_arr(idx);  % LV volume.
P_LV_cycle = P_LV_arr(idx);  % LV pressure.
V_RV_cycle = V_RV_arr(idx);  % RV volume.
P_RV_cycle = P_RV_arr(idx);  % RV pressure.

% LV PV Loop.
figure(1);
clf
hold on

% Plot main PV loop.
h1 = plot(V_LV_cycle, P_LV_cycle, 'r', 'LineWidth', 2);

% Labelling
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
title('LV PV Loop')

set(gca,'FontSize',20)
grid on
legend(h1, 'LV', 'Location', 'Northwest')

% Axis limits match to article.
xlim([40 160])
ylim([0 140])


% RV PV Loop.
figure(2);
clf
hold on

% Plot main PV loop
h2 = plot(V_RV_cycle, P_RV_cycle, 'b', 'LineWidth', 2);

% Labelling
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
title('RV PV Loop')

set(gca,'FontSize',20)
grid on
legend(h2, 'RV', 'Location', 'Northwest')

% Axis limits match to article
xlim([50 300])
ylim([2 22])

%% Figures for metrics analysis.

% LV Volume [mL] vs. time [sec].
figure(3);
clf
hold on
h3 = plot(t,V_LV_arr, 'r', 'LineWidth', 2);
xlabel('time (sec)');
ylabel('LV Volume (mL)');
title('LV Volume [mL] vs. time [sec]')

% RV Volume [mL] vs. time [sec].
figure(4);
clf
hold on
h4 = plot(t,V_RV_arr, 'b', 'LineWidth', 2);
xlabel('time (sec)');
ylabel('RV Volume (mL)');
title('RV Volume [mL] vs. time [sec]')

% LV Pressure [mmHg] vs. time [sec].
figure(5);
clf
hold on
h5 = plot(t,P_LV_arr, 'r', 'LineWidth', 2);
xlabel('time (sec)');
ylabel('LV Pressure (mmHg)')
title('LV Pressure [mmHg] vs. time [sec]')

% RV Pressure [mmHg] vs. time [sec].
figure(6);
clf
hold on
h6 = plot(t,P_RV_arr, 'b', 'LineWidth', 2);
xlabel('time (sec)');
ylabel('RV Pressure (mmHg)')
title('RV Pressure [mmHg] vs. time [sec]')

% LV Elastance [mmHg/mL] vs. time [sec].
figure(7);
clf
hold on
h7 = plot(t,E_LV_arr, 'r', 'LineWidth', 2);
xlabel('time (sec)');
ylabel('LV Elastance (mmHg/mL)')
title('LV Elastance [mmHg/mL] vs. time [sec]')

% RV Elastance [mmHg/mL] vs. time [sec].
figure(8);
clf
hold on
h8 = plot(t,E_RV_arr, 'b', 'LineWidth', 2);
xlabel('time (sec)');
ylabel('RV Elastance (mmHg/mL)')
title('RV Elastance [mmHg/mL] vs. time [sec]')


%% Model Function Code for ODE15s.

function dxdt = model(t, x, pars)

% Unpack states from init (now x).
V_LV = x(1);
V_RV = x(2);
V_LA = x(3);
V_RA = x(4);

Psa = x(5);
Psv = x(6);
Ppa = x(7);
Ppv = x(8);

tmod = mod(t, pars.T); % Calculate tmod for elastance equations.

% LV Elastance & Pressure. g1 = contractile phase, g2 = relaxation phase.
nm_LV = max(tmod - pars.tonset_LV,0); % Clamp the numerator to be above 0 to avoid complex numbers.
g1_LV = (nm_LV/pars.tau1_LV).^pars.m1_LV;
g2_LV = (nm_LV/pars.tau2_LV).^pars.m2_LV;

E_LV = pars.k_LV * ((g1_LV) / (1 + g1_LV)) * (1 / (1 + g2_LV)) + pars.Emin_LV; % Elastance of left ventricle.
P_LV = E_LV * (V_LV - pars.Vp0_LV); % Pressure of left ventricle.


% RV Elastance & Pressure. g1 = contractile phase, g2 = relaxation phase.
nm_RV = max(tmod - pars.tonset_RV,0); % Clamp the numerator to be above 0 to avoid complex numbers.
g1_RV = (nm_RV/pars.tau1_RV).^pars.m1_RV;
g2_RV = (nm_RV/pars.tau2_RV).^pars.m2_RV;

E_RV = pars.k_RV * ((g1_RV) / (1 + g1_RV)) * (1 / (1 + g2_RV)) + pars.Emin_RV; % Elastance of right ventricle.
P_RV = E_RV * (V_RV - pars.Vp0_RV); % Pressure of right ventricle.


% LA Elastance & Pressure. g1 = contractile phase, g2 = relaxation phase.
nm_LA = max(tmod - pars.tonset_LA,0); % Clamp the numerator to be above 0 to avoid complex numbers.
g1_LA = (nm_LA/pars.tau1_LA).^pars.m1_LA;
g2_LA = (nm_LA/pars.tau2_LA).^pars.m2_LA;

E_LA = pars.k_LA * ((g1_LA) / (1 + g1_LA)) * (1 / (1 + g2_LA)) + pars.Emin_LA; % Elastance of left atria.
P_LA = E_LA * (V_LA - pars.Vp0_LA); % Pressure of left atria.


% RA Elastance & Pressure. g1 = contractile phase, g2 = relaxation phase.
nm_RA = max(tmod - pars.tonset_RA,0); % Clamp the numerator to be above 0 to avoid complex numbers.
g1_RA = (nm_RA/pars.tau1_RA).^pars.m1_RA;
g2_RA = (nm_RA/pars.tau2_RA).^pars.m2_RA;

E_RA = pars.k_RA * ((g1_RA) / (1 + g1_RA)) * (1 / (1 + g2_RA)) + pars.Emin_RA; % Elastance of right atria.
P_RA = E_RA * (V_RA - pars.Vp0_RA); % Pressure of right atria.


%% Valve Flows (Forward)

Q_mi = max(0, (P_LA - P_LV) ./ pars.Rvalve_LA); % Mitral valve.
Q_ao = max(0, (P_LV - Psa) ./ pars.Rvalve_LV); % Aortic valve.
Q_ti = max(0, (P_RA - P_RV) ./ pars.Rvalve_RA); % Tricuspid inflow.
Q_po = max(0, (P_RV - Ppa) ./ pars.Rvalve_RV); % Pulmonary outflow.

%% Systemic / Pulmonary Flows

Q_sa = (Psa - Psv) ./ (pars.Rsa + pars.Rsvb);  % Systemic artery flow.
Q_sv = (Psv - P_RA) ./ pars.Rsv; % Systemic venous flow.
Q_pa = (Ppa - Ppv) ./ (pars.Rpa + pars.Rpvb); % Pulmonary artery flow.
Q_pv = (Ppv - P_LA) ./ pars.Rpv; % Pulmonary venous flow.

%% ODEs

% Change in volume. dV/dt = Qin - Qout.
dV_LV = Q_mi - Q_ao; % Change in left ventricular volume.
dV_RV = Q_ti - Q_po; % Change in right ventricular volume.
dV_LA = Q_pv - Q_mi; % Change in left atrial volume.
dV_RA = Q_sv - Q_ti; % Change in rigth atrial volume.

% Change in pressures. C * dP/dt = Qin - Qout.
dP_sa = (Q_ao - Q_sa) / pars.Csa; % Change in systemic arterial pressure.
dP_sv = (Q_sa - Q_sv) / pars.Csv; % Change in systemic venous pressure.
dP_pa = (Q_po - Q_pa) / pars.Cpa; % Change in pulmonary arterial pressure.
dP_pv = (Q_pa - Q_pv) / pars.Cpv; % Change in pulmonary venous pressure.

dxdt = [dV_LV; dV_RV; dV_LA; dV_RA; dP_sa; dP_sv; dP_pa; dP_pv]; 

end