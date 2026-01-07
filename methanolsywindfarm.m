%% Wind Speed Nominal Interval Definition
wind_speed_nominal = 0:0.5:25; % Wind speed bins from 0 to 25 m/s in 0.5 m/s steps
%figure(99)
%plot(wind_speed_nominal, P_out);
%grid on
%xlabel('Wind Speed (m/s)');
%ylabel('Power (kW)');
%title('Power Curve Vestas V117-4.0 MW');

%% Process Real-Time Wind Speed Data
T = readtable('wind speed.xlsx','Range','A:A');
First_column = T{:,1};
WS1=First_column;%
WS_h = (WS1(1:2:end) + WS1(2:2:end)) / 2; % Half-hourly average
WS_mapped = round(WS_h * 2) / 2; % Round to nearest 0.5 m/s
WS_mapped(WS_mapped > max(wind_speed_nominal)) = 0; % Clip to max nominal range

%% Time Series Setup
year_selected = 2022;
start_time = datetime(year_selected, 1, 1, 0, 0, 0);
end_time = datetime(year_selected, 12, 31, 23, 0, 0);
timeseries1 = start_time:hours(1):end_time;

%% Power Generation Calculation
Y= readtable('wind speed.xlsx','Range','F2:F52');
Power_column = Y{:,1};
P_out=Power_column ;%vector power of the wind turbine in fuction of the wind speed
P_generated = interp1(wind_speed_nominal, P_out, WS_mapped, 'linear', 0); % kW
Ur = (P_generated / max(P_generated))*100; % Utilization rate (normalized)

% Efficiency curve application based on utilization rate
L_s = 2.2; % 2.2% minimum threshold
soglia_15 = 15; % 15% threshold
Fy = 2030;
eta = zeros(size(Ur));

idx1 = Ur >= L_s & Ur < soglia_15;
idx2 = Ur >= soglia_15;
idx3=Ur<=L_s;

eta(idx1) = ((0.00005*Ur(idx1).^5 - 0.0061*Ur(idx1).^4 + 0.2372*Ur(idx1).^3 - 4.2014*Ur(idx1).^2 + 36.675*Ur(idx1) - 62.87) * (1 + 0.0025)^(Fy - 2020)) ;
eta(idx2) = ((-0.149*Ur(idx2) + 74.977) * (1 + 0.0025)^(Fy - 2020)) ;
eta(idx3)=0;
P_netta = P_generated ;% Net power after efficiency (kW)

%% Hydrogen Production Using Dynamic Efficiency Only
kWh_per_kg_H2 = 33.3; % PEM conversion factor (kWh/kg)
optimal_turbines = 2.7;
P_farm = P_netta * optimal_turbines; % Total wind farm net power (kW)
E_generated = P_farm / 1e3; % MWh
E_tot = sum(E_generated);
E_tot_cf = sum((P_generated .* optimal_turbines) / 1e3);
eta_perc=eta./100;

P_PEM_kW = max(P_farm); % Sizing PEM system by peak power (kW)
m_hydrogen = P_farm.*eta_perc./ kWh_per_kg_H2; % Hydrogen production (kg/h)
m_H2required = 107.7420; % Required hydrogen flow (kg/h)

%% Hydrogen Storage Reserve Initialization
reserve = zeros(8760, 1);
reserve(1) = 5000; % Initial H2 reserve in kg

for t = 1:length(P_generated)
    surplus = m_hydrogen(t) - m_H2required;
    if t == 1
        reserve(t) = surplus + reserve(t);
    else
        reserve(t) = reserve(t - 1) + surplus;
    end
    reserve1 = reserve(1:1200); % First 6552 hours (Jan to Sept)
    storage_size = max(reserve1); % Max required storage

    reserve(t) = min(reserve(t), storage_size);
    reserve(t) = max(reserve(t), 0);
end

H_total = m_hydrogen + reserve;
reserve_tot = sum(reserve);

fprintf('Optimal Number of Turbines: %d turbines\n', optimal_turbines);
fprintf('Minimum Required Hydrogen Storage Capacity: %.2f kg\n', storage_size);

if any(H_total < m_H2required)
    fprintf('The hydrogen production is not sufficient\n');
else
    fprintf('The hydrogen production is sufficient\n');
end

H_history = sum(H_total);
cf = E_tot_cf / (8760 * optimal_turbines * 4);
fprintf('The capacity factor of the plant is: %.2f \n', cf);

%% Plot Hydrogen Reserve and Production
figure(1)
plot(timeseries1, reserve, 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Storage Level (kg)');
title('Hydrogen Tank Level Over Time');
hold on
plot(timeseries1, H_total, 'g', 'LineWidth', 1.5);
yline(m_H2required, 'r', 'Required Hydrogen');
hold off;

%% PEM Cost Calculation
specific_cost_PEM = 800; % EUR/kW
cost_PEM1 = (specific_cost_PEM * P_PEM_kW) * 1.1;
cost_PEM2 = 113 * P_PEM_kW; % Replacement 2032
cost_PEM3 = 68 * P_PEM_kW;  % Replacement 2040
cost_PEM = cost_PEM1 + cost_PEM2 + cost_PEM3;

%% Wind Farm CAPEX
h_hub = 91.5;
R = 117 / 2;
P_nom = 4000;
c_found = (303.24 * h_hub * (R^2) * pi * 0.4037) / 1000;
c_transp = ((1.581E-5 * P_nom^2 - 0.0375 * P_nom + 54.7) * P_nom) / 1000;
c_instal = (1.965 * (h_hub^2 * R)^1.1736) / 1000;
c_civil = ((2.17E-6 * P_nom^2 - 0.0145 * P_nom + 69.54) * P_nom) / 1000;
c_eng = ((9.94E-4 * P_nom + 20.31) * P_nom) / 1000;
BOS = c_found + c_transp + c_instal + c_civil + c_eng;
TCC = 4200000;
c_el_interf = ((3.49E-6 * P_nom^2 - 0.0221 * P_nom + 109.7) * P_nom) / 1000;
c_el_install = 0.65 * c_el_interf * optimal_turbines;
CAPEX_turb = optimal_turbines * (BOS + TCC) + c_el_install;

%% Wind Farm OPEX
eta_elec = 0.95;
O_and_M = 0.007 * E_tot;
decomissioning_cost = 0.0107 * P_nom * optimal_turbines;
land_fill_cost = 0.00108 * E_tot;
OPEX = ((land_fill_cost + (1 / eta_elec) * (decomissioning_cost + O_and_M)) / E_tot) * 1000;

%% Hydrogen Storage Cost
op_cost = 0.32 * 0.92;
maint_cost = 0.05 * 0.92;
d = 0.1;
t = 40;
LCOS = 1.5 * 0.92;
cushion_perc = 0.5;
reserve = storage_size + storage_size * cushion_perc;
cf = 0.8;
level_TCC = LCOS * reserve + op_cost + maint_cost;
Hydrogen_stroage_cost = ((cf * level_TCC * ((1 + d)^t - 1)) / (d * (1 + d)^t)) / 1000;

%% Compressor Cost (Oxygen)
n = 1.4;
p1H2 = 30e5;
p2O2 = 150e5;
C0_O2 = 2.327;
R_O2 = 259.81;
T_amb = 298.15;
eta_poly=0.8;
m_O2comprex = max(m_hydrogen) * 8;
m_O2_dot = m_O2comprex / 3600;
P1_O2comp = ((m_O2_dot * R_O2 * T_amb * n / (eta_poly * (n - 1))) * ((p2O2 / p1H2)^((n - 1)/n) - 1)) / 1000;
P2_O2comp = P1_O2comp / 0.85;
P3_O2comp = P2_O2comp / 0.93;
P_O2comp = P3_O2comp / 0.97;
betaO2 = 5;
Comprex_cost_O2 = C0_O2 * (m_O2comprex * log(betaO2))^0.65;

%% Cost Summary
m_meOH=5000;%ton/year
mol_meOH=m_meOH*1000000/32.04;%mol
m_CO2=mol_meOH*44.01/1000000;%ton/year
Cost_PEM_scaled = cost_PEM / 1000;%KEUR
CO2_PLANT = 1051.6;%KEUR
Comprex_system_cost = 1747.2 - 256.5627 + Comprex_cost_O2;%KEUR 
MEOH_plant_cost = 1422.0;%KEUR
wind_farm = CAPEX_turb / 1000;%KEUR
data = [Comprex_system_cost, Cost_PEM_scaled, CO2_PLANT, MEOH_plant_cost, wind_farm, Hydrogen_stroage_cost];
percentuali = data / sum(data) * 100;

figure;
pie(data);
legend('Compression System Cost', 'Hydrogen System', 'CO2 plant cost', 'Methanol sythesis plant', 'Wind Farm', 'Hydrogen Storage', 'Location', 'southeastoutside');
title('Fixed Costs Distribution');

%% Revenue & LCOH
oxygen_price = 150;
oxygen_production = m_O2comprex * 8760 / 1000;
REVENUE_oxygenY = oxygen_production * oxygen_price;

r = 0.08;
FCR = (r * (1 + r)^25) / ((1 + r)^25 - 1);
LCOE_wf = (FCR * CAPEX_turb) / E_tot + OPEX;
el_cost_pem = (LCOE_wf/1000)*(mean(eta_perc)/kWh_per_kg_H2)*((sum(m_hydrogen)*8760));% eur/y cost of electricity for the production of hydrogen%LCOE_wf * sum(P_farm) / 1000;
I0 = ( Comprex_system_cost + Cost_PEM_scaled + MEOH_plant_cost + wind_farm + Hydrogen_stroage_cost+CO2_PLANT) / 1000;
fixedOM = 0.03 * (I0-wind_farm) + OPEX;
E_c_tot = 1.5352e+07 - (404.6583 * 8760) + sum(P_O2comp * 8760);
EL_comprex = (E_c_tot * 1e-3) * LCOE_wf;
totalvariable = el_cost_pem + fixedOM + EL_comprex;

NPV0 = -I0 * 1e6;
k = 0.08;
n_payback = 10;
tlife = 1:20;
t1 = 1:10;
t2 = 11:20;
alphafactor = (1 - (1 + k)^(-n_payback)) / k;
CREDIT_CO2=85*m_CO2;%eur/y
CAPEX = (-NPV0) / alphafactor;
CC = CAPEX;
OC = el_cost_pem +fixedOM + EL_comprex;
RO = REVENUE_oxygenY;
m_meOH=5000;%ton/year
E=CREDIT_CO2;
LCOM= (sum((CC+OC-RO)./(1+k).^t1)+sum((OC-RO)./(1+k).^t2))/(sum(m_meOH./((1+k).^tlife)));

figure(6)
plot (timeseries1,eta)
grid on
xlabel('Time')
ylabel('Efficiency')
title('Efficiency of the PEM plant durnig the year')
REVENUE_MEOHy= 1220.1*m_meOH;%eur/y
cash_flowy=(REVENUE_MEOHy+REVENUE_oxygenY-el_cost_pem-fixedOM-EL_comprex);% yaar cash flows as income
NPV0=-I0*10^6;  % Start with the initial investment (negative)
%%
%FIND THE CAPEX WITH A INTEREST OF 8% AND A PAY BACK PERIOD OF 10 YEARS
k=0.08;%interest rate fixed 8%
n_payback=10;%payback of investement
%total capex of the investment
%atualized_cash_flowy=cash_flowy/(1 + k)^t;%atualization of the cash flows
%=sum(atualized_cash_flowy);%
alphafactor=(1-(1+k)^(-10))/k;
CAPEX=((-NPV0)/alphafactor);%total yearly capex  of the investment
TOTAL_INVESTMENT=sum(CAPEX./(1+k).^t1);
NEW_CASH_FLOW1=cash_flowy-CAPEX;
NEW_CASH_FLOW2=cash_flowy;
tot_cash_flow=sum((NEW_CASH_FLOW1)./(1 + k).^t1)+sum((NEW_CASH_FLOW2)./(1 + k).^t2);
NPV=tot_cash_flow;
fprintf('The Net Present Value (NPV) is: %.2f\n', NPV);
PB=TOTAL_INVESTMENT/(cash_flowy);
fprintf('The PB time (PB) is: %.2f\n', PB);
%%DPB 
cumulative_discounted_cash_flow=0;
DPP = 0;
for t = 1:60

    % Calcolo del valore attuale del flusso di cassa per l'anno t
    discounted_cash_flow = cash_flowy / (1 + k)^t;
    
    % Aggiungi il flusso di cassa scontato al cumulativo
    cumulative_discounted_cash_flow = cumulative_discounted_cash_flow + discounted_cash_flow;
    
    % Verifica se l'investimento iniziale è recuperato
    if cumulative_discounted_cash_flow >= TOTAL_INVESTMENT
        % Calcola la frazione dell'anno in cui l'investimento viene recuperato
        DPP = t - 1 + (TOTAL_INVESTMENT - (cumulative_discounted_cash_flow - discounted_cash_flow)) / discounted_cash_flow;
        break;
    end
end

% Visualizzazione del risultato
if DPP > 0
    fprintf('the Discounted Payback Period (DPP) is  %.2f anni.\n', DPP);
else
    fprintf('the investment can not be payed in the terms .\n');
end

