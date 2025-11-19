%the following code simulate the operations of an ammonia synthesis process, since the sythesis process is powered by the grid
%the electricity average price can be selected or imported. To import use the third column of X ECONOMIC ANALISIS AND PROJECT importing it as a column vector 
%name the imported data as el_dailycost
m_NH3=36500 ; %ton/y size methanol sinthesis process
EXM_MEOH=4000; % ton/y example power plant of reference 
EXCOST=109000;%keur CAPEX methanol synthesis reactor
cf_reactor=0.9;%capa  city factor of the methanol system process
%%find the streams of CO2 and H2 equivalent";
"evaluation of flow rate required of N2 and H2";
mol_NH3=m_NH3*1000000/17.03;%mol/y
mol_N2=mol_NH3/2;%mol/y
mol_H2=(mol_NH3*3)/2; %mol/y
m_H2=mol_H2*2.016/1000000; %ton/year
m_N2=mol_N2*28/1000000; %ton/year
m_N2d=m_N2/365;%ton/day
mh_H2=m_H2*1000/(8760);%;kg/h
mh_N2=m_N2*1000/(8760);%kg/h
ms_N2=mh_N2/3600;
molh_H2=mol_H2/(8760);%;mol/h
molh_N2=mol_N2/(8760);%mol/h
%%gasses condition at starting of the synthesis process
p1H2=30*10^5;%starting prex hydrogen given from PEM(Pa)
p1N2=6*10^5;%starting pressure of nitrogen after air separation
p2O2=150*10^5;%final prex of oxygen Pa (150 baar)
p_final=157;% final prex(Pa, 65 bar)
T_1H2=353.15;%initial temperature hydrogen after PEM (k)
T_1N2=298.15;% Temperature (K, ambient)
V=1;% assumed volume of the mixing chamber m^3 
residence_t=5;%residence time in the mixing chamber of the gas mixture (s)
R=8.314;%J/mol*K cost. of gasses
R_N2=296.4;% J/kg*K
R_O2=259.81 ;% specific constant for O2 (J/kg·K)
T_amb=298.15;%K
betaO2=5;
n= 1.3;% politropic index
eta=0.75;%efficency
cp_H2=14.570;% kJ/kg*k at starting condition from refprop
cp_N2=1.0495;%kJ/kg*K at starting conditionns from refprop
C0_H2=36.858;%base parametere for hydrogen
C0_O2=2.327;%base parametere for Oxygen
rc_N2=2.24;%compression ratio nitrogen
%NITROGEN COMPRESSOR
%1 stage NITROGEN  
T_1N2;%entering temperature 
p1N2;%entering prssure in first stage is equal to the mix pressure
p2N2=p1N2*rc_N2;%final pressure of stage 
T_2N2=T_1N2*(p2N2/p1N2)^((n-1)/n);%must be under 403 K
P_1N=((ms_N2*R_N2*T_1N2*n/(eta*(n- 1)))*((p2N2/p1N2)^((n- 1)/n)-1))/1000;%kW
%intercooler nitrogen 
T2a_N2=T_1N2;
Q_coolN2=cp_N2*ms_N2*(T_2N2-T2a_N2);%KW
%2 stage nitrogen 
T2a_N2;%enetering temperature 
p2N2;%entering prssure in first stage is equal to the mix pressure
p3N2=p1H2;%final pressure of stage 
T_3N2=T2a_N2*(p3N2/p2N2)^((n-1)/n);%must be under 403 K
P_2N=((ms_N2*R_N2*T2a_N2*n/(eta*(n- 1)))*((p3N2/p2N2)^((n- 1)/n)-1))/1000;%kW


%%MIXING PROCESS
T_mix=((mh_H2*cp_H2*T_1H2)+(mh_N2*cp_N2*T_3N2))/(mh_H2*cp_H2+mh_N2*cp_N2);% final mixing temperature(K)
mol_mix=((molh_H2+molh_N2)/(60*60))*residence_t;%total number of moles in the mixing chamber
p_mix=p1H2;%final mixing pressure (bar)
H2_density=2.0196;% density of hydrogen at mixing equilibrum kg/m^3
N2_density=28.352;%density of nitrogen at mixing equilibrum kg/m^3
vol_flow=mh_H2/H2_density+mh_N2/N2_density;%volumetric flow rate after mixing m^3/h
ms_mix=(mh_N2+mh_H2)/3600;%kg/s
m_mix_h=mh_N2+mh_H2;%kg/h
V_n2=mh_N2/N2_density;%m^3/h

%% PRECOOLING
T_2=313;
cp_mix=(molh_H2/(molh_H2+molh_N2))*cp_H2+(molh_N2/(molh_H2+molh_N2))*cp_N2;%kJ/kg*k
Q_cool=cp_mix*ms_mix*(T_mix-T_2);%KW


%% compression 
rc=2.1;%compression ratio
n=1.4;%Specific Heat Ratio 
R_mix=976.36; %J/kg*K
eta=0.8;%politropic efficency 

%1 stage 
T_2;%enetering temperature 
p_2=p_mix;%entering prssure in first stage is equal to the mix pressure
p_3=p_2*rc;%final pressure of stage 
T_3=T_2*(p_3/p_mix)^((n-1)/n);%must be under 403 K
P_1=((ms_mix*R_mix*T_2*n/(eta*(n- 1)))*((p_3/p_2)^((n- 1)/n)-1))/1000;%kW

%1 intercooler
T_3a=313;%K asumed temperature (40)
Q_cool1=cp_mix*ms_mix*(T_3-T_3a);

%2 stage 
T_3a;%enetering temperature 
p_3;%entering prssure in first stage is equal to the mix pressure
p_4=p_3*rc;%final pressure of stage 
T_4=T_3a*(p_4/p_3)^((n-1)/n);%must be under 403 K
P_2=((ms_mix*R_mix*T_3a*n/(eta*(n- 1)))*((p_4/p_3)^((n- 1)/n)-1))/1000;%kW

%2 intercooler
T_4a=313;%K asumed temperature (40)
Q_cool2=cp_mix*ms_mix*(T_4-T_4a);

%3 stage 
T_4a;%enetering temperature 
p_4;%entering prssure in first stage is equal to the mix pressure
p_5=157*10^5;%final pressure of stage 
T_5=T_3a*(p_5/p_4)^((n-1)/n);%must be under 403 K
P_3=((ms_mix*R_mix*T_4a*n/(eta*(n- 1)))*((p_5/p_4)^((n- 1)/n)-1))/1000;%kW
P_fluid=(P_1+P_2+P_3+P_1N+P_2N);
P_shaft=P_fluid/0.93; %total power required for the compression of the mixture from p_mix to p=157 bar using mech. eff
P_comprex_result=P_shaft/0.97;%el eff
fprintf(' the power required from the  comprex is %.2f kW\n', P_shaft);
%compression cost
a=8400;
b=3100;
n_cost=0.6;
Comprex_cost=(a+b*P_comprex_result^n_cost)/1000;%keur value found using empiric function as theCHEMICAL ENGINEERING DESIGN Principles, Practice andEconomics of Plant andProcess Design GAVIN TOWLER RAY SINNOTT
%% OXYGEN COMPRESSOR
m_O2comprex=mh_H2*8;% Flux in kg/h ,molar ratio of electrolisis is 8
m_O2_dot =m_O2comprex/3600;%kg/s
P1_O2comp=((m_O2_dot*R_O2*T_amb*n/(eta*(n- 1)))*((p2O2/p1H2)^((n- 1)/n)-1))/1000;%kW
P2_O2comp=P1_O2comp;
P_O2comp=P2_O2comp/0.93;%using also the mech eff
P_O2elcomp=P_O2comp/0.97;%using the el eff
Comprex_cost_O2=C0_O2*(m_O2comprex*log(betaO2))^0.65;%Keur
%COMPRESSION SYSTEM
P_c_system=P_O2comp+P_shaft;%KW
%% recycling compression 
m_recy=(m_mix_h*0.8)/3600;% flow rate of the mix that must recycolate in the system kg/s
p_recy_start=145;%bar pression oafter the separation unit 
p_recy_end=157;%bar pression at the end of the compression
T_recy=268.13;%K temperature after the separation phase 
P1_recy_compr=((m_recy*R_mix*T_recy*n/(eta*(n- 1)))*((p_recy_end/p_recy_start)^((n- 1)/n)-1))/1000;%kW
P2_recy_compr=P1_recy_compr/0.85;%taking into acount the isoentropic eff
P_recy_compr=P2_recy_compr/0.93;%taking into account the mech. eff.
P_recy_compr_el=P_recy_compr/0.97;%using el eff
Comprex_cost_recy=(a+b*P_recy_compr^n_cost)/1000;%keur
P_c_tot=P_comprex_result+P_recy_compr_el+P_O2elcomp;%power required by the compressors in the ammonia loop kW
E_c_tot=P_c_tot*8760;%annual energy requred by the ammonia compressors kWh

%% FIXED COSTS

%% COMPRESSION SISTEM COST
Comprex_system_cost=Comprex_cost_O2+Comprex_cost+Comprex_cost_recy;%keur

%% AMMONIA REACTOR COST FUNCTION 
AF=2.75;%alloy factor cost
C_fixed=63000;%reference cost of cpntrol system
LM=2.3;% labor menagment cost
CEPCI_19=607.5;%chemical plant cost index (CEPCI 2010 =1000)
C_j=110000;%cost of the of the component j
s_reg=3.4;%pressure regulation factor for reactor that works at 150 bar
s_j=20;%reference size 
n_j=0.52;%scaling factor
Reactor_cost=(C_fixed+LM*(CEPCI_19/1000)*AF*C_j*(s_reg/s_j)^n_j)/1000;%chemical cost function keur

%% AMMONIA LOOP WITOUT COMPRESSIORS AND REACTORS 
C_ref=8360460;%eur 
size_ref=300;%ton/d
other_cost=(C_ref*(100/size_ref)^0.65)/1000;% keur scaling down the cost


%% NITROGEN PRODUCTION COST 
m_Nref=250;%ton/d
ref_N2_cost=6922430;%eur
N2_plant_cost=(ref_N2_cost*(m_N2d/m_Nref)^0.65)/1000;%nitrogen plant cost keur

%% HYDROGEN SYSTEM COSTS
PEM_cons=55.8;%KWh/kg
Power_PEM=PEM_cons*mh_H2;%KW
ref_H2cost=6.3787*10^6;%eur
ref_H2power=6012;%KW
Cost_PEM_scaled=(ref_H2cost*(Power_PEM/ref_H2power)^0.85)/1000;%keur

% EVALUATION OF COST
data = [Reactor_cost, Comprex_system_cost, Cost_PEM_scaled,other_cost, N2_plant_cost,];
percentuali = data / sum(data) * 100;

% Pie chart for costs distribution
figure;
%subplot(1, 2, 1); % Subplot for side-by-side plots
pie(data);
legend('Reactor Cost','Compression System Cost', 'Hydrogen System','ammonia loop', 'Nitrogen', 'Location', 'southeastoutside');
title('fixed Costs Distribution');

% Histogram for comparing costs
%subplot(1, 2, 2);
figure
bar(data);
set(gca, 'XTickLabel', {'Reactor Cost', 'Compression System', 'Hydrogen System','ammonia loop','Nitrogen Cost'}, 'XTickLabelRotation',45);
ylabel('Cost (k€)');
title('Cost Comparison');
grid on;

%%  EVALUATION OF VARIABLE COSTS ARE MAINLY RELATED TO ELECTRICITY 
%OPEX
%electricity_cost=mean(el_dailycost);%eur/MWh
%electricity_cost=90.95;%eur/mwh third simulation value, intermidiate between 60 and average
%electricity_cost=85.05;%eur/MWh
%electricity_cost=30;%eur/MWh
%electricity_cost=85;
electricity_cost=121.9;
%electricity_cost=60;%eur/MWh
%electricity_cost=50;% historical data 4th simulation
el_cost_pem=(electricity_cost/1000)*(PEM_cons)*(m_H2*1000);% eur/y cost of electricity for the production of hydrogen
I0=(Reactor_cost+ Comprex_system_cost+ Cost_PEM_scaled+ N2_plant_cost+other_cost)/10^3;
fixedOM=0.03*I0;%fixed o&m costs from literature eur/y
EL_comprex=(E_c_tot)*(10^-3)*electricity_cost;%eur
totalvariable=el_cost_pem+fixedOM+EL_comprex;


%% 
t=20;%plant lifetime
t1=1:10;
t2=11:20;
%% OXYGEN & METHANOL REVENUE
oxygen_price=150; %eur/ton
oxygen_production=m_O2comprex*8760/1000;%ton/y
REVENUE_oxygenY=oxygen_production*oxygen_price;%eur/y
REVENUE_oxygen_perNH3production=REVENUE_oxygenY/m_NH3;%eur/tonN2
REVENUE_NH3y=600*m_NH3;%eur/y

%% INITIAL INVESTMENT COST
I0=(Reactor_cost+ Comprex_system_cost+ Cost_PEM_scaled+ N2_plant_cost+other_cost)/10^3;% total cost Meur
fprintf('the total cost of the system is %.2f Meur\n',I0)
cash_flowy=(REVENUE_NH3y+REVENUE_oxygenY-el_cost_pem-fixedOM-EL_comprex);% yaar cash flows as income
NPV0=-I0*10^6;  % Start with the initial investment (negative)
%%
%% FIND THE CAPEX WITH A INTEREST OF 8% AND A PAY BACK PERIOD OF 10 YEARS
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
%% DPB 
cumulative_discounted_cash_flow=0;
DPP = 0;
for t = 1:60

    % evaluation of the actualised cash flow for the year t 
    discounted_cash_flow = cash_flowy / (1 + k)^t;
    
    % add the cash flow to the comulative
    cumulative_discounted_cash_flow = cumulative_discounted_cash_flow + discounted_cash_flow;
    
    % verify if the initial investment has been paid
    if cumulative_discounted_cash_flow >= TOTAL_INVESTMENT
        % evaluate when during the year the investment has been paid off
        DPP = t - 1 + (TOTAL_INVESTMENT - (cumulative_discounted_cash_flow - discounted_cash_flow)) / discounted_cash_flow;
        break;
    end
end

% visualisation of the result
if DPP > 0
    fprintf('the Discounted Payback Period (DPP) is  %.2f anni.\n', DPP);
else
    fprintf('the investment can not be payed in the terms .\n');
end
%%
%LEVELIZED COST OF AMMONIA 
tlife=1:20;
CC=CAPEX;
OC=el_cost_pem+fixedOM+EL_comprex;
RO=REVENUE_oxygenY;

LCONH3=(sum((CC+OC-RO)./(1+k).^t1)+sum((OC-RO)./(1+k).^t2))/(sum(m_NH3./((1+k).^tlife)));
fprintf('The LCONH3 is: %.2f eur/ton\n', LCONH3);

% Pie chart for costs distribution
data1=[CAPEX,totalvariable];
figure;
pie(data1);
legend('fixed cost','variable cost');
title('comparison between fixed and variable costs');
% Pie chart for costs distribution
data2=[el_cost_pem,EL_comprex,fixedOM];
figure;
pie(data2);
legend('PEM','comprx','fixed O&M');
title('comparison between variable costs');
%plot of the cashflows
ARRAYCASH_FLOW=[NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 ];
cumulative_cash_flow = cumsum(ARRAYCASH_FLOW);
years=1:20;
figure;
plot(years, cumulative_cash_flow, 'r-o');
grid on;
xlabel('year'); 
ylabel('comulative cash flow'); % Etichetta sull'asse y
%text(5,70000000,'the electricity price is fixed at 60 eur/MWh','FontSize',8, 'Color', 'black');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% energy efficency ammonia sythesis pp
m_NH3s=(m_NH3*1000)/(8760*60*60);%ammonia mass flow rate kg/s
m_H2s=(m_H2*1000)/(8760*60*60); %hydrogen flow rate kg/s
LHV_NH3=18500;%kJ/kg ammonia lower heating value 
LHV_H2=120000;%Kj/kg
m_N2s=(m_N2*1000)/(8760*60*60);%kg/s
eta_sythesis=(m_NH3s*LHV_NH3)/((m_H2s*LHV_H2)+P_comprex_result+P_recy_compr_el);%efficency of the ammonia sythesis process
%energy requrement AMMONIA SYTEHSIS
EN_SYTinput=((m_H2s*LHV_H2)+P_comprex_result+P_recy_compr_el)*3.6;
EN_SYToutput=(m_NH3s*LHV_NH3)*3.6;
EN_SYTlosses=EN_SYTinput-EN_SYToutput;


%%
%energy efficency PEM pp
eta_PEM=(m_H2s*LHV_H2)/(Power_PEM+P_O2elcomp);
%energy requrement PEM 
EN_PEM_input=(Power_PEM+P_O2elcomp)*3.6;
EN_PEMoutput=EN_PEM_input*eta_PEM;

EN_PEMlosses=EN_PEM_input-EN_PEMoutput;
%%
%energy efficency ASU (N2)
spec_cons_N2=1.5156;%MJ/kg 
eta_N2=((m_NH3s*LHV_NH3/1000)-(m_N2s*spec_cons_N2))/(m_NH3s*(LHV_NH3/1000));
%energy requrement ASU
EN_ASUinput=(m_N2s*spec_cons_N2*1000)*3.6;
EN_ASUoutput=EN_ASUinput*eta_N2;
EN_ASUlosses=EN_ASUinput-EN_ASUoutput;
%%
%round efficency    
eta_tot=eta_sythesis*eta_PEM*eta_N2; 


   
