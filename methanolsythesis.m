%the following code simulate the operations of an methanol synthesis process, since the sythesis process is powered by the grid
%the electricity average price can be selected or imported. To import use the third column of X ECONOMIC ANALISIS AND PROJECT importing it as a column vector 
%name the imported data as el_dailycost
m_meOH=5000; %ton/y size methanol sinthesis process
EXM_MEOH=4000; % ton/y example power plant of reference 
EXCOST=1230;%keur CAPEX methanol synthesis reactor
m_co2_reference=688;%kg/h co2 that enter in the system in the reference plant data used to find the cost of the co2 plant 
co2_reference_system_cost=965971.058;%eur capex for ccs in reference plant
cf_reactor=0.9;%capacity factor of the methanol system process
%% find the streams of CO2 and H2 equivalent";
"evaluation of flow rate required of CO_2 and H2 using the first chemical reaction";
mol_meOH=m_meOH*1000000/32.04;%mol
mol_H2=mol_meOH*3; %mol
m_H2=mol_H2*2.016/1000000; %ton/year
m_CO2=mol_meOH*44.01/1000000; %ton/year
"cost of co2 capture using ammine in a distillery from IRENA range 12-22$/ton";
CO2_cost=20*0.92; %eur/ton assuming an ammine process in worst case, with transf $->eur
H2_cost=4000; % eur/ton LCOH deriving from  Quang work
meOH_cost_aprox=CO2_cost*1.375+0.1875*H2_cost; %approximation of the cost of methanol
%% compressors ;
"evaluation of power requirements of methanol synthesis (compression of CO2 and H2) " + ...
"most critical p_ambient to p_req ";
% Parameters 
p1=1.01325e5;% starting prex (Pa, atomosph.)
p1H2=30*10^5;%starting prex hydrogen given from PEM(Pa)
p2O2=150*10^5;%final prex of oxygen Pa (150 baar)
p2=65e5;% final prex(Pa, 65 bar)
p2CO2=20e5;%final prex(Pa, 20 bar)
T_amb=298;% Temperature (K, ambient)
R_CO2=188.9;% specific constant for CO2 (J/kg·K)
R_H2=4124;% specific constant for H2 (J/kg·K)
R_O2=259.81 ;% specific constant for O2 (J/kg·K)
n= 1.3;% politropic index
m_CO2comprex=m_CO2*1000/8760; % Flux in kg/h
m_H2comprex=m_H2*1000/8760;% Flux in kg/h
m_O2comprex=m_H2comprex*8;% Flux in kg/h ,molar ratio of electrolisis is 8
eta=0.75;%efficency
C0_H2=36.858;%base parametere for hydrogen
C0_CO2=2.651;%base parametere for co2
C0_O2=2.327;%base parametere for Oxygen
betaH2=p2/p1H2;%compressio ratio H2 
betaCO2=p2CO2/p1;%compressio ratio CO2
betaO2=5;
%% CO2 COMPREX;
mco2_dot =m_CO2comprex/3600; % from kg/h to kg/s
%power of comprex co2 KW
P1_CO2comp=((mco2_dot*R_CO2*T_amb*n/(eta*(n- 1)))*((p2CO2/p1)^((n- 1)/n)-1))/1000;%kW
P3_CO2comp=P1_CO2comp/0.93;%mechanical eff
P_CO2comp=P3_CO2comp/0.97;%el. eff
T_2co2=T_amb*((30*10^5)/p1)^((n-1)/n);%temperature of co2 exit of compressor
fprintf(' the power required from the CO2 comprex is %.2f kW\n', P_CO2comp);

%% CO2 PUMP
eta_pump=0.85;%assumed pump efficency
DP=45e05;%delta prex required by the pump Pa
R_liqCO2=1.73;% kg/m³
rho_CO2=(p2CO2/(R_CO2*T_amb));%density of co2 at 20 bar
m_CO2pump=mco2_dot/rho_CO2;%flux m^3/s
PpCO2=(DP*m_CO2pump)/(1000*eta_pump);%power KW of CO2 pump
Cost_pump_CO2=(1.417*10^6*(PpCO2/1000)+0.09*10^6)/1000;%keur
fprintf(' the power required from theCO2 pump is %.2f kW\n',PpCO2); 


%% O2 COMPREX
m_O2_dot =m_O2comprex/3600;%kg/s
P1_O2comp=((m_O2_dot*R_O2*T_amb*n/(eta*(n- 1)))*((p2O2/p1H2)^((n- 1)/n)-1))/1000;%kW
P2_O2comp=P1_O2comp;
P3_O2comp=P2_O2comp/0.93;%mechanical eff
P_O2comp=P3_O2comp/0.97;%el. eff
fprintf(' the power required from the O2 comprex is %.2f kW\n', P_O2comp);


%% H2 COMPREX 
mH2_dot =m_H2comprex/3600; % from ton/y to kg/s
%%gasses condition at starting of the synthesis process
p1H2=30*10^5;%starting prex hydrogen given from PEM(Pa)
p2H2=65*10^5;%final prex hydrogen
T_1H2=353.15;%initial temperature hydrogen after PEM (k)
eta=0.75;%efficency
cp_H2=14.570;% kJ/kg*k at starting condition from refprop
%%compression 
rc=1.5;%compression ratio
n=1.3;%Specific Heat Ratio 
eta_poly=0.8;%politropic efficency 
%%1 stage 
T_1H2;%enetering temperature 
p_1=p1H2;%entering prssure in first stage is equal to the mix pressure
p_2=p1H2*rc;%final pressure of stage 
T_2=T_1H2*(p_2/p_1)^((n-1)/n);%must be under 403 K
P_1=((mH2_dot*R_H2*T_2*n/(eta_poly*(n- 1)))*((p_2/p_1)^((n- 1)/n)-1))/1000;%kW
%1 intercooler
T_2a=350;%K asumed temperature (40)
Q_cool1=cp_H2*mH2_dot*(T_2-T_2a);
%2 stage 
T_2a;%enetering temperature 
p_2;%entering prssure in first stage is equal to the mix pressure
p_3=p2H2;%final pressure of stage 
T_3=T_2a*(p_3/p_2)^((n-1)/n);%must be under 403 K
P_2=((mH2_dot*R_H2*T_2a*n/(eta_poly*(n- 1)))*((p_3/p_2)^((n- 1)/n)-1))/1000;%kW
%total power hydrogen compressor
P1_H2comp=P_1+P_2 ; %H2 comprex power KW
P2_H2comp=P1_H2comp/0.93;%mechanical eff
P_H2comp=P2_H2comp/0.97;%el. eff
fprintf(' the power required from the H2 comprex is %.2f kW\n', P_H2comp);
%% compressors costs
Comprex_cost_CO2=C0_CO2*(m_CO2comprex*log(betaCO2))^0.65;%Keur
Comprex_cost_H2=C0_H2*(m_H2comprex*log(betaH2))^0.65;%Keur
Comprex_cost_O2=C0_O2*(m_O2comprex*log(betaO2))^0.65;%Keur

%the total power required for compression taking also the electricity for
%the oxygen compression into account 
P_tot_comp=P_CO2comp+P_H2comp+PpCO2+P_O2comp;
fprintf(' the power required from the total compression system is %.2f kW\n',P_tot_comp);
Tot_cost_comprex_system=(Comprex_cost_H2+Comprex_cost_CO2+Cost_pump_CO2+Comprex_cost_O2);
fprintf('cost of the compression system is %.2f keur\n',Tot_cost_comprex_system);
comprex_system_cost=(Tot_cost_comprex_system)*1000*1.2;%eur + 20% to take into account possible costs increseing
%%
%% cost of the methanol synthesis reactor
"cost of the methanol synthesis process";
Cplant=EXCOST*(m_meOH/EXM_MEOH)^0.65;%keur
fprintf('cost of the methanol synthesis reactor is %.2f keur\n',Cplant);
methanol_reactor_cost=Cplant*1000;%eur
%%
%% HYDROGEN SYSTEM COSTS
PEM_cons=55.8;%KWh/kg
Power_PEM=PEM_cons*m_H2comprex;%KW
specific_cost_PEM=800;%eur/KW average value from papers
cost_PEM1=(specific_cost_PEM*Power_PEM)*1.1;%eur additional 10% for additional costs
specific_cost_PEM2032=113;%eur/KW assuming 60000h (8Y) of PEM life
specific_cost_PEM2040=68;%%eur/KW assuming 60000h (8Y) of PEM life
cost_PEM2=specific_cost_PEM2032*Power_PEM;%investment cost made for replace stack in 2032
cost_PEM3=specific_cost_PEM2040*Power_PEM;%investment cost made for replace stack in 2040
cost_PEM=cost_PEM1+cost_PEM2+cost_PEM3;% eur total cost for the hydrogen generation system
fprintf('cost of the hydrogen production system is %.2f keur\n',(cost_PEM/1000));
%%
%% co2 plant costs
cost_co2_plant=(co2_reference_system_cost*(m_CO2comprex/m_co2_reference)^0.65);%eur
CO2_calculated_cost=(cost_co2_plant+0.03+cost_co2_plant)/(m_CO2*20);
%%
%% rapresentation of fixed costs;
data=[comprex_system_cost methanol_reactor_cost cost_PEM cost_co2_plant];
percentuali = data / sum(data)*100;
labels = {'compression system cost', 'methanol reactor cost', 'hydrogen system cost','carbon dioxide cost' };
pie(data);
legend('compression system cost', 'methanol reactor cost', 'hydrogen system','carbon dioxide', 'Location','southeastoutside');
title('fixed cost analysis');
%%
%% cost of the methanol
%comprex_system_ycost=comprex_system_cost/20; 
%methanol_reactor_ycost=methanol_reactor_cost/20;
%cost_h2pertonmeoh=cost_PEM/m_meOH;%hydrogen cost per ton of methanol
%cost_co2pertonmeoh=cost_co2_plant/m_meOH;%CO2 COST PER ton of methanol
%system_cost=(methanol_reactor_ycost+comprex_system_ycost)/m_meOH;%impact of the system in the cost of methanol
% bar graphic
%data=[cost_h2pertonmeoh;cost_co2pertonmeoh;system_cost];
%figure;
%bar(1, data, 'stacked');  % Cration of the bar "data"
%colormap(jet(length(data)));  % cration of different colours for each segment
%title('costs on the methanol production for ton ');
%ylabel('eur/ton');
%legend({'hydrogen cost', 'carbon dioxide cost', 'sytem cost'},'Location','southeast');
%grid on
%%
%% OPEX
%electricity_cost=mean(el_dailycost);%eur/MWh
electricity_cost=121.9;%eur/MWh
%electricity_cost=90.95;%eur/mwh third simulation value, intermidiate between 60 and average
%electricity_cost=60;% historical data 4th simulation
el_cost_pem=(electricity_cost/1000)*(PEM_cons)*(m_H2*1000);% eur/y cost of electricity for the production of hydrogen
fixedOM=891000;%fixed o&m costs from literature eur/y
el_cost_comprex=P_tot_comp*(8700/1000)*electricity_cost;
%rapresentation of variable costs;
figure
data1=[el_cost_pem fixedOM el_cost_comprex];%variable costs 
percentuali= data / sum(data1)*100;
%labels = {'PEM', 'fixed O&M', 'compressors'};
pie(data1);
legend('PEM', 'fixed O&M', 'compressors');
title('variable cost analysis');


%% 
t=20;%plant lifetime
t1=1:10;
t2=11:20;
%% OXYGEN & METHANOL REVENUE
oxygen_price=150; %eur/ton
oxygen_production=m_O2comprex*8760/1000;%ton/y
REVENUE_oxygenY=oxygen_production*oxygen_price;%eur/y
REVENUE_oxygen_permeohproduction=REVENUE_oxygenY/m_meOH;%eur/tonMeOH
REVENUE_MeOHy=900*m_meOH;%eur/y
CREDIT_CO2=85*m_CO2;%eur/y

%% INITIAL INVESTMENT COST
I0=(comprex_system_cost+methanol_reactor_cost+cost_PEM+cost_co2_plant)/10^6;% total cost Meur
fprintf('the total cost of the system is %.2f Meur\n',I0)
cash_flowy=(REVENUE_MeOHy+REVENUE_oxygenY-el_cost_pem-fixedOM+CREDIT_CO2-el_cost_comprex);% yaar cash flows as income
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

    % evaluation of the actualised value of the cash flow for year t
    discounted_cash_flow = cash_flowy / (1 + k)^t;
    
    % Add the actualised cash flow to the comulative 
    cumulative_discounted_cash_flow = cumulative_discounted_cash_flow + discounted_cash_flow;
    
    % verifiy if the initial investmen has been repaid
    if cumulative_discounted_cash_flow >= TOTAL_INVESTMENT
        % evaluate the fraction of the year where the investment is
        % fullfilled
        DPP = t - 1 + (TOTAL_INVESTMENT - (cumulative_discounted_cash_flow - discounted_cash_flow)) / discounted_cash_flow;
        break;
    end
end

% result visualisation 
if DPP > 0
    fprintf('the Discounted Payback Period (DPP) is  %.2f anni.\n', DPP);
else
    fprintf('the investment can not be payed in the terms .\n');
end
%%
%LEVELIZED COST OF METHANOL
tlife=1:20;
CREDIT_CO2=85*m_CO2;%eur/y
CC=CAPEX;
OC=el_cost_pem+fixedOM+el_cost_comprex;
RO=REVENUE_oxygenY;
E=CREDIT_CO2;
LCOM= (sum((CC+OC-RO-E)./(1+k).^t1)+sum((OC-RO-E)./(1+k).^t2))/(sum(m_meOH./((1+k).^tlife)));
fprintf('The LCOM is: %.2f eur/ton\n', LCOM);

%%comparison between vfixed and variable cost
FIXED_COST=comprex_system_cost+methanol_reactor_cost+cost_PEM+cost_co2_plant;%sum of the fixed costs ù
VARIABLE_COST=el_cost_pem+fixedOM+el_cost_comprex;%sum of the variable costs
figure
data3=[CAPEX VARIABLE_COST];%variable costs 
percentuali = data3 / sum(data1)*100;
%labels = {'FIXED COSTS','VARIABLE COSTS'};
pie(data3);
legend('FIXED COSTS','VARIABLE COSTS');
title('FIXED & VARIABLE COSTS');
% Histogram for comparing costs
figure
bar(data);
set(gca, 'XTickLabel', {'compression system cost', 'methanol reactor cost', 'hydrogen system cost','carbon dioxide cost'}, 'XTickLabelRotation',45);
ylabel('Cost (€)');
title('Cost Comparison');
grid on;
%plot of the cashflows
ARRAYCASH_FLOW1=[NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW1 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 NEW_CASH_FLOW2 ];
cumulative_cash_flow = cumsum(ARRAYCASH_FLOW1);
years=1:20;
figure;
plot(years, cumulative_cash_flow, '-o');
grid on;
xlabel('year'); 
ylabel('comulative cash flow'); % Etichetta sull'asse y

%text(5,40000000,'the electricity price is fixed at 70 eur/MWh','FontSize',8, 'Color', 'black');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROUND TRIP EFFICIENCY CALCULATION 
% energy efficency Methanol sythesis pp
m_MEOHs=(m_meOH*1000)/(8760*60*60);%ammonia mass flow rate kg/s
m_H2s=(m_H2*1000)/(8760*60*60); %hydrogen flow rate kg/s
LHV_MEOH=19700;%kJ/kg ammonia lower heating value 
LHV_H2=120000;%Kj/kg
eta_sythesis=(m_MEOHs*LHV_MEOH)/((m_H2s*LHV_H2)+P_CO2comp+P_H2comp+PpCO2);%efficency of the ammonia sythesis process
%energy requrement AMMONIA SYTEHSIS
EN_SYTinput=((m_H2s*LHV_H2)+P_CO2comp+P_H2comp+PpCO2)*3.6;
EN_SYToutput=(m_MEOHs*LHV_MEOH)*3.6;
EN_SYTlosses=EN_SYTinput-EN_SYToutput;


%%
%energy efficency PEM pp
eta_PEM=(m_H2s*LHV_H2)/(Power_PEM+P_O2comp);
%energy requrement PEM 
EN_PEM_input=(Power_PEM+P_O2comp)*3.6;
EN_PEMoutput=EN_PEM_input*eta_PEM;

EN_PEMlosses=EN_PEM_input-EN_PEMoutput;


%%
%energy efficency ASU
spec_cons_co2=3000;%KJ/kg 
m_CO2s=mco2_dot;
eta_CO2=((m_MEOHs*LHV_MEOH)-(m_CO2s*spec_cons_co2))/(m_MEOHs*LHV_MEOH);
%energy requrement ASU
EN_ASUinput=(m_CO2s*spec_cons_co2)*3.6;
EN_ASUoutput=EN_ASUinput*eta_CO2;
EN_ASUlosses=EN_ASUinput-EN_ASUoutput;


%%
%round efficency 
eta_tot=eta_sythesis*eta_PEM*eta_CO2;




   

