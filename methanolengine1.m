% daily data
%timeseries =[timeseries];  % day in the year
year_selected = 2023;

% timevector
start_time = datetime(year_selected, 1, 1, 0, 0, 0); %  (1 jen alle 00:00)
end_time = datetime(year_selected, 12, 31, 23, 0, 0); % (31 Dic alle 23:00)
timeseries = start_time:hours(1):end_time; % Vettore di timestamp con intervalli di un'ora
prices =[el_dailycost];  % daily prices


results2 = []; % Array to memorize the results

%for percentile_90 =prctile(prices,10:10:90);  % change value between 10 to 90

%evaluate the 90° percentile of prices
percentile_90 = prctile(prices,90);
%percentile_90=155.52;
% arrey that rapresent the engione status 
%motor_on = prices > percentile_90;  % engin is on if the prices is higher than the 90%
for j = 1:length(prices)
    % check if the two consequtive values are higher than the treshold
    if prices(j) >percentile_90  && prices(j+1) > percentile_90
        motor_on(j) = 1; % if it is true activate motor on 
    else 
        motor_on(j)=0;
    end
end

% evaluation of the energy pruduction assuming a daily continuos operation
% of the engine 
energy = 9.73*motor_on;  % Energy = 9.73 when engine is on, 0 when is off

% comulative energy in time 
cumulative_energy = cumsum(energy);
%evaluation of monthly energy production 
months = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];  % Days in each month
% sum of energy in each month
starting_index = 1;
energy_month = zeros(1, 12); % Array monthly energy
for month = 1:12
    %calculate monthly hours
    monthly_hours = days_in_month(month) * 24;

    % evaluate final index for current month
    final_index = starting_index + monthly_hours - 1;
    % Ensure final_index does not exceed the length of the energy array
    final_index = min(final_index, length(energy));
   
    % Sum energy of current month
    energy_month(month) = sum(energy(starting_index:final_index));

    % update indes
    starting_index = final_index + 1;
end

% graphic creation
figure;
bar(energy_month);
xlabel('month');
ylabel('energy production MWh');
title('monthly energy generation');
grid on;

% days in whic the engine is working
figure;
plot(timeseries(motor_on == 1), cumulative_energy(motor_on == 1), 'r');
xlabel('month');
ylabel('energy MWh');
title('COMULATIVE ENERGY GENERATED');
grid on;
%switch on dates
display(timeseries(motor_on==1))
%storage level of the methanol tank
storagecapac=5000; %capacity of storage in tons assumed 
maxload_cons_ton=340.2/1000; % consumption of methanol on ton/MWh 0.791 assumed methanol density;
pilotoilcons=8.45/1000;%consmption of pilot oil ton/MWh 
consumption=energy*maxload_cons_ton;% daily consumption 
pilotconsumptionD=energy*pilotoilcons;%daily oil consumption
storagelvl=storagecapac;

% Initial storage level
storagelvl = storagecapac;
levels = zeros(1, length(timeseries)); % Array for storing level over time
levels(1) = storagelvl; % Starting level

% Loop through each hour and calculate methanol consumption
for i = 1:length(timeseries)
    storagelvl = storagelvl - consumption(i); % Update storage level

    % Ensure storage level does not drop below zero
    if storagelvl < 0
        storagelvl = 0;
    end
    
    % Store the current storage level
    levels(i) = storagelvl;
end

figure;
line(timeseries,levels);
xlabel('month');
ylabel('storage level (ton)');
title('reserve of methanol on time');
set(gca, 'XTickLabel', {'Gen', 'Feb', 'Mar', 'Apr', 'Mag', 'Giu', 'Lug', 'Ago', 'Set', 'Ott', 'Nov', 'Dic'});
grid on;
%engine capacity factor
total_en=cumulative_energy(end);%MWh/y
cf=(total_en/(9.7*8760));
fprintf('the capacity factor of this plant is %.2f \n',cf);

%%
% evaluatinig the economical aspect
t=20;%plant lifetime
t1=1:10;
t2=11:25;
k=0.08;
I0_ENGINE=-6*10^6;%engine cost, eur
OPEX=(0.03*I0_ENGINE)/20;%operation and mantaineance costs, eur
tot_consum=storagecapac-storagelvl;%ton/y 
tot_pilotoil_cons=total_en*pilotoilcons;%total consumption of pilot oil (MWh/y)*(ton/mwh)-> ton/y
methanol_cost=360.81;% eur/ton
pilotoil_cost=900*0.92;%eur/ton (transformation from USD to eur)
rowenergy=energy';
revenue_electricityh=rowenergy.*el_dailycost;%eur/h
total_electricity_revenue=sum(revenue_electricityh);
totalmethanol_cost=methanol_cost*tot_consum;%eur/y
totpilotoil_cost=tot_pilotoil_cons*pilotoil_cost;%eur/y
alphafactor=(1-(1+k)^(-10))/k;
NPV0=I0_ENGINE;
CAPEX=((-NPV0)/alphafactor);%total yearly capex  of the investment
TOTAL_INVESTMENT=sum(CAPEX./(1+k).^t1);
cash_flowy=(total_electricity_revenue-totalmethanol_cost-OPEX-totpilotoil_cost);%annual cash flow
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
    fprintf('the Discounted Payback Period (DPP) is of  %.2f anni.\n', DPP);
else
    fprintf('the investment is not recoverd before the time.\n');
end
%LEVELIZED COST OF METHANOL aTO HAVE A POSITIVE NPV
tlife=1:25;
CC=CAPEX;
OC=OPEX; 
RO=total_electricity_revenue;
PO=totpilotoil_cost;%COST pilot oil year
LCOM=  (sum((-CC-OC+RO-PO)./(1+k).^t1)+sum((-OC+RO-PO)./(1+k).^t2))/(sum(tot_consum./((1+k).^tlife)));
fprintf('The LCOM is: %.2f eur/ton\n', LCOM);
 % results2 = [results2; percentile_90, LCOM, NPV,cf];
% %end
% % table with the explanation
% results_table = array2table(results2, 'VariableNames', {'Percentile_90', 'LCOM', 'NPV','cf'});
% 
% % save table in file CSV 
% % writetable(results_table, 'results2.xlsx', 'Sheet', 1, 'WriteVariableNames', true);
% 
% disp('Simulazioni completate, risultati salvati in results.csv');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% engine efficency
LHV_MEOH=5.42;%MWh/ton (19.5 mj/kg)
LHV_OIL=11.8;%MWh/ton
eta_engine=(9.73)/((9.73*0.3402*LHV_MEOH)+(0.0084*LHV_OIL*9.73));
