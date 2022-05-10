% Output: Economic Dispatch with UC Decisions (with Reserve)
% February 22, 2022
% Arnab Sur

function results = Output_EDUC


%mkdir('Output_NY22_Aug2019');

define_constants

		
savefilep = {'ed_actual_ny22_fr'};

dstart = 02;                                  % start date
nts = 1;                                       % # of days(a week) 

Total_cost_educ_week = zeros(1,nts);
Total_dcost_educ_week = zeros(1,nts);
Total_rcost_educ_week = zeros(1,nts);
Total_wind_week = zeros(1,nts);
Total_gen_educ_week = zeros(1,nts);
Total_demand_educ_week = zeros(1,nts);

for ctn = dstart:dstart+nts-1
% Loading Data
%data_ed_max = sprintf('scedfr_ny22sh_2019_08_%2.2i_%2.2i', ctn, 23);
%load(data_ed_max);
data_ed_actual = sprintf('scedfr_ny22sh_2019_08_%2.2i_%2.2i', ctn, 1);
load(data_ed_actual);


%data_ed = sprintf('output_dirscucfrseqsh_ny_%4.4i.mat',ctn);
%load(data_uc_ed);

nt = mdo.idx.nt;                                     % # of time periods 

% Indexing Fuel Type Generators
ngen = size(mdo.mpc.gen(:,GEN_BUS));                                 % # of generators 
Index_ng = find(contains(mdo.mpc.genfuel, 'ng'));
Index_rfo = find(contains(mdo.mpc.genfuel, 'rfo'));
Index_nuclear = find(contains(mdo.mpc.genfuel, 'nuclear'));
Index_coal = find(contains(mdo.mpc.genfuel, 'coal'));
Index_hydro = find(contains(mdo.mpc.genfuel, 'hydro'));
Index_other = find(contains(mdo.mpc.genfuel, 'other'));
Index_wood = find(contains(mdo.mpc.genfuel, 'wood'));
Index_wind = find(contains(mdo.mpc.genfuel, 'wind'));
Index_dfo = find(contains(mdo.mpc.genfuel, 'dfo'));
Index_solar = find(contains(mdo.mpc.genfuel, 'solar'));
Index_unknown = find(contains(mdo.mpc.genfuel, 'unknown'));

% Creating Generator per Zone Matrix
nzone = max(mdo.mpc.bus(:,BUS_AREA));                              % # of zones
BusCols = [mdo.mpc.bus(:,BUS_I), mdo.mpc.bus(:,BUS_AREA)];   
GenCol = mdo.mpc.gen(:,1);
BusTable = array2table(BusCols,'VariableNames',{'Bus','Zone'});
GenTable = array2table(GenCol, 'VariableNames', {'Bus'});
ZonalTable = innerjoin(GenTable, BusTable );
Gen_zone = table2array(ZonalTable);


% Dispatched Load
Disp_educ = [];
for i = 1:ngen
    for j = 1:nt
        Disp_educ(i,j) = mdo.results.ExpectedDispatch(i,j);
    end 
end
 
Total_gen_day = sum(sum(Disp_educ(:,:)));
Total_gen_educ_week(ctn-dstart+1) = Total_gen_day;

% Total Demand
demand = zeros(1,nt);
for j =1:nt
    demand(1,j) = sum(mdo.flow(j).mpc.bus(:,PD));
end
Total_demand_day = sum(demand(1,:));
Total_demand_educ_week(ctn-dstart+1) = Total_demand_day;

% Dispatch Cost 
Disp_educ_cost = [];
for j = 1:nt
    for i = 1:ngen
        Disp_educ_cost(i,j) = Disp_educ(i,j)*mdo.mpc.gencost(i,COST); 
    end
end

% Total Dispatch Cost 
Total_dcost_day_educ = sum(sum(Disp_educ_cost(:,:)));
Total_dcost_educ_week(ctn-dstart+1) = Total_dcost_day_educ;

% UC Decision
Decision_educ = [];
for i = 1:ngen
    for j = 1:nt
        Decision_educ(i,j) = mdo.UC.CommitSched(i,j);
    end 
end 


% Reserve Quantity 
Reserve_educ = [];
for i = 1:nt
    Reserve_educ(:,i) = mdo.flow(i).mpc.reserves.R;
end 

% Reserve Cost
Reserve_cost_educ = [];
for j = 1:nt
    for i = 1:ngen
        Reserve_cost_educ(i,j) = Reserve_educ(i,j)*mdo.FixedReserves(j).cost(i);
    end
end 

% Total Reserve Cost
Total_rcost_day_educ = sum(sum(Reserve_cost_educ(:,:)));
Total_rcost_educ_week(ctn-dstart+1) = Total_rcost_day_educ;

% Total Cost
Total_cost_educ = Disp_educ_cost + Reserve_cost_educ;
Total_cost_day_educ = sum(sum(Total_cost_educ(:,:)));
Total_cost_educ_week(ctn-dstart+1) = Total_cost_day_educ;



% Dispatched Load, Dispatch Cost -- Fuel Type
% Nuclear -- Load
Disp_educ_Nuclear = [];
for j = 1:nt
   Disp_educ_Nuclear(j,1) = sum(Disp_educ(Index_nuclear,j));
end  

% Nuclear -- Cost
Disp_educ_cost_Nuclear = [];
for j = 1:nt
   Disp_educ_cost_Nuclear(j,1) = sum(Disp_educ_cost(Index_nuclear,j));
end  

% Nuclear -- Avg Cost
Avg_educ_cost_Nuclear = [];
for j = 1:nt
    Avg_educ_cost_Nuclear(j,1) = Disp_educ_cost_Nuclear(j,1)/ Disp_educ_Nuclear(j,1);
end  

% Nuclear -- Reserve Quantity
Reserve_educ_Nuclear = [];
for j = 1:nt
   Reserve_educ_Nuclear(j,1) = sum(Reserve_educ(Index_nuclear,j));
end  

% Nuclear -- Reserve Cost
Reserve_cost_educ_Nuclear = [];
for j = 1:nt
   Reserve_cost_educ_Nuclear(j,1) = sum(Reserve_cost_educ(Index_nuclear,j));
end  


% NG -- Load
Disp_educ_Ng = [];
for j = 1:nt
   Disp_educ_Ng(j,1) = sum(Disp_educ(Index_ng,j));
end  

% NG -- Cost
Disp_educ_cost_Ng = [];
for j = 1:nt
   Disp_educ_cost_Ng(j,1) = sum(Disp_educ_cost(Index_ng,j));
end  

% NG -- Avg Cost
Avg_educ_cost_Ng = [];
for j = 1:nt
    Avg_educ_cost_Ng(j,1) = Disp_educ_cost_Ng(j,1)/ Disp_educ_Ng(j,1);
end  

% Ng -- Reserve Quantity
Reserve_educ_Ng = [];
for j = 1:nt
   Reserve_educ_Ng(j,1) = sum(Reserve_educ(Index_ng,j));
end  

% Ng -- Reserve Cost
Reserve_cost_educ_Ng = [];
for j = 1:nt
   Reserve_cost_educ_Ng(j,1) = sum(Reserve_cost_educ(Index_ng,j));
end  


% Dfo -- Load
Disp_educ_Dfo = [];
for j = 1:nt
   Disp_educ_Dfo(j,1) = sum(Disp_educ(Index_dfo,j));
end  


% Rfo -- Load
Disp_educ_Rfo = [];
for j = 1:nt
   Disp_educ_Rfo(j,1) = sum(Disp_educ(Index_rfo,j));
end  

% Rfo -- Cost
Disp_educ_cost_Rfo = [];
for j = 1:nt
   Disp_educ_cost_Rfo(j,1) = sum(Disp_educ_cost(Index_rfo,j));
end 

% Rfo -- Avg Cost
Avg_educ_cost_Rfo = [];
for j = 1:nt
    Avg_educ_cost_Rfo(j,1) = Disp_educ_cost_Rfo(j,1)/ Disp_educ_Rfo(j,1);
end  

% Rfo -- Reserve Quantity
Reserve_educ_Rfo = [];
for j = 1:nt
   Reserve_educ_Rfo(j,1) = sum(Reserve_educ(Index_rfo,j));
end  

% Rfo -- Reserve Cost
Reserve_cost_educ_Rfo = [];
for j = 1:nt
   Reserve_cost_educ_Rfo(j,1) = sum(Reserve_cost_educ(Index_rfo,j));
end  

% Hydro -- Load
Disp_educ_Hydro = [];
for j = 1:nt
   Disp_educ_Hydro(j,1) = sum(Disp_educ(Index_hydro,j));
end  

% Hydro -- Cost
Disp_educ_cost_Hydro = [];
for j = 1:nt
   Disp_educ_cost_Hydro(j,1) = sum(Disp_educ_cost(Index_hydro,j));
end  

% Hydro -- Avg Cost
Avg_educ_cost_Hydro = [];
for j = 1:nt
    Avg_educ_cost_Hydro(j,1) = Disp_educ_cost_Hydro(j,1)/ Disp_educ_Hydro(j,1);
end  

% Hydro -- Reserve Quantity
Reserve_educ_Hydro = [];
for j = 1:nt
   Reserve_educ_Hydro(j,1) = sum(Reserve_educ(Index_hydro,j));
end  

% Hydro -- Reserve Cost
Reserve_cost_educ_Hydro = [];
for j = 1:nt
   Reserve_cost_educ_Hydro(j,1) = sum(Reserve_cost_educ(Index_hydro,j));
end  

% Coal -- Load
Disp_educ_Coal = [];
for j = 1:nt
   Disp_educ_Coal(j,1) = sum(Disp_educ(Index_coal,j));
end  

% Coal -- Cost
Disp_educ_cost_Coal = [];
for j = 1:nt
   Disp_educ_cost_Coal(j,1) = sum(Disp_educ_cost(Index_coal,j));
end  

% Coal -- Avg Cost
Avg_educ_cost_Coal = [];
for j = 1:nt
    Avg_educ_cost_Coal(j,1) = Disp_educ_cost_Coal(j,1)/ Disp_educ_Coal(j,1);
end  

% Coal -- Reserve Quantity
Reserve_educ_Coal = [];
for j = 1:nt
   Reserve_educ_Coal(j,1) = sum(Reserve_educ(Index_coal,j));
end  

% Coal -- Reserve Cost
Reserve_cost_educ_Coal = [];
for j = 1:nt
   Reserve_cost_educ_Coal(j,1) = sum(Reserve_cost_educ(Index_coal,j));
end  

% Wood -- Load
Disp_educ_Wood = [];
for j=1:nt
   Disp_educ_Wood(j,1) = sum(Disp_educ(Index_wood,j));
end  

% Wood -- Cost
Disp_educ_cost_Wood = [];
for j=1:nt
   Disp_educ_cost_Wood(j,1) = sum(Disp_educ_cost(Index_wood,j));
end  

% Wood -- Avg Cost
Avg_educ_cost_Wood = [];
for j = 1:nt
    Avg_educ_cost_Wood(j,1) = Disp_educ_cost_Wood(j,1)/ Disp_educ_Wood(j,1);
end  

% Wood -- Reserve Quantity
Reserve_educ_Wood = [];
for j = 1:nt
   Reserve_educ_Wood(j,1) = sum(Reserve_educ(Index_wood,j));
end  

% Wood -- Reserve Cost
Reserve_cost_educ_Wood = [];
for j = 1:nt
   Reserve_cost_educ_Wood(j,1) = sum(Reserve_cost_educ(Index_wood,j));
end  

% Wind -- Load
Disp_educ_Wind = [];
for j=1:nt
   Disp_educ_Wind(j,1) = sum(Disp_educ(Index_wind,j));
end  

% Total Wind Dispatch
Total_wind_day = sum(Disp_educ_Wind);
Total_wind_week(ctn-dstart+1) = Total_wind_day;

% Wind -- Cost
Disp_educ_cost_Wind = [];
for j = 1:nt
   Disp_educ_cost_Wind(j,1) = sum(Disp_educ_cost(Index_wind,j));
end  

% Wind -- Avg Cost
Avg_educ_cost_Wind = [];
for j = 1:nt
    Avg_educ_cost_Wind(j,1) = Disp_educ_cost_Wind(j,1)/ Disp_educ_Wind(j,1);
end  

% Wind -- Reserve Quantity
Reserve_educ_Wind = [];
for j = 1:nt
   Reserve_educ_Wind(j,1) = sum(Reserve_educ(Index_wind,j));
end  

% Wind -- Reserve Cost
Reserve_cost_educ_Wind = [];
for j = 1:nt
   Reserve_cost_educ_Wind(j,1) = sum(Reserve_cost_educ(Index_wind,j));
end  

% Solar -- Load
Disp_educ_Solar = [];
for j = 1:nt
   Disp_educ_Solar(j,1) = sum(Disp_educ(Index_solar,j));
end  
% Total Solar Dispatch
Total_solar_day = sum(Disp_educ_Solar);
Total_solar_week(ctn-dstart+1) = Total_solar_day;


% Other -- Load
Disp_educ_Other = [];
for j = 1:nt
   Disp_educ_Other(j,1) = sum(Disp_educ(Index_other,j));
end  

% Other -- Cost
Disp_educ_cost_Other = [];
for j=1:nt
   Disp_educ_cost_Other(j,1) = sum(Disp_educ_cost(Index_other,j));
end  

% Other -- Avg Cost
Avg_educ_cost_Other = [];
for j = 1:nt
    Avg_educ_cost_Other(j,1) = Disp_educ_cost_Other(j,1)/ Disp_educ_Other(j,1);
end  

% Other -- Reserve Quantity
Reserve_educ_Other = [];
for j = 1:nt
   Reserve_educ_Other(j,1) = sum(Reserve_educ(Index_other,j));
end  

% Other -- Reserve Cost
Reserve_cost_educ_Other = [];
for j = 1:nt
   Reserve_cost_educ_Other(j,1) = sum(Reserve_cost_educ(Index_other,j));
end  

% Unknown -- Load
Disp_educ_Unknown = [];
for j = 1:nt
   Disp_educ_Unknown(j,1) = sum(Disp_educ(Index_unknown,j));
end  

% Net Demand
Net_demand = [];
for j = 1:nt
   Net_demand(j,1) = demand(1,j) - (Disp_educ_Wind(j) + Disp_educ_Solar(j));
end 


% Dispatched Load, Dispatch Cost -- Zonal
% Dispatched Load--Zonal
Disp_educ_zonal = zeros(nzone,nt);
for i = 1:nt
    for k = 1:nzone
        index = find(Gen_zone(:,2) == k);
        for j = 1: length(index)
            Disp_educ_zonal(k,i) = Disp_educ_zonal(k,i) + Disp_educ(index(j),i);
        end
    end 
end

% Didpatch Cost--Zonal
Disp_educ_cost_zonal = zeros(nzone,nt);
for i=1:nt
    for k = 1:nzone
        index = find(Gen_zone(:,2) == k);
        for j = 1: length(index)
            Disp_educ_cost_zonal(k,i) = Disp_educ_cost_zonal(k,i) + Disp_educ_cost(index(j),i);
        end
    end 
end

% Average Cost-- Zonal
Average_educ_cost_zonal = zeros(nzone,nt);
for i = 1:nzone
    for j = 1:nt
        Average_educ_cost_zonal (i,j) = Disp_educ_cost_zonal(i,j)/Disp_educ_zonal(i,j);
    end
end

% Reserve Quantity--Zonal
Reserve_educ_zonal = zeros(nzone,nt);
for i = 1:nt
    for k = 1:nzone
        index = find(Gen_zone(:,2) == k);
        for j = 1: length(index)
            Reserve_educ_zonal(k,i) = Reserve_educ_zonal(k,i) + Reserve_educ(index(j),i);
        end
    end 
end

% Reserve Cost--Zonal
Reserve_cost_educ_zonal = zeros(nzone,nt);
for i=1:nt
    for k = 1:nzone
        index = find(Gen_zone(:,2) == k);
        for j = 1: length(index)
            Reserve_cost_educ_zonal(k,i) = Reserve_cost_educ_zonal(k,i) + Reserve_cost_educ(index(j),i);
        end
    end 
end


% UC Decision -- Fuel Type
%Nuclear
Decision_educ_Nuclear = [];
for j = 1:nt
    for i = 1:length(Index_nuclear)
    Decision_educ_Nuclear(i,j) = Decision_educ(Index_nuclear(i),j);
    end 
end

% NG
Decision_educ_Ng = [];
for j = 1:nt
    for i = 1:length(Index_ng)
    Decision_educ_Ng(i,j) = Decision_educ(Index_ng(i),j);
    end 
end


% RFO
Decision_educ_Rfo = [];
for j = 1:nt
    for i = 1:length(Index_rfo)
    Decision_educ_Rfo(i,j) = Decision_educ(Index_rfo(i),j);
    end 
end


% Hydro 
Decision_educ_Hydro = [];
for j = 1:nt
    for i = 1:length(Index_hydro)
    Decision_educ_Hydro(i,j) = Decision_educ(Index_hydro(i),j);
    end 
end


% Coal 
Decision_educ_Coal = [];
for j = 1:nt
    for i = 1:length(Index_coal)
    Decision_educ_Coal(i,j) = Decision_educ(Index_coal(i),j);
    end 
end

% Wood 
Decision_educ_Wood = [];
for j = 1:nt
    for i = 1:length(Index_wood)
    Decision_educ_Wood(i,j) = Decision_educ(Index_wood(i),j);
    end 
end


% Wind
Decision_educ_Wind = [];
for j = 1:nt
    for i = 1:length(Index_wind)
    Decision_educ_Wind(i,j) = Decision_educ(Index_wind(i),j);
    end 
end


% Other
Decision_educ_Other = [];
for j = 1:nt
    for i = 1:length(Index_other)
    Decision_educ_Other(i,j) = Decision_educ(Index_other(i),j);
    end 
end
 

cd('Output_NY22_Aug2019')

% Saving the Results 
%save(sprintf('%s_2019_08_%2.2i_%2.2i','output_', savefilep{1}, ctn, jsel),"Disp_educ","Disp_educ_cost","Decision_educ",...
    %"Disp_educ_Nuclear","Disp_educ_cost_Nuclear","Decision_educ_Nuclear","Disp_educ_Ng","Disp_educ_cost_Ng",...
   % "Decision_educ_Ng","Disp_educ_Rfo","Disp_educ_cost_Rfo","Decision_educ_Rfo","Disp_educ_Hydro",...
   % "Disp_educ_cost_Hydro","Decision_educ_Hydro","Disp_educ_Coal","Disp_educ_cost_Coal","Decision_educ_Coal",...
  %  "Disp_educ_Wood","Disp_educ_cost_Wood","Decision_educ_Wood","Disp_educ_Wind","Disp_educ_cost_Wind",...
  %  "Decision_educ_Wind","Disp_educ_Other","Disp_educ_cost_Other","Decision_educ_Other","Disp_educ_zonal",...
  %  "Disp_educ_cost_zonal", "Average_educ_cost_zonal", "Avg_educ_cost_Nuclear", "Avg_educ_cost_Ng",...
  %  "Avg_educ_cost_Rfo","Avg_educ_cost_Hydro","Avg_educ_cost_Coal","Avg_educ_cost_Wood","Avg_educ_cost_Wind",...
  %  "Avg_educ_cost_Other","Reserve_cost_educ","Reserve_educ","Total_cost_educ","Reserve_cost_educ_zonal",...
  %  "Reserve_educ_zonal","Reserve_cost_educ_Other","Reserve_educ_Other","Reserve_cost_educ_Wind",...
  %  "Reserve_educ_Wind","Reserve_cost_educ_Wood","Reserve_educ_Wood","Reserve_cost_educ_Coal",...
  %  "Reserve_educ_Coal","Reserve_cost_educ_Hydro","Reserve_educ_Hydro","Reserve_cost_educ_Rfo",...
  %  "Reserve_educ_Rfo","Reserve_cost_educ_Ng","Reserve_educ_Ng","Reserve_cost_educ_Nuclear","Reserve_educ_Nuclear",...
  %  "Total_cost_day_educ","Total_wind_day","Total_dcost_day_educ","Total_rcost_day_educ","Total_gen_day",...
   % "Total_demand_day"); 





% Writing Excel Sheets
currentFolder = sprintf('Excel_educ_actual_%2.2i',ctn);
%currentFolder = sprintf('Excel_educ_max_%2.2i',ctn);
mkdir(currentFolder)
cd(currentFolder)
% Dispatched Load Table
Disp_educ_ft = [transpose(mdo.mpc.genfuel); transpose(num2cell(Disp_educ))];  % Required to be same data type
T = cell2table(Disp_educ_ft);
writetable(T,'Disp_educ_fr.csv');

% Didpatch Cost Table
Disp_educ_cost_ft = [transpose(mdo.mpc.genfuel);  transpose(num2cell(Disp_educ_cost))];
T = cell2table(Disp_educ_cost_ft);
writetable(T,'Disp_educ_cost_fr.csv');

% Hourly and Net Demand Table
D = [transpose(demand) Disp_educ_Wind Disp_educ_Solar Net_demand];

D_title = array2table(D, 'VariableNames', {'Demand', 'Wind', 'Solar', 'Net Demand'});
writetable(D_title, 'Hourly_demand.csv');

% Reserve Quantity Table
Reserve_educ_ft = [transpose(mdo.mpc.genfuel); transpose(num2cell(Reserve_educ))];  
T = cell2table(Reserve_educ_ft);
writetable(T,'Reserve_educ_fr.csv');

% Reserve Cost Table
Reserve_cost_educ_ft = [transpose(mdo.mpc.genfuel); transpose(num2cell(Reserve_cost_educ))];
T = cell2table(Reserve_cost_educ_ft);
writetable(T,'Reserve_cost_educ_fr.csv');

% Total Cost Table
Total_cost_educ_ft = [transpose(mdo.mpc.genfuel);  transpose(num2cell(Total_cost_educ))];
T = cell2table(Total_cost_educ_ft);
writetable(T,'Total_cost_educ_fr.csv');


% Decision Table
%Decision_educ_ft = [mdo.mpc.genfuel num2cell(Decision_educ)];
%T = cell2table(Decision_educ_ft);
%writetable(T,'Decision_educ_fr.csv');


% Dispatched Load Table -- Fuel Type
A = [Disp_educ_Nuclear  Disp_educ_Ng Disp_educ_Dfo Disp_educ_Rfo Disp_educ_Hydro  Disp_educ_Coal Disp_educ_Wood  ...
     Disp_educ_Wind Disp_educ_Solar Disp_educ_Other Disp_educ_Unknown];

T_ft = array2table(A, 'VariableNames', {'Nuclear', 'NG', 'DFO', 'RFO', 'Hydro','Coal', 'Wood', 'Wind', ...
    'Solar', 'Other', 'Unknown'});
writetable(T_ft, 'Disp_educ_fr_fueltype.csv');

% Dispatch Cost Table -- Fuel Type
C = [Disp_educ_cost_Nuclear Avg_educ_cost_Nuclear Disp_educ_cost_Ng Avg_educ_cost_Ng Disp_educ_cost_Rfo ... 
    Avg_educ_cost_Rfo Disp_educ_cost_Hydro Avg_educ_cost_Hydro Disp_educ_cost_Coal Avg_educ_cost_Coal ...
    Disp_educ_cost_Wood Avg_educ_cost_Wood Disp_educ_cost_Wind Avg_educ_cost_Wind Disp_educ_cost_Other ...
    Avg_educ_cost_Other];

T_cost_ft = array2table(C, 'VariableNames', {'Nuclear', 'Avg_Nuclear', 'NG', 'Avg_NG', 'RFO', 'Avg_RFO',...
    'Hydro', 'Avg_Hydro', 'Coal', 'Avg_Coal', 'Wood', 'Avg_wood', 'Wind', 'Avg_Wind', 'Other', 'Avg_Other'});
writetable(T_cost_ft, 'Disp_educ_cost_fr_fueltype.csv');

% Reserve Quantity Table -- Fuel Type
 R = [Reserve_educ_Nuclear  Reserve_educ_Ng  Reserve_educ_Rfo Reserve_educ_Hydro  Reserve_educ_Coal  Reserve_educ_Wood  ...
     Reserve_educ_Wind  Reserve_educ_Other];

R_ft = array2table(R, 'VariableNames', {'Nuclear', 'NG', 'RFO',  'Hydro','Coal', 'Wood', 'Wind', 'Other'});
writetable(R_ft, 'Reserve_educ_fr_fueltype.csv');

% Reserve Cost Table -- Fuel Type
R_cost = [Reserve_cost_educ_Nuclear Reserve_cost_educ_Ng Reserve_cost_educ_Rfo Reserve_cost_educ_Hydro ...
    Reserve_cost_educ_Coal Reserve_cost_educ_Wood Reserve_cost_educ_Wind Reserve_cost_educ_Other];

R_cost_ft = array2table(R_cost, 'VariableNames', {'Nuclear', 'NG', 'RFO', 'Hydro', 'Coal', 'Wood', 'Wind', 'Other'});
writetable(R_cost_ft, 'Reserve_cost_educ_fr_fueltype.csv');


% Dispatched Load Table -- Zonal
T_zone = array2table(transpose(Disp_educ_zonal));
T_zone.Properties.VariableNames(1:nzone) = {'Zone1', 'Zone2', 'Zone3', 'Zone4', 'Zone5', 'Zone6',....
                       'Zone7', 'Zone8', 'Zone9', 'Zone10', 'Zone11'};
writetable(T_zone, 'Disp_educ_fr_zonal.csv');

% Dispatch Cost Table -- Zonal
T_cost_zone = array2table(transpose(Disp_educ_cost_zonal));
T_cost_zone.Properties.VariableNames(1:nzone) = {'Zone1', 'Zone2', 'Zone3', 'Zone4', 'Zone5', 'Zone6',....
                       'Zone7', 'Zone8', 'Zone9', 'Zone10', 'Zone11'};
writetable(T_cost_zone, 'Disp_educ_cost_fr_zonal.csv');

% Reserve Quantity Table -- Zonal
R_zone = array2table(transpose(Reserve_educ_zonal));
R_zone.Properties.VariableNames(1:nzone) = {'Zone1', 'Zone2', 'Zone3', 'Zone4', 'Zone5', 'Zone6',....
                       'Zone7', 'Zone8', 'Zone9', 'Zone10', 'Zone11'};
writetable(R_zone, 'Reserve_educ_fr_zonal.csv');

% Reserve Cost Table -- Zonal
R_cost_zone = array2table(transpose(Reserve_cost_educ_zonal));
R_cost_zone.Properties.VariableNames(1:nzone) = {'Zone1', 'Zone2', 'Zone3', 'Zone4', 'Zone5', 'Zone6',....
                       'Zone7', 'Zone8', 'Zone9', 'Zone10', 'Zone11'};
writetable(R_cost_zone, 'Reserve_cost_educ_fr_zonal.csv');

% Average Cost Table -- Zonal
T_avg_zone = array2table(transpose(Average_educ_cost_zonal));
T_avg_zone.Properties.VariableNames(1:nzone) = {'Zone1', 'Zone2', 'Zone3', 'Zone4', 'Zone5', 'Zone6',....
                       'Zone7', 'Zone8', 'Zone9', 'Zone10', 'Zone11'};
writetable(T_avg_zone, 'Average_educ_fr_zonal.csv');


% Decision Table -- Fuel Type
%csvwrite('Decision_educ_fr_Nuclear.csv', Decision_educ_Nuclear);
%csvwrite('Decision_educ_fr_Ng.csv', Decision_educ_Ng);
%csvwrite('Decision_educ_fr_Rfo.csv', Decision_educ_Rfo);
%csvwrite('Decision_educ_fr_Hydro.csv', Decision_educ_Hydro);
%csvwrite('Decision_educ_fr_Coal.csv', Decision_educ_Coal);
%csvwrite('Decision_educ_fr_Wood.csv', Decision_educ_Wood);
%csvwrite('Decision_educ_fr_Wind.csv', Decision_educ_Wind);
%csvwrite('Decision_educ_fr_Other.csv', Decision_educ_Other);

cd ..
cd ..
end

cd('Output_NY22_Aug2019')
currentFolder = sprintf('Excel_educ_actual_%2.2i',ctn);
%currentFolder = sprintf('Excel_educ_max_%2.2i',ctn);
mkdir(currentFolder)
cd(currentFolder)
 %Weekly Cost Table 
C_week = transpose([Total_gen_educ_week; Total_demand_educ_week; Total_wind_week; Total_solar_week;...
         Total_dcost_educ_week; Total_rcost_educ_week; Total_cost_educ_week]);
C_week = array2table(C_week,'VariableNames', {'Tot Gen(MW)', 'Tot Demand(MW)','Wind Gen(MW)','Solar Gen(MW)'...
                    'Disp Cost($)','Res Cost($)','Tot Cost($)'}, 'RowNames', {'Day1'});
writetable(C_week, 'Overall_results_week.csv','WriteRowNames',1); 

% Weekly Wind Dispatch    
%W_dis = transpose(Total_wind_week);
%W_dis = array2table(W_dis,'VariableNames', {'Wind Gen'}, 'RowNames',...
 %                          {'Day1', 'Day2', 'Day3', 'Day4', 'Day5', 'Day6','Day7'});
%writetable(W_dis, 'Wind_dispatch_week.csv','WriteRowNames',1); 

cd ..
cd ..

