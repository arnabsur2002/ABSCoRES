% Output: Multiple ED with UC Decisions (with Reserves)
% Arnab Sur
% March 16,2022

mkdir('Output_NY22_Aug2019_mod');
define_constants

savefilepscen = {'ed_scen_ny22_fr'};		
savefilep = {'ed_ny22_fr'};

dstart = 02;                                  % start date
nts = 1;                                       % # of days(a week) 
nscen = 50;
%nzone = 11; 

Total_cost_educ_week = zeros(1,nts);
Total_dcost_educ_week = zeros(1,nts);
Total_rcost_educ_week = zeros(1,nts);
Total_wind_week = zeros(1,nts);
Total_gen_educ_week = zeros(1,nts);
Total_demand_educ_week = zeros(1,nts);



for ctn = dstart:dstart+nts-1
% Loading Data
for jsel = [2:37,39:51];
data_ed_scen = sprintf('scedfr_ny22sh_2019_08_%2.2i_%2.2i', ctn, jsel);
load(data_ed_scen);


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
Index_solar = find(contains(mdo.mpc.genfuel, 'solar'));
Index_unknown = find(contains(mdo.mpc.genfuel, 'unknown'));
Index_dfo = find(contains(mdo.mpc.genfuel, 'dfo'));



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
 
Total_gen_day_scen(jsel) = sum(sum(Disp_educ(:,:)));
%Total_gen_educ_week(ctn-dstart+1) = Total_gen_day;


% Dispatch Cost 
Disp_educ_cost = [];
for j = 1:nt
    for i = 1:ngen
        Disp_educ_cost(i,j) = Disp_educ(i,j)*mdo.mpc.gencost(i,COST); 
    end
end

% Total Dispatch Cost 
Total_dcost_day_educ_scen(jsel) = sum(sum(Disp_educ_cost(:,:)));
%Total_dcost_educ_week(ctn-dstart+1) = Total_dcost_day_educ;

% Total and Hourly Demand
demand = zeros(1,nt);
for j =1:nt
    demand(1,j) = sum(mdo.flow(j).mpc.bus(:,PD));
end
Total_demand_day_scen(jsel) = sum(demand(1,:));
%Total_demand_educ_week(ctn-dstart+1) = Total_demand_day;



for j = 1:nt
    Demand_hourly_scen(jsel,j) = demand(1,j);
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
Total_rcost_day_educ_scen(jsel) = sum(sum(Reserve_cost_educ(:,:)));
%Total_rcost_educ_week(ctn-dstart+1) = Total_rcost_day_educ;


% Total Cost
Total_cost_educ = Disp_educ_cost + Reserve_cost_educ;
Total_cost_day_educ_scen(jsel) = sum(sum(Total_cost_educ(:,:)));
%Total_cost_educ_week(ctn-dstart+1) = Total_cost_day_educ;


% Dispatched Load, Dispatch Cost -- Fuel Type
% Nuclear -- Total and Hourly Load
Disp_educ_Nuclear = [];
for j = 1:nt
   Disp_educ_Nuclear(j,1) = sum(Disp_educ(Index_nuclear,j));
end  

Disp_educ_day_Nuclear = sum(Disp_educ_Nuclear);

for j = 1:nt
    Disp_educ_hourly_scen_Nuclear(jsel,j) = Disp_educ_Nuclear(j);
end

% Nuclear -- Cost
Disp_educ_cost_Nuclear = [];
for j = 1:nt
   Disp_educ_cost_Nuclear(j,1) = sum(Disp_educ_cost(Index_nuclear,j));
end  

Disp_educ_cost_day_Nuclear = sum(Disp_educ_cost_Nuclear);
Disp_educ_cost_day_scen_Nuclear(jsel) = Disp_educ_cost_day_Nuclear;

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

Reserve_educ_day_Nuclear = sum(Reserve_educ_Nuclear);
Reserve_educ_day_scen_Nuclear(jsel) = Reserve_educ_day_Nuclear;

% Nuclear -- Reserve Cost
Reserve_cost_educ_Nuclear = [];
for j = 1:nt
   Reserve_cost_educ_Nuclear(j,1) = sum(Reserve_cost_educ(Index_nuclear,j));
end  

Reserve_cost_educ_day_Nuclear = sum(Reserve_cost_educ_Nuclear);
Reserve_cost_educ_day_scen_Nuclear(jsel) = Reserve_cost_educ_day_Nuclear;

% Ng -- Total and Hourly Load
Disp_educ_Ng = [];
for j = 1:nt
   Disp_educ_Ng(j,1) = sum(Disp_educ(Index_ng,j));
end  

Disp_educ_day_Ng = sum(Disp_educ_Ng);
for j = 1:nt
    Disp_educ_hourly_scen_Ng(jsel,j) = Disp_educ_Ng(j);
end

% Ng -- Cost
Disp_educ_cost_Ng = [];
for j = 1:nt
   Disp_educ_cost_Ng(j,1) = sum(Disp_educ_cost(Index_ng,j));
end  

Disp_educ_cost_day_Ng = sum(Disp_educ_cost_Ng);
Disp_educ_cost_day_scen_Ng(jsel) = Disp_educ_cost_day_Ng;

% Ng -- Avg Cost
Avg_educ_cost_Ng = [];
for j = 1:nt
    Avg_educ_cost_Ng(j,1) = Disp_educ_cost_Ng(j,1)/ Disp_educ_Ng(j,1);
end  

% Ng -- Reserve Quantity
Reserve_educ_Ng = [];
for j = 1:nt
   Reserve_educ_Ng(j,1) = sum(Reserve_educ(Index_ng,j));
end  

Reserve_educ_day_Ng = sum(Reserve_educ_Ng);
Reserve_educ_day_scen_Ng(jsel) = Reserve_educ_day_Ng;

% Ng -- Reserve Cost
Reserve_cost_educ_Ng = [];
for j = 1:nt
   Reserve_cost_educ_Ng(j,1) = sum(Reserve_cost_educ(Index_ng,j));
end  

Reserve_cost_educ_day_Ng = sum(Reserve_cost_educ_Ng);
Reserve_cost_educ_day_scen_Ng(jsel) = Reserve_cost_educ_day_Ng;


% Dfo -- Total and Hourly Load
Disp_educ_Dfo = [];
for j = 1:nt
   Disp_educ_Dfo(j,1) = sum(Disp_educ(Index_dfo,j));
end  

Disp_educ_day_Dfo = sum(Disp_educ_Dfo);
for j = 1:nt
    Disp_educ_hourly_scen_Dfo(jsel,j) = Disp_educ_Dfo(j);
end

% Dfo -- Cost
Disp_educ_cost_Dfo = [];
for j = 1:nt
   Disp_educ_cost_Dfo(j,1) = sum(Disp_educ_cost(Index_dfo,j));
end  

Disp_educ_cost_day_Dfo = sum(Disp_educ_cost_Dfo);
Disp_educ_cost_day_scen_Dfo(jsel) = Disp_educ_cost_day_Dfo;

% Dfo -- Avg Cost
Avg_educ_cost_Dfo = [];
for j = 1:nt
    Avg_educ_cost_Dfo(j,1) = Disp_educ_cost_Dfo(j,1)/ Disp_educ_Dfo(j,1);
end  

% Dfo -- Reserve Quantity
Reserve_educ_Dfo = [];
for j = 1:nt
   Reserve_educ_Dfo(j,1) = sum(Reserve_educ(Index_dfo,j));
end  

Reserve_educ_day_Dfo = sum(Reserve_educ_Dfo);
Reserve_educ_day_scen_Dfo(jsel) = Reserve_educ_day_Dfo;

% Dfo -- Reserve Cost
Reserve_cost_educ_Dfo = [];
for j = 1:nt
   Reserve_cost_educ_Dfo(j,1) = sum(Reserve_cost_educ(Index_dfo,j));
end  

Reserve_cost_educ_day_Dfo = sum(Reserve_cost_educ_Dfo);
Reserve_cost_educ_day_scen_Dfo(jsel) = Reserve_cost_educ_day_Dfo;


% Rfo -- Total and Hourly Load
Disp_educ_Rfo = [];
for j = 1:nt
   Disp_educ_Rfo(j,1) = sum(Disp_educ(Index_rfo,j));
end  

Disp_educ_day_Rfo = sum(Disp_educ_Rfo);
for j = 1:nt
    Disp_educ_hourly_scen_Rfo(jsel,j) = Disp_educ_Rfo(j);
end

% Rfo -- Cost
Disp_educ_cost_Rfo = [];
for j = 1:nt
   Disp_educ_cost_Rfo(j,1) = sum(Disp_educ_cost(Index_rfo,j));
end  

Disp_educ_cost_day_Rfo = sum(Disp_educ_cost_Rfo);
Disp_educ_cost_day_scen_Rfo(jsel) = Disp_educ_cost_day_Rfo;

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

Reserve_educ_day_Rfo = sum(Reserve_educ_Rfo);
Reserve_educ_day_scen_Rfo(jsel) = Reserve_educ_day_Rfo;

% Rfo -- Reserve Cost
Reserve_cost_educ_Rfo = [];
for j = 1:nt
   Reserve_cost_educ_Rfo(j,1) = sum(Reserve_cost_educ(Index_rfo,j));
end  

Reserve_cost_educ_day_Rfo = sum(Reserve_cost_educ_Rfo);
Reserve_cost_educ_day_scen_Rfo(jsel) = Reserve_cost_educ_day_Rfo;


% Coal -- Total and Hourly Load
Disp_educ_Coal = [];
for j = 1:nt
   Disp_educ_Coal(j,1) = sum(Disp_educ(Index_coal,j));
end  

Disp_educ_day_Coal = sum(Disp_educ_Coal);
for j = 1:nt
    Disp_educ_hourly_scen_Coal(jsel,j) = Disp_educ_Coal(j);
end

% Coal -- Cost
Disp_educ_cost_Coal = [];
for j = 1:nt
   Disp_educ_cost_Coal(j,1) = sum(Disp_educ_cost(Index_coal,j));
end  

Disp_educ_cost_day_Coal = sum(Disp_educ_cost_Coal);
Disp_educ_cost_day_scen_Coal(jsel) = Disp_educ_cost_day_Coal;

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

Reserve_educ_day_Coal = sum(Reserve_educ_Coal);
Reserve_educ_day_scen_Coal(jsel) = Reserve_educ_day_Coal;

% Coal -- Reserve Cost
Reserve_cost_educ_Coal = [];
for j = 1:nt
   Reserve_cost_educ_Coal(j,1) = sum(Reserve_cost_educ(Index_coal,j));
end  

Reserve_cost_educ_day_Coal = sum(Reserve_cost_educ_Coal);
Reserve_cost_educ_day_scen_Coal(jsel) = Reserve_cost_educ_day_Coal;


% Wood -- Total and Hourly Load
Disp_educ_Wood = [];
for j = 1:nt
   Disp_educ_Wood(j,1) = sum(Disp_educ(Index_wood,j));
end  

Disp_educ_day_Wood = sum(Disp_educ_Wood);
for j = 1:nt
    Disp_educ_hourly_scen_Wood(jsel,j) = Disp_educ_Wood(j);
end

% Wood -- Cost
Disp_educ_cost_Wood = [];
for j = 1:nt
   Disp_educ_cost_Wood(j,1) = sum(Disp_educ_cost(Index_wood,j));
end  

Disp_educ_cost_day_Wood = sum(Disp_educ_cost_Wood);
Disp_educ_cost_day_scen_Wood(jsel) = Disp_educ_cost_day_Wood;

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

Reserve_educ_day_Wood = sum(Reserve_educ_Wood);
Reserve_educ_day_scen_Wood(jsel) = Reserve_educ_day_Wood;

% Wood -- Reserve Cost
Reserve_cost_educ_Wood = [];
for j = 1:nt
   Reserve_cost_educ_Wood(j,1) = sum(Reserve_cost_educ(Index_wood,j));
end  

Reserve_cost_educ_day_Wood = sum(Reserve_cost_educ_Wood);
Reserve_cost_educ_day_scen_Wood(jsel) = Reserve_cost_educ_day_Wood;



% Wind -- Total and Hourly Load
Disp_educ_Wind = [];
for j = 1:nt
   Disp_educ_Wind(j,1) = sum(Disp_educ(Index_wind,j));
end  

Disp_educ_day_Wind(jsel) = sum(Disp_educ_Wind(:));

for j = 1:nt
    Disp_educ_hourly_scen_Wind(jsel,j) = Disp_educ_Wind(j);
end

% Wind -- Cost
Disp_educ_cost_Wind = [];
for j = 1:nt
   Disp_educ_cost_Wind(j,1) = sum(Disp_educ_cost(Index_wind,j));
end  

Disp_educ_cost_day_Wind = sum(Disp_educ_cost_Wind);
Disp_educ_cost_day_scen_Wind(jsel) = Disp_educ_cost_day_Wind;

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

Reserve_educ_day_Wind = sum(Reserve_educ_Wind);
Reserve_educ_day_scen_Wind(jsel) = Reserve_educ_day_Wind;

% Wind -- Reserve Cost
Reserve_cost_educ_Wind = [];
for j = 1:nt
   Reserve_cost_educ_Wind(j,1) = sum(Reserve_cost_educ(Index_wind,j));
end  

Reserve_cost_educ_day_Wind = sum(Reserve_cost_educ_Wind);
Reserve_cost_educ_day_scen_Wind(jsel) = Reserve_cost_educ_day_Wind;


% Solar -- Total and Hourly Load
Disp_educ_Solar = [];
for j = 1:nt
   Disp_educ_Solar(j,1) = sum(Disp_educ(Index_solar,j));
end  

Disp_educ_day_Solar(jsel) = sum(Disp_educ_Solar(:));

for j = 1:nt
    Disp_educ_hourly_scen_Solar(jsel,j) = Disp_educ_Solar(j);
end


% Solar -- Cost
Disp_educ_cost_Solar = [];
for j = 1:nt
   Disp_educ_cost_Solar(j,1) = sum(Disp_educ_cost(Index_solar,j));
end  

Disp_educ_cost_day_Solar = sum(Disp_educ_cost_Solar);
Disp_educ_cost_day_scen_Solar(jsel) = Disp_educ_cost_day_Solar;

% Solar -- Avg Cost
Avg_educ_cost_Solar = [];
for j = 1:nt
    Avg_educ_cost_Solar(j,1) = Disp_educ_cost_Solar(j,1)/ Disp_educ_Solar(j,1);
end  

% Solar -- Reserve Quantity
Reserve_educ_Solar = [];
for j = 1:nt
   Reserve_educ_Solar(j,1) = sum(Reserve_educ(Index_solar,j));
end  

Reserve_educ_day_Solar = sum(Reserve_educ_Solar);
Reserve_educ_day_scen_Solar(jsel) = Reserve_educ_day_Solar;

% Solar -- Reserve Cost
Reserve_cost_educ_Solar = [];
for j = 1:nt
   Reserve_cost_educ_Solar(j,1) = sum(Reserve_cost_educ(Index_solar,j));
end  

Reserve_cost_educ_day_Solar = sum(Reserve_cost_educ_Solar);
Reserve_cost_educ_day_scen_Solar(jsel) = Reserve_cost_educ_day_Solar;



% Hydro -- Total and Hourly Load
Disp_educ_Hydro = [];
for j = 1:nt
   Disp_educ_Hydro(j,1) = sum(Disp_educ(Index_hydro,j));
end  

Disp_educ_day_Hydro = sum(Disp_educ_Hydro);
for j = 1:nt
    Disp_educ_hourly_scen_Hydro(jsel,j) = Disp_educ_Hydro(j);
end

% Hydro -- Cost
Disp_educ_cost_Hydro = [];
for j = 1:nt
   Disp_educ_cost_Hydro(j,1) = sum(Disp_educ_cost(Index_hydro,j));
end  

Disp_educ_cost_day_Hydro = sum(Disp_educ_cost_Hydro);
Disp_educ_cost_day_scen_Hydro(jsel) = Disp_educ_cost_day_Hydro;

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

Reserve_educ_day_Hydro = sum(Reserve_educ_Solar);
Reserve_educ_day_scen_Hydro(jsel) = Reserve_educ_day_Hydro;

% Hydro -- Reserve Cost
Reserve_cost_educ_Hydro = [];
for j = 1:nt
   Reserve_cost_educ_Hydro(j,1) = sum(Reserve_cost_educ(Index_hydro,j));
end  

Reserve_cost_educ_day_Hydro = sum(Reserve_cost_educ_Hydro);
Reserve_cost_educ_day_scen_Hydro(jsel) = Reserve_cost_educ_day_Hydro;


% Other -- Total and Hourly Load
Disp_educ_Other = [];
for j = 1:nt
   Disp_educ_Other(j,1) = sum(Disp_educ(Index_other,j));
end  

Disp_educ_day_Other = sum(Disp_educ_Other);
for j = 1:nt
    Disp_educ_hourly_scen_Other(jsel,j) = Disp_educ_Other(j);
end

% Other -- Cost
Disp_educ_cost_Other = [];
for j = 1:nt
   Disp_educ_cost_Other(j,1) = sum(Disp_educ_cost(Index_other,j));
end  

Disp_educ_cost_day_Other = sum(Disp_educ_cost_Other);
Disp_educ_cost_day_scen_Other(jsel) = Disp_educ_cost_day_Other;

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

Reserve_educ_day_Other = sum(Reserve_educ_Other);
Reserve_educ_day_scen_Other(jsel) = Reserve_educ_day_Other;

% Other -- Reserve Cost
Reserve_cost_educ_Other = [];
for j = 1:nt
   Reserve_cost_educ_Other(j,1) = sum(Reserve_cost_educ(Index_other,j));
end  

Reserve_cost_educ_day_Other = sum(Reserve_cost_educ_Other);
Reserve_cost_educ_day_scen_Other(jsel) = Reserve_cost_educ_day_Other;


% Unknown -- Total and Hourly Load
Disp_educ_Unknown = [];
for j = 1:nt
   Disp_educ_Unknown(j,1) = sum(Disp_educ(Index_unknown,j));
end  

Disp_educ_day_Unknown = sum(Disp_educ_Unknown);
for j = 1:nt
    Disp_educ_hourly_scen_Unknown(jsel,j) = Disp_educ_Unknown(j);
end

% Unknown -- Cost
Disp_educ_cost_Unknown = [];
for j = 1:nt
   Disp_educ_cost_Unknown(j,1) = sum(Disp_educ_cost(Index_unknown,j));
end  

Disp_educ_cost_day_Unknown = sum(Disp_educ_cost_Unknown);
Disp_educ_cost_day_scen_Unknown(jsel) = Disp_educ_cost_day_Unknown;

% Unknown -- Avg Cost
Avg_educ_cost_Unknown = [];
for j = 1:nt
    Avg_educ_cost_Unknown(j,1) = Disp_educ_cost_Unknown(j,1)/ Disp_educ_Unknown(j,1);
end  

% Unknown -- Reserve Quantity
Reserve_educ_Unknown = [];
for j = 1:nt
   Reserve_educ_Unknown(j,1) = sum(Reserve_educ(Index_unknown,j));
end  

Reserve_educ_day_Unknown = sum(Reserve_educ_Unknown);
Reserve_educ_day_scen_Unknown(jsel) = Reserve_educ_day_Unknown;

% Unknown -- Reserve Cost
Reserve_cost_educ_Unknown = [];
for j = 1:nt
   Reserve_cost_educ_Unknown(j,1) = sum(Reserve_cost_educ(Index_unknown,j));
end  

Reserve_cost_educ_day_Unknown = sum(Reserve_cost_educ_Unknown);
Reserve_cost_educ_day_scen_Unknown(jsel) = Reserve_cost_educ_day_Unknown;

end

% Average of the day across Scenarios
Average_gen_day_scen = sum(Total_gen_day_scen(:))/nscen;
Average_demand_day_scen = sum(Total_demand_day_scen(:))/nscen;
Average_educ_day_Wind = sum(Disp_educ_day_Wind(:))/nscen;
Average_educ_day_Solar = sum(Disp_educ_day_Solar(:))/nscen;
Average_dcost_day_educ_scen = sum(Total_dcost_day_educ_scen(:))/nscen;
Average_rcost_day_educ_scen = sum(Total_rcost_day_educ_scen(:))/nscen;
Average_tcost_day_educ_scen = sum(Total_cost_day_educ_scen(:))/nscen;

% Average Demand across Scenarios
for j = 1:nt
    Average_demand_hourly_scen(1,j) = sum(Demand_hourly_scen(:,j))/nscen;
end

% Average Load across scenarios -- Nuclear 
for j = 1:nt
    Avg_educ_hourly_scen_Nuclear(j) = sum(Disp_educ_hourly_scen_Nuclear(:,j))/nscen;
end

% Average Load across scenarios -- Ng 
for j = 1:nt
    Avg_educ_hourly_scen_Ng(j) = sum(Disp_educ_hourly_scen_Ng(:,j))/nscen;
end

% Average Load across scenarios -- Dfo 
for j = 1:nt
    Avg_educ_hourly_scen_Dfo(j) = sum(Disp_educ_hourly_scen_Dfo(:,j))/nscen;
end

% Average Load across scenarios -- Rfo 
for j = 1:nt
    Avg_educ_hourly_scen_Rfo(j) = sum(Disp_educ_hourly_scen_Rfo(:,j))/nscen;
end

% Average Load across scenarios -- Coal 
for j = 1:nt
    Avg_educ_hourly_scen_Coal(j) = sum(Disp_educ_hourly_scen_Coal(:,j))/nscen;
end

% Average Load across scenarios -- Wood 
for j = 1:nt
    Avg_educ_hourly_scen_Wood(j) = sum(Disp_educ_hourly_scen_Wood(:,j))/nscen;
end

% Average Load across scenarios -- Hydro 
for j = 1:nt
    Avg_educ_hourly_scen_Hydro(j) = sum(Disp_educ_hourly_scen_Hydro(:,j))/nscen;
end

% Average Load across scenarios -- Wind
for j = 1:nt
    Avg_educ_hourly_scen_Wind(j) = sum(Disp_educ_hourly_scen_Wind(:,j))/nscen;
end

% Average Load across scenarios -- Solar 
for j = 1:nt
    Avg_educ_hourly_scen_Solar(j) = sum(Disp_educ_hourly_scen_Solar(:,j))/nscen;
end

% Average Load across scenarios -- Other 
for j = 1:nt
    Avg_educ_hourly_scen_Other(j) = sum(Disp_educ_hourly_scen_Other(:,j))/nscen;
end

% Average Load across scenarios -- Unknown 
for j = 1:nt
    Avg_educ_hourly_scen_Unknown(j) = sum(Disp_educ_hourly_scen_Unknown(:,j))/nscen;
end

% Average Net Demand across Scenarios
for j = 1:nt
    Average_net_demand_hourly_scen(1,j) = Average_demand_hourly_scen(1,j)- (Avg_educ_hourly_scen_Wind(j)...
        + Avg_educ_hourly_scen_Solar(j));
end

% Writing Excel Files
cd('Output_NY22_Aug2019_mod')
currentFolder = sprintf('Excel_educ_average_%2.2i',ctn);
mkdir(currentFolder)
cd(currentFolder)
% Average Load Table across Scenarios -- Fuel Type
A = [transpose(Avg_educ_hourly_scen_Nuclear) transpose(Avg_educ_hourly_scen_Ng) transpose(Avg_educ_hourly_scen_Dfo)...
    transpose(Avg_educ_hourly_scen_Rfo) transpose(Avg_educ_hourly_scen_Hydro) transpose(Avg_educ_hourly_scen_Coal)...
    transpose(Avg_educ_hourly_scen_Wood) transpose(Avg_educ_hourly_scen_Wind) transpose(Avg_educ_hourly_scen_Solar)...
    transpose(Avg_educ_hourly_scen_Other) transpose(Avg_educ_hourly_scen_Unknown)];

T_ft = array2table(A, 'VariableNames', {'Nuclear', 'NG', 'DFO', 'RFO',  'Hydro','Coal', 'Wood', 'Wind', ...
    'Solar', 'Other', 'Unknown'});
writetable(T_ft, 'Disp_educ_fr_fueltype.csv');

% Average Demand Table across Scenarios
D = [transpose(Average_demand_hourly_scen) transpose(Avg_educ_hourly_scen_Wind) transpose(Avg_educ_hourly_scen_Solar)...
    transpose(Average_net_demand_hourly_scen)];

D_title = array2table(D, 'VariableNames', {'Demand_avg', 'Wind_avg', 'Solar_avg', 'Net Demand_avg'});
writetable(D_title, 'Hourly_demand.csv');

cd ..
cd ..

end

cd('Output_NY22_Aug2019_mod')
%currentFolder = sprintf('Excel_educ_actual_%2.2i',ctn);
%currentFolder = sprintf('Excel_educ_max_%2.2i',ctn);
currentFolder = sprintf('Excel_educ_average_%2.2i',ctn);
mkdir(currentFolder)
cd(currentFolder)
 %Weekly Cost Table 
C_week = transpose([Average_gen_day_scen; Average_demand_day_scen; Average_educ_day_Wind; Average_educ_day_Solar;...
         Average_dcost_day_educ_scen; Average_rcost_day_educ_scen; Average_tcost_day_educ_scen]);
C_week = array2table(C_week,'VariableNames', {'Avg Gen(MW)', 'Avg Demand(MW)','Wind Gen(MW)','Solar Gen(MW)'...
                    'Disp Cost($)','Res Cost($)','Avg Cost($)'}, 'RowNames', {'Day1'});
writetable(C_week, 'Overall_results_week.csv','WriteRowNames',1);

