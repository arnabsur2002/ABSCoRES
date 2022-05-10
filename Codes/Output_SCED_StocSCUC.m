% Output: Multiple ED with Stochastic UC Decisions (with Reserves)
% Arnab Sur
% May 5, 2022

mkdir('Output_Stoc_Aug2019');
define_constants

savefilepscen = {'ed_scen_ny22_fr'};		
savefilep = {'ed_ny22_fr'};

dstart = 02;                                   % start date
nts = 1;                                       % # of days 
nscen = 50;                                    % # of scenarios

for ctn = dstart:dstart+nts-1
% Loading Data
for jsel = [2:51];
data_ed_scen = sprintf('scedfr_ny22sh_2019_08_%2.2i_%2.2i', ctn, jsel);
load(data_ed_scen);


data_uc = sprintf('scucfr_ny22sh_2019_08_%2.2i(dispatchable)',ctn);
load(data_uc);

nt = mdo.idx.nt;   
ngen = 605;

% Dispatchable Load 
load_disp = [];
for t = 1:nt
    load_disp(1,t) = sum(total_load(mdo.flow(t).mpc, 'all', 'DISPATCHABLE'));
end
load_disp_day(jsel-1,:) = load_disp;

% Fixed Load
load_fixed = [];
for t = 1:nt
    load_fixed(1,t) = sum(total_load(mdo.flow(t).mpc, 'all', 'FIXED'));
end

load_fixed_day(jsel-1,:) = load_fixed;

% Sum of Dispatchable and Fixed Load
load_both = [];
for t = 1:nt
    load_both(1,t) = sum(total_load(mdo.flow(t).mpc, 'all', 'BOTH'));
end

load_both_day(jsel-1,:) = load_both;

% Generation Amount
generation = [];
for t = 1:nt
    generation(1,t) = sum(mdo.results.ExpectedDispatch(1:ngen,t));
end 

generation_day(jsel-1,:) = generation;

% Amount of Lost Load
shed_load = [];
for t = 1:nt
    shed_load(1,t) = sum(loadshed(mdo.flow(t).mpc.gen));
end 

shed_load_day(jsel-1,:) = shed_load;


% For the Verification: qty1 = load_both
qty1 = [];
for t = 1:nt
    qty1(1,t) = generation(1,t) + shed_load(1,t);
end

% Net Dispatchable Load
dispatchable_net = [];
for t = 1:nt
    dispatchable_net(1,t) = sum(mdo.results.ExpectedDispatch(ngen+1:end,t));
end

dispatchable_net_day(jsel-1,:) = dispatchable_net;

% For the Verification: qty2 = - dispatchable_net
qty2 = [];
for t = 1:nt
    qty2(1,t) = load_disp(1,t) - shed_load(1,t);
end


% Financial Quantities: welfare_net = sum(dispatchable_net)*9000 + sum(generation_cost)

% Net Welfare
welfare_net = [];
for t = 1:nt
    welfare_net(:, t) = totcost(mdo.mpc.gencost, mdo.results.ExpectedDispatch(:,t));
end 

welfare_net_day(jsel-1,:,:) = welfare_net;

% Generation Cost
generation_cost = [];
for t = 1:nt
    generation_cost(1, t) = sum(welfare_net(1:ngen,t));
end 

generation_cost_day(jsel-1,:) = generation_cost;

% Total Amount of Reserve Used
reserve_used = [];
for i = 1:ngen
    for t = 1:nt
        reserve_used(i,t) = mdouc.results.ExpectedDispatch(i,t) - mdo.results.ExpectedDispatch(i,t);
    end
end

reserve_used_day(jsel-1,:,:) = reserve_used;

end


% Average of the day across Scenarios
avg_load_disp = [];
for t = 1:nt
    avg_load_disp(1,t) = sum(load_disp_day(:,t))/nscen;
end 

avg_load_fixed = [];
for t = 1:nt
    avg_load_fixed(1,t) = sum(load_fixed_day(:,t))/nscen;
end 

avg_load_both = [];
for t = 1:nt
    avg_load_both(1,t) = sum(load_both_day(:,t))/nscen;
end 

avg_generation = [];
for t = 1:nt
    avg_generation(1,t) = sum(generation_day(:,t))/nscen;
end 

avg_shed_load = [];
for t = 1:nt
    avg_shed_load(1,t) = sum(shed_load_day(:,t))/nscen;
end 

avg_dispatchable_net = [];
for t = 1:nt
    avg_dispatchable_net(1,t) = sum(dispatchable_net_day(:,t))/nscen;
end 

avg_welfare_net = [];
for i = 1:size(welfare_net_day,2)
    for t = 1:nt
        avg_welfare_net(i,t) = sum(welfare_net_day(:,i,t))/nscen;
    end
end 

avg_generation_cost = [];
for t = 1:nt
    avg_generation_cost(1,t) = sum(generation_cost_day(:,t))/nscen;
end 


% Reserve Up/Down across Scenarios
reserve_up = [];
for i = 1:ngen
    for t = 1:nt
        if  max(reserve_used_day(:,i,t)) >= 0;
            reserve_up(i,t) = max(reserve_used_day(:,i,t));
        else reserve_up(i,t) = 0;
        end 
    end
end 

reserve_down = [];
for i = 1:ngen
    for t = 1:nt
        if  min(reserve_used_day(:,i,t)) < 0;
            reserve_down(i,t) = min(reserve_used_day(:,i,t));
        else reserve_down(i,t) = 0;
        end 
    end
end 

cd('Output_Stoc_Aug2019')
currentFolder = sprintf('Excel_educ_average_%2.2i',ctn);
mkdir(currentFolder)
cd(currentFolder)

% Average Load, Generation and Lost Load Table across Scenarios
D = [transpose(avg_load_disp) transpose(avg_load_fixed) transpose(avg_load_both) transpose(avg_generation)...
    transpose(avg_shed_load)];

D_title = array2table(D, 'VariableNames', {'Dispatchable_avg', 'Fixed_avg', 'Both_avg', 'Generation_avg', 'Loadshedding_avg'});
writetable(D_title, 'Hourly_Average.csv');

% Average Net Welfare
writematrix(transpose(avg_welfare_net),'Welfare_hourly.csv');

% Reserve Up/Down
writematrix(transpose(reserve_up),'Reserve_up.csv');
writematrix(transpose(reserve_down),'Reserve_down.csv');


cd ..
cd ..

end 
