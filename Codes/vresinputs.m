% calculate transition matrices
% currently using percentage, move to check whether percentage or absolute amount is being used
% 2022.01.15
% Alberto J. Lamadrid L. & Arnab Sur

%addpath('~/gdrive/tex_files/fund_arpae19sh/data_code/scuc_data/infeasible_runs/');
%input_dir = '~/gdrive/tex_files/fund_arpae19sh/data_code/scuc_data/infeasible_runs/';
%addpath(input_dir);

% location to save information
%output_dir = '~/gdrive/abscores_data/profiles/NYISO/';						% all outputs
%addpath(output_dir);
%fn_prefix = output_dir;	      % functions

savefilep = {'vresinputs_ny22sh'};

dprefm = '2019_08';
define_constants
percentile = [0.15, 0.45, 1];				% percentiles for classification into bins
nj = 3;															% number of states desired

% load profiles
nsim = 50;													% number of simulations/scenarios created
ctn = 02;
jsel = 1;
%addpath(sprintf('%sload_profiles/', input_dir))
%addpath(sprintf('%swind_profiles/', input_dir))

% Demand profiles
% load info for some dimension checks
nameprofd = sprintf('profd%s_%2.2i_%d', dprefm, ctn, jsel);
profilesd = getprofiles(nameprofd);
profinfod = zeros(nsim, size(profilesd.values, 1), size(profilesd.values, 3));
nt = size(profilesd.values, 1);							% default : 24

nameprofv = sprintf('profw%s_%2.2i_%d', dprefm, ctn, jsel);
profilesv = getprofiles(nameprofv);
profinfov = zeros(nsim, size(profilesv.values, 1), size(profilesv.values, 3));

% inputs (nsim x nt x nz) number of scenarios, time periods, zones
for jsel = 1:nsim
	nameprofd = sprintf('profd%s_%2.2i_%d', dprefm, ctn, jsel);
	profilesd = getprofiles(nameprofd);
	prfdt = squeeze(profilesd.values);
	profinfod(jsel, :, :) = prfdt;

	nameprofv = sprintf('profw%s_%2.2i_%d', dprefm, ctn, jsel);
	profilesv = getprofiles(nameprofv);
	prfvt = squeeze(profilesv.values);
	profinfov(jsel, :, :) = prfvt;
	
end

% sorting criteria
%             [d] emand
%             [v] RES
%             [n] et load
%             if vres is desired sorting criteria, pass
%             [0]
%             [1]
%             [0]
sortc = [0;1;0];


% calculate total over the system
profd_sum = sum(profinfod, 3);
profv_sum = sum(profinfov, 3);
profn_sum = profd_sum - profv_sum;

% calculate ranges for series
profd_range = max(profd_sum, [], 1) - min(profd_sum, [], 1);
profv_range = max(profv_sum, [], 1) - min(profv_sum, [], 1);
profn_range = max(profn_sum, [], 1) - min(profn_sum, [], 1);

% sorting criteria, have by default the max range, unless something else specified
if max(sortc, [], 1) == 1;    % do we need to assign a value here??   
	sortc0 = sortc;
else
	sortc0 = [max(profd_range, [], 2); max(profv_range, [], 2); max(profn_range, [], 2)];
end

% add ranking keys to observations form scenarios
[marange, maxloc] = max(sortc0, [], 1);
if maxloc == 1												% max due to demand
	sortcf = [profd_sum, [1:1:nsim]'];
elseif maxloc == 2									% max due to VRES
	sortcf = [profv_sum, [1:1:nsim]'];
else																% max  due to net load
	sortcf = [profn_sum, [1:1:nsim]'];
end

% sort by each time period, save rankings from sorting 
% sortcflist [nsim x nt+1 x nt])
% Last dimension corresponds to the order observed when period t is used as sorting criteria (1:nsim)

% sortkey [nsim x nt]
% Sortkey stores the order (1:nsim) in the last dimension of sortcflist 

sortcflist = [];
sortkey = [];
for ct = 1:nt
	sortcflist(:, :, ct) = sortrows(sortcf, ct);
	sortkey(:, ct) = sortcflist(:, nt+1, ct);
end

cdfp = [0, percentile]' * nsim;

% determine bins to use according to break desired and target number of states (nj)
% bins [nsim x 1]
bins=[];
for ct = 1:nj
	lbi = floor(cdfp(ct+1)-cdfp(ct));
	bins = [bins; ct * ones(lbi, 1)];
end
defc= nsim - size(bins, 1);					% safety valve in case missing bins
bins = [bins; max(bins)*ones(defc, 1)];
 
% create new variable with sorted key and bin information in alternative manner interparsing the bins used
% sortkeybin [nsim, nt*2] (e..g., 1000 x 48)
sortkeybin = zeros(nsim, nt*2);
for ct = 1:nt
	sortkeybin(:, (ct-1)*2+1) = sortkey(:, ct);
	sortkeybin(:, ct*2) = bins;
end 

% compute transitions of bins from hour to hour
% sortedbin has two digit transition information for [nsim x (nt-1)]
sortedbin=[];
for ct2 = 1:size(sortkeybin, 2)/2 - 1
	for ct1 = 1:nsim
		sortedbin(ct1, ct2) = str2num([	num2str(sortkeybin(ct1, ct2*2)) ...
			num2str(sortkeybin(find(sortkeybin(:, ct2*2+1) == sortkeybin(ct1, (ct2-1)*2+1)), ...
			(ct2 + 1)*2))]);
	end
end

% define all potential transitions [nj x nj]
tc = [];
for ct1 =1:nj
	for ct2 = 1:nj
		tc(ct1, ct2) = str2num([num2str(ct1), num2str(ct2)]);
	end
end

% transcount computes how many observations are corresponding to each transition [nj x nj x nt]
transcount = [];
for ct3 = 1:size(sortedbin, 2)
	for ct1 = 1:nj
		for ct2 =1:nj
			transcount(ct1, ct2, ct3) = length(find(sortedbin(:, ct3)==tc(ct1, ct2)));
		end
	end
end

% transmat calculates the transition probabilities based on transcount, by making sum of row equal to 1 (nj-1 x nj-1 x nt-1) (markovian)
transmat = [];
for ct1 = 1:size(sortedbin, 2)
	for ct2 = 1:nj
		transmat(:, ct2, ct1) = transcount(:, ct2, ct1)./sum(transcount(:, :, ct1), 2);
	end
end

transmat2d = zeros(size(sortedbin, 2)*(nj+1), nj);
for ct1 = 1:size(sortedbin, 2)
	transmat2d((ct1-1)*(nj+1)+1:ct1*(nj+1)-1, :) = transmat(:, :, ct1);
end


% compute mean centroids for average VRES (nt x nj)

avgcentroids=[];
temp=[];
for ct1=1:nt
	temp=sortcflist(:, ct1,  ct1);
	for ct2=1:nj
	  lbi = floor(cdfp(ct2)+1:cdfp(ct2+1));
		avgcentroids(ct1, ct2) = mean(temp(lbi));
	end
end	

% centroids for VRES, compute mean centroids for each site (nt x nj x nzv)

nzv = size(profinfov, 3);

vres_scenarios = [];
temp=[];
% apply centroids of total wind back to each site wind
for ct2 = 1:nzv
	for ct1 = 1:nt
		for ct3 = 1:nj
		  lbi = floor(cdfp(ct3)+1:cdfp(ct3+1));
			vres_scenarios(ct1, ct3, ct2) = mean(profinfov(sortkey(lbi, ct1), ct1, ct2));		
		end
	end	
end	


% centroids for demand, compute mean centroids for each site (nt x nj x nzd)

nzd = size(profinfod, 3);

demand_scenarios = [];
temp=[];
% apply centroids of total demand to each zone
for ct2 = 1:nzd
	for ct1 = 1:nt
		for ct3 = 1:nj
		  lbi = floor(cdfp(ct3)+1:cdfp(ct3+1));
			demand_scenarios(ct1, ct3, ct2) = mean(profinfod(sortkey(lbi, ct1), ct1, ct2));		
		end
	end	
end	

save(sprintf('%s_%s_%2.2i', savefilep{1}, dprefm, ctn), 'transmat', 'demand_scenarios', 'vres_scenarios');

% write .m files to use in optimization/simulation
datei = '2022.04.25';
%% -- transtion prob matrix
for t = 1:size(transmat, 3)
  transraw(:, :, t) = transpose(transmat(:, :, t));
end
transrawa = zeros(size(transraw));
transrawa(:,:,2:nt) = transraw(:, :, 1:nt-1);
transrawa(:,:,1) = transraw(:, :, 1);
namefile = sprintf('transmat%s_%2.2i', dprefm, ctn);
settransmatsh(transrawa, datei, namefile);

%% -- wind scenarios
dp = vres_scenarios;
% datei same as transmat
namefile = sprintf('profw%s_%2.2i', dprefm, ctn);
startj = 1;
hourb = 0;
typt = 'w';
nzs = size(dp, 3);
col = PMAX;
chgtype = CT_REL;
table = CT_TGEN;
%namefilepath = output_dir;
setprofinfosh(dp, datei, namefile, startj, hourb, typt, nzs)

% THIS PART WILL NOT WORK YET, BUT ADDED HERE AS EXAMPLE FOR THE CASE OF DEMAND
%% -- Demand scenarios
dp = demand_scenarios;
% datei same as transmat
namefile = sprintf('profd%s_%2.2i', dprefm, ctn);
startj = 1;
hourb = 0;
typt = 'l';
nzs = size(dp, 3);
col = CT_LOAD_ALL_PQ;
chgtype = CT_REP;
table = CT_TAREALOAD;
%namefilepath = output_dir;
setprofinfosh(dp, datei, namefile, startj, hourb, typt, nzs, col, chgtype, table);