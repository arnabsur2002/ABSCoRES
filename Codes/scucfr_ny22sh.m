% Prototype for SCUC run different days
% 2022.03.11
% Alberto J. Lamadrid L. & Arnab Sur 

%addpath('~/gdrive/tex_files/fund_arpae19sh/data_code/scuc/');
%input_dir = '~/gdrive/abscores_data/profiles/NYISO/';
%addpath(input_dir);

% location to save information
%fn_prefix = '~/gdrive/tex_files/fund_arpae19sh/data_LLNL/';		% functions
%output_dir = '~/gdrive/abscores_data/runs/NYISO/';						% all outputs
%addpath(output_dir);

% warm start for uc, standardize this
%warmset = load('output_ny_uc_2019_08_02');
%warmset = load('scucfr_ny22sh_uc_2019_08_02_01');
%cd('Decision_Matrix_Aug')
%warmset = load('scucfr_ny22sh_2019_08_02_01');
%cd ..
%dprefix = '/Users/ajlamadrid/gdrive/tex_files/fund_arpae19sh/data_code/scuc/';  % local
%dprefix = output_dir;  % local

%if ~exist(output_dir, 'dir')
    mkdir('SCUC_Aug');
%end

%nts = 7;														% number of time cycles
nts = 1;														% number of time cycles
dstart = 2;													% day started

define_constants

dprefm = '2019_08';									% matlab creates error when using '-' instead of '_'

% define options runs
mpopt = mpoption('verbose', 2, 'out.all', 1);
mpopt = mpoption(mpopt, 'opf.violation', 1e-2);
mpopt = mpoption(mpopt, 'most.dc_model', 1, 'most.security_constraints', 0);
mpopt = mpoption(mpopt, 'most.solver', 'GUROBI');
%mpopt = mpoption(mpopt, 'most.skip_prices', 1);
mpopt = mpoption(mpopt, 'most.uc.run', 0);

savefilet = 'scucfr_ny22sh';				% name of text file for tracking
savefilepuc = {'scucfr_ny22sh'};			% array of names to save cases
savefilep = {'scedfr_ny22sh'};	% array of names to save cases
%mpcase_name = {'NYISO03082022v2f'};	% case to run
%mpcase_name = {'NYISO03082022v2d'};
%mpcase_name = {'NYISO03082022v2d15'};
mpcase_name = {'NYISO03082022v2d_revised'};
nshifts = nts;											% number of runs

tev = zeros(nshifts, 1);						% time elapsed vector

% upload a profile to get dimensions of run (nt)
ctn = dstart;
jsel = 1;
nameprofd = sprintf('profd%s_%2.2i_%d', dprefm, ctn, jsel);
profiletr = getprofiles(nameprofd);

% information that does not depend on case
nt = size(profiletr.values, 1);			% typically 24, sometimes 28
% -- transmat --
for t = 1:nt
  transmat{1, t} = 1;              	% create deterministic transition probability matrix
end

contab = ex_contab();

diary([savefilet '.txt']);
tic
jsel = 1;														% scenario selected, 1 is base, cycle over these
for ctn = dstart:dstart+nts-1
	fprintf('==================================\n')
	fprintf('Period: %s-%2.2i\n', dprefm, ctn);
	fprintf('==================================\n')
	% case information
	mpcase = loadcase(mpcase_name{1});
	mpcase = addfuelname(mpcase);			% add idxs for fuels	
 	xgd = loadxgendata('xgd_ny22_ucsh', mpcase);

	% demand profiles
%	nameprofd = sprintf('profd%s_%2.2i', dprefm, ctn);
	nameprofd = sprintf('profd%s_%2.2i_%d', dprefm, ctn, jsel);
	profiles = getprofiles(nameprofd);
	
	% wind profiles
%	nameprofw = sprintf('profw%s_%2.2i', dprefm, ctn);
	nameprofw = sprintf('profw%s_%2.2i_%d', dprefm, ctn, jsel);
 	profiles = getprofiles(nameprofw, profiles, mpcase.iwind);

 	% setup info for SCED runs
 	mdi = loadmd(mpcase, transmat, xgd, [], contab, profiles);
 	
 	%% creation of fixed reserve requirements
	reserves = fixreq_NewYork_bau22f;
% reserves = fixreqscucf(mpcase);
	reserves.qty(mpcase.iwind) = 0;
	rfact = 1;											% debug, reserve factor to support convergence, aim for unit multiplier (1)
	reserves.req = rfact * reserves.req;

	for t = 1:nt
		mdi.FixedReserves(t, 1, 1) = reserves; % include fixed zonal reserves
	end
 	
 	% run case with UC
 	mpoptuc = mpoption(mpopt, 'most.uc.run', 1);
	%mdi.UC.CommitKey(find(warmset.mdouc.UC.CommitSched == 0)) = 0;
	%mdi.UC.CommitKey(find(warmset.mdouc.UC.CommitSched == 1 )) = 2;    
	mdouc = most(mdi, mpoptuc);

 	tev(ctn, 1) = toc;
    
    cd('SCUC_Aug')
 	save(sprintf('%s_%s_%2.2i_%2.2i(dispatchable_revised)', savefilepuc{1}, dprefm, ctn, jsel), 'mdouc');
 	cd ..
    
 	% setup info for SCED runs
 	mdiuc = mdi;
 	clear('mdi');
 	mdi = loadmd(mpcase, transmat, xgd, [], contab, profiles);
 	
 	%% modifying fixed reserve requirements
 	reservesd = reserves;
	rfact = .25;											% reserve factor, SCED runs aim for unit multiplier (1)
	reservesd.req = rfact * reservesd.req;
	
	for t = 1:nt
		mdi.FixedReserves(t, 1, 1) = reservesd; % include fixed zonal reserves
	end
 	
 	%mdi.UC.CommitKey = mdouc.UC.CommitKey;
% 	mdi.UC.CommitSched = mdouc.UC.CommitSched;
  mdi.UC.CommitSched(find(mdouc.UC.CommitSched == 0)) = 0 ;
	mdi.UC.CommitSched(find(mdouc.UC.CommitSched == 1)) = 1 ;
 	 	
 	mdo = most(mdi, mpopt);
 	tev(ctn, 1) = toc;
 	
    cd('SCUC_Aug')
 	save(sprintf('%s_%s_%2.2i_%2.2i(dispatchable_revised)', savefilep{1}, dprefm, ctn, jsel), 'mdo');    % save structure with results
 	cd ..
    
 	clearvars('mpcase', 'xgd', 'iwind', 'iess', 'mdi', 'mdo');
 	fprintf('----------------------------------\n')
 	fprintf('Time Run: %d\n', tev(ctn, 1));
 	fprintf('----------------------------------\n')
 	
	% calculation economic variables 	
%   optd.savename = sprintf('%s-data_c%s_%2.2i', savefilep{1}, dprefm, ctn);
%   optd.savepath = output_dir;
%   optd.tframe = nt;
%   nOut = 90;
%   rev = load(sprintf('%s_%s_%2.2i_%2.2i', savefilep{1}, dprefm, ctn, jsel));
%   [outArgs{1:nOut}] = data_mpsd(rev.mdo, optd);
%     
% 	nOut = 123;
% 	[outArgs{1:nOut}] = data_fxres(rev.mdo);

%	checkfroutsh

	% plots
%   optp.savefile = sprintf('%s-plotsr%s_%2.2i', savefilep{1}, dprefm, ctn);
%   savedvars = optd.savename;
%   optp.savepath = output_dir;
%   optp.profiles = [];
%   doplots_hpdcf(rev.mdo, savedvars, optp);

 	
end	

fprintf('==================================\n')
fprintf('Time Run Total: %d\n', tev);
fprintf('==================================\n')
diary off;