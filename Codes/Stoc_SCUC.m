% Stochastic SCUC
% Arnab Sur
% April, 2022




mkdir('Stochastic_SCUC_Aug');
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
%mpcase_name = {'NYISO03082022v2d'};	% case to run
%mpcase_name = {'NYISO03082022v2f'};
%mpcase_name = {'NYISO03082022v2d_revised'};
mpcase_name = {'NYISO03082022v2d_revised2'};

nshifts = nts;											% number of runs

tev = zeros(nshifts, 1);						% time elapsed vector

% upload a profile to get dimensions of run (nt)
ctn = dstart;
%jsel = 1;
nameprofd = sprintf('profd%s_%2.2i', dprefm, ctn);
profiletr = getprofiles(nameprofd);

% information that does not depend on case
nt = size(profiletr.values, 1);			% typically 24, sometimes 28
nj = size(profiletr.values, 2);

% -- transmat --
%nametrans = sprintf('transmat%s_%2.2i', dprefm, ctn);
%transmat = nametrans;


contab = ex_contab();

diary([savefilet '.txt']);
tic
%jsel = 1;														% scenario selected, 1 is base, cycle over these
for ctn = dstart:dstart+nts-1
    % transmat
    nametrans = sprintf('transmat%s_%2.2i', dprefm, ctn);
    transmat = nametrans;
	
    fprintf('==================================\n')
	fprintf('Period: %s-%2.2i\n', dprefm, ctn);
	fprintf('==================================\n')
	% case information
	mpcase = loadcase(mpcase_name{1});
	mpcase = addfuelname(mpcase);			% add idxs for fuels	
 	xgd = loadxgendata('xgd_ny22_ucsh', mpcase);

	% demand profiles
	nameprofd = sprintf('profd%s_%2.2i', dprefm, ctn);
	profiles = getprofiles(nameprofd);
	
	% wind profiles
	nameprofw = sprintf('profw%s_%2.2i', dprefm, ctn);
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
        for j = 1:nj
            mdi.FixedReserves(t, j, 1) = reserves;             % include fixed zonal reserves
        end
    end
 	
 	% run case with UC
 	mpoptuc = mpoption(mpopt, 'most.uc.run', 1);
	%mdi.UC.CommitKey(find(warmset.mdouc.UC.CommitSched == 0)) = 0;
	%mdi.UC.CommitKey(find(warmset.mdouc.UC.CommitSched == 1 )) = 2;    
	mdouc = most(mdi, mpoptuc);

 	tev(ctn, 1) = toc;

    cd('Stochastic_SCUC_Aug')
 	save(sprintf('%s_%s_%2.2i(dispatchable_revised2)', savefilepuc{1}, dprefm, ctn), 'mdouc');
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