%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'HybridPC';
M_.dynare_version = '4.5.4';
oo_.dynare_version = '4.5.4';
options_.dynare_version = '4.5.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('HybridPC.log');
M_.exo_names = 'eps_nu';
M_.exo_names_tex = 'eps\_nu';
M_.exo_names_long = 'eps_nu';
M_.exo_names = char(M_.exo_names, 'eps_cp');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_cp');
M_.exo_names_long = char(M_.exo_names_long, 'eps_cp');
M_.endo_names = 'ppi';
M_.endo_names_tex = 'ppi';
M_.endo_names_long = 'ppi';
M_.endo_names = char(M_.endo_names, 'y_gap');
M_.endo_names_tex = char(M_.endo_names_tex, 'y\_gap');
M_.endo_names_long = char(M_.endo_names_long, 'y_gap');
M_.endo_names = char(M_.endo_names, 'mpshock');
M_.endo_names_tex = char(M_.endo_names_tex, 'mpshock');
M_.endo_names_long = char(M_.endo_names_long, 'mpshock');
M_.endo_names = char(M_.endo_names, 'cp');
M_.endo_names_tex = char(M_.endo_names_tex, 'cp');
M_.endo_names_long = char(M_.endo_names_long, 'cp');
M_.endo_names = char(M_.endo_names, 'mp');
M_.endo_names_tex = char(M_.endo_names_tex, 'mp');
M_.endo_names_long = char(M_.endo_names_long, 'mp');
M_.endo_names = char(M_.endo_names, 'AUX_ENDO_LAG_1_1');
M_.endo_names_tex = char(M_.endo_names_tex, 'AUX\_ENDO\_LAG\_1\_1');
M_.endo_names_long = char(M_.endo_names_long, 'AUX_ENDO_LAG_1_1');
M_.endo_partitions = struct();
M_.param_names = 'aalpha';
M_.param_names_tex = 'aalpha';
M_.param_names_long = 'aalpha';
M_.param_names = char(M_.param_names, 'betta');
M_.param_names_tex = char(M_.param_names_tex, 'betta');
M_.param_names_long = char(M_.param_names_long, 'betta');
M_.param_names = char(M_.param_names, 'ggamma_b');
M_.param_names_tex = char(M_.param_names_tex, 'ggamma\_b');
M_.param_names_long = char(M_.param_names_long, 'ggamma_b');
M_.param_names = char(M_.param_names, 'pphi1');
M_.param_names_tex = char(M_.param_names_tex, 'pphi1');
M_.param_names_long = char(M_.param_names_long, 'pphi1');
M_.param_names = char(M_.param_names, 'pphi2');
M_.param_names_tex = char(M_.param_names_tex, 'pphi2');
M_.param_names_long = char(M_.param_names_long, 'pphi2');
M_.param_names = char(M_.param_names, 'kkappa');
M_.param_names_tex = char(M_.param_names_tex, 'kkappa');
M_.param_names_long = char(M_.param_names_long, 'kkappa');
M_.param_names = char(M_.param_names, 'rho_cp');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_cp');
M_.param_names_long = char(M_.param_names_long, 'rho_cp');
M_.param_names = char(M_.param_names, 'rho_mp');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_mp');
M_.param_names_long = char(M_.param_names_long, 'rho_mp');
M_.param_names = char(M_.param_names, 'dSmp');
M_.param_names_tex = char(M_.param_names_tex, 'dSmp');
M_.param_names_long = char(M_.param_names_long, 'dSmp');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 6;
M_.param_nbr = 9;
M_.orig_endo_nbr = 5;
M_.aux_vars(1).endo_index = 6;
M_.aux_vars(1).type = 1;
M_.aux_vars(1).orig_index = 2;
M_.aux_vars(1).orig_lead_lag = -1;
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('HybridPC_static');
erase_compiled_function('HybridPC_dynamic');
M_.orig_eq_nbr = 5;
M_.eq_nbr = 6;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 1 6 12;
 2 7 0;
 0 8 0;
 3 9 0;
 4 10 0;
 5 11 0;]';
M_.nstatic = 1;
M_.nfwrd   = 0;
M_.npred   = 4;
M_.nboth   = 1;
M_.nsfwrd   = 1;
M_.nspred   = 5;
M_.ndynamic   = 5;
M_.equations_tags = {
  1 , 'name' , 'New Keynesian Phillips Curve eq. (22)' ;
  2 , 'name' , 'Dynamic IS Curve eq. (23)' ;
  3 , 'name' , 'cost push shock' ;
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(6, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(9, 1);
M_.NNZDerivatives = [20; -1; -1];
M_.params( 1 ) = (-1);
aalpha = M_.params( 1 );
M_.params( 2 ) = 0.3;
betta = M_.params( 2 );
M_.params( 3 ) = 0.6;
ggamma_b = M_.params( 3 );
M_.params( 6 ) = 0.4;
kkappa = M_.params( 6 );
M_.params( 4 ) = 1.2;
pphi1 = M_.params( 4 );
M_.params( 5 ) = (-0.4);
pphi2 = M_.params( 5 );
M_.params( 7 ) = 0;
rho_cp = M_.params( 7 );
M_.params( 8 ) = 0;
rho_mp = M_.params( 8 );
M_.params( 9 ) = 0.1;
dSmp = M_.params( 9 );
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 1;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 0;
options_.nocorr = 1;
options_.nomoments = 1;
options_.order = 1;
options_.periods = 1000;
var_list_ = char();
info = stoch_simul(var_list_);
save('HybridPC_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('HybridPC_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('HybridPC_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('HybridPC_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('HybridPC_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('HybridPC_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('HybridPC_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
