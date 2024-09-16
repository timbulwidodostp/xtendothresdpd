
/*
This code is proto-type code for a STATA command, which estimates dynamic panel 
threshold regression with endogeneity, based on Seo & Shin (2016, JoE)

Notes:
1. xtset should be done before running this. Moreover, variables must be 
   sorted by (1) panel variable and (2) time variable beforehand.
2. Data should be strongly balanced panel.
3. Outputs are coefficients estimates and asymptotic S.E.
4. Inputs should be put as y q x, where q is the threshold variable.
5. 'moremata' library is required.
*/

/*
// Housekeeping
mata mata clear
capture program drop xthenreg
*/ //Not needed here

program define xthenreg, eclass
	syntax varlist(numeric) [if] [in] 									///
							[, ENDOgenous(string) EXOgenous (string)	/// 
							 inst(string) kink static 					///
							 grid_num(integer 100) trim_rate(real 0.4) 	/// 
							 h_0(real 1.5) boost(real 0)]
	quietly xtset
	marksample touse
	
	
	// Check if input data are strongly balanced panel.
	if ("`r(balanced)'" != "strongly balanced"){
	display as error "Error: The given data are not strongly balanced."
	exit
	}
	
	// Check if MOREMATA library is installed.
	mata t_vec = (0, 1, 2, 3)'
	capture mata mm_quantile(t_vec, 1, 0.5)
	
	if (_rc == 0) {
		mata mata drop t_vec
	}
	else {
	display as error "Error: Please install MOREMATA module first." _newline ///
	`"You can install the module with command "ssc install moremata"."'
	exit
	}
	
	
	// ** Model setting to estimate ** //
	if ("`endogenous'" == "") { 
	// Is there an endogenous variable among covariates? 
		local k_endo 0
		local var_all "`varlist'"
	}
	else {
		local k_endo: word count `endogenous'
		local var_all "`varlist' `endogenous'"
	}
	
	if ("`exogenous'" == "") { 
	// Is there an exogeneous variable among covariates? 
		local k_exo 0
	}
	else {
		local k_exo: word count `exogenous'
		local var_all "`var_all' `exogenous'"
	}
	
	if ("`kink'" == ""){ 
	// Is the model kink? (1 if yes)
		local flag_kink 0
	}
	else{
		local flag_kink 1
	}
	
	if ("`static'" == ""){	
	// Is the model static?
		local flag_static 0
	}
	else{
		local flag_static 1
	}
	
	if ("`inst'" == "") { 
	// Is there an additional IV?
		local flag_inst 0
		local inst_list .
	}
	else {
		local flag_inst 1
		local inst_list "`inst'"
	}
	
	if (`boost' == 0){
	// Should we do bootstrapping for linearity test?
		local flag_b 0
		local B 0
	}
	else {
		local flag_b 1
		local B `boost'
	}
	
	mata: GMM_estimates = .
	mata: cov_mat = .
	mata: boots_p = . // Wald stat. results
	mata: CI_95 = .; N = .; T = .
	mata: maincode_0719("`var_all'", `k_endo', `k_exo', `grid_num', `trim_rate', `h_0',  ///
						`flag_static', `flag_kink', `flag_inst', "`inst_list'", ///
						"`touse'", `flag_b', `B', GMM_estimates, cov_mat, 		///
						boots_p, CI_95, N, T)
	
	
	mata: st_matrix("b", GMM_estimates') // b should be (1*p). So include transpose operator.
	mata: st_matrix("V", cov_mat)
	mata: st_matrix("CI_95", CI_95)
	mata: st_numscalar("N", N); st_numscalar("T", T); st_numscalar("boots_p", boots_p)
	
	mata mata drop GMM_estimates cov_mat boots_p CI_95 N T
	// For housekeeping
	
	//// Arrange depvar, qx (06/25) ////
	local num_vl: word count `varlist' 
	// # of variables y, q, x
	
	if ("`static'" == ""){
		local indepvar = "L.y" // if dynamic, L.y is the first depvar
	}
	else {
		local indepvar
	}
	
	forvalues ii2 = 3/`num_vl' {
		local tmp_v : word `ii2' of `varlist'
		local indepvar = "`indepvar'" + " " + "`tmp_v'"
	}
	
	local depvar: word 1 of `varlist'
	local qx: word 2 of `varlist'
	// y and Threshold Variable
	
	
	
	if ("`kink'" == "" & "`static'" == "") {
	// Dynamic & not kink
		local coef_list = "Lag_y_b" // To label coefficients (the first one should be L.y)
		local num_v: word count `var_all' // # of variables including y, q, x, z
		
		forvalues ii = 3/`num_v' {
			local temp_var : word `ii' of `var_all'
			local coef_list = "`coef_list'" + " " + "`temp_var'" + "_b"
		}
		
		local coef_list = "`coef_list'" + " " + "cons_d" + " " + "Lag_y_d"
		
		forvalues ii = 3/`num_v' {
			local temp_var : word `ii' of `var_all'
			local coef_list = "`coef_list'" + " " + "`temp_var'" + "_d"
		}
		
		local coef_list = "`coef_list'" + " " + "r"
		matrix rownames V = `coef_list'
		matrix colnames V = `coef_list'
		matrix colnames b = `coef_list'
		
		
		// Store Results in e() //
		
		ereturn post b V, depname("`depvar'") 
		ereturn display
		
		matrix rownames CI_95 = `coef_list'
		matrix colnames CI_95 = "LB" "UB"
		ereturn matrix CI = CI_95
		
		ereturn scalar N = N
		ereturn scalar T = T
		ereturn scalar boots_p = boots_p
		ereturn scalar grid = `grid_num'
		ereturn scalar trim = `trim_rate'
		ereturn scalar bs = `B'
		
		ereturn local indepvar `indepvar'
		ereturn local depvar `depvar'
		ereturn local qx `qx'
		ereturn local zx `inst'
		
	}
	else if ("`kink'" == "") {
	// Static & not kink
		local coef_list // To label coefficients (No need to include L.y)
		local num_v: word count `var_all' // # of variables including y, q
		
		forvalues ii = 3/`num_v' {
			local temp_var : word `ii' of `var_all'
			local coef_list = "`coef_list'" + " " + "`temp_var'" + "_b"
		}
		
		local coef_list = "`coef_list'" + " " + "cons_d"
		
		forvalues ii = 3/`num_v' {
			local temp_var : word `ii' of `var_all'
			local coef_list = "`coef_list'" + " " + "`temp_var'" + "_d"
		}
		
		local coef_list = "`coef_list'" + " " + "r"
		matrix rownames V = `coef_list'
		matrix colnames V = `coef_list'
		matrix colnames b = `coef_list'
		
		
		// Store Results in e() //
		
		ereturn post b V, depname("`depvar'") 
		ereturn display
		
		matrix rownames CI_95 = `coef_list'
		matrix colnames CI_95 = "LB" "UB"
		ereturn matrix CI = CI_95
		
		ereturn scalar N = N
		ereturn scalar T = T
		ereturn scalar boots_p = boots_p
		ereturn scalar grid = `grid_num'
		ereturn scalar trim = `trim_rate'
		ereturn scalar bs = `B'
		
		ereturn local indepvar `indepvar'
		ereturn local depvar `depvar'
		ereturn local qx `qx'
		ereturn local zx `inst'	
	}
	
	else if ("`static'" == ""){
	// Dynamic kink case
		local coef_list = "Lag_y_b" // To label coefficients (the first one should be L.y)
		local num_v: word count `var_all' // # of variables including y, q
		
		forvalues ii = 3/`num_v' {
			local temp_var : word `ii' of `var_all'
			local coef_list = "`coef_list'" + " " + "`temp_var'" + "_b"
		}
		
		local coef_list = "`coef_list'" + " " + "kink_slope" 
		
		local coef_list = "`coef_list'" + " " + "r"
		matrix rownames V = `coef_list'
		matrix colnames V = `coef_list'
		matrix colnames b = `coef_list'
		

		// Store Results in e() //
		
		ereturn post b V, depname("`depvar'") 
		ereturn display
		
		matrix rownames CI_95 = `coef_list'
		matrix colnames CI_95 = "LB" "UB"
		ereturn matrix CI = CI_95
		
		ereturn scalar N = N
		ereturn scalar T = T
		ereturn scalar boots_p = boots_p
		ereturn scalar grid = `grid_num'
		ereturn scalar trim = `trim_rate'
		ereturn scalar bs = `B'
		
		ereturn local indepvar `indepvar'
		ereturn local depvar `depvar'
		ereturn local qx `qx'
		ereturn local zx `inst'
	}
	
	else{
	// Static kink case
		local coef_list // To label coefficients (No need to include L.y)
		local num_v: word count `var_all' // # of variables including y, q
		
		forvalues ii = 3/`num_v' {
			local temp_var : word `ii' of `var_all'
			local coef_list = "`coef_list'" + " " + "`temp_var'" + "_b"
		}
		
		local coef_list = "`coef_list'" + " " + "kink_slope" 
		
		local coef_list = "`coef_list'" + " " + "r"
		matrix rownames V = `coef_list'
		matrix colnames V = `coef_list'
		matrix colnames b = `coef_list'
		

		// Store Results in e() //
		
		ereturn post b V, depname("`depvar'") 
		ereturn display
		
		matrix rownames CI_95 = `coef_list'
		matrix colnames CI_95 = "LB" "UB"
		ereturn matrix CI = CI_95
		
		ereturn scalar N = N
		ereturn scalar T = T
		ereturn scalar boots_p = boots_p
		ereturn scalar grid = `grid_num'
		ereturn scalar trim = `trim_rate'
		ereturn scalar bs = `B'
		
		ereturn local indepvar `indepvar'
		ereturn local depvar `depvar'
		ereturn local qx `qx'
		ereturn local zx `inst'
	}
	
end

mata:

void maincode_0719(string scalar vname, real scalar k_endo,			///
					real scalar k_exo, real scalar grid_num, 		///
					real scalar trim_rate, real scalar h_0, 		///
					real scalar flag_static, real scalar flag_kink,	///
					real scalar flag_inst, string scalar inst_list, ///
					string scalar touse,							/// 
					real scalar flag_b, real scalar B,				///
					GMM_estimates, cov_mat, boots_p, CI_95, N, T)
{
// ***(Step 1): Data construction*** //

// #### (Start) Updated 18/02/16 ####
pan_var = st_global("r(panelvar)")
time_var = st_global("r(timevar)")
PT_mat = st_data(., pan_var + " " + time_var, touse)
// (# of touse-selected * 2) matrix.

N = .
T = .
N_ori = .
ind_non_ms = .
index_non_ms(N, T, PT_mat, N_ori, ind_non_ms)


n_del = strofreal(N_ori - N)
if (N_ori - N > 0) {
	display(n_del + " sample(s) are ignored further due to missing values")
}
// #### (End) Updated 18/02/16 ####


data_all = st_data(., vname, touse)
data_all = select(data_all, ind_non_ms)
//First col. = y, second col. = q, other col.'s = x
//"touse" restricts the sample selection

y_mat = data_all[., 1]
y_mat = colshape(y_mat, T)
// It is (N*T) matrix containing y_{i, t} 

q_mat = data_all[., 2]
q_mat = colshape(q_mat, T)
// It is (N*T) matrix containing q_{i, t} 

k_1 = cols(data_all) - 1
// # of independent variables (including y_{i, t-1})

x1_mat = (J(N, 1, 1), y_mat[., 1..T-1]) 
x_mat = x1_mat
//Note that error occurs if N is not integer. (If sample restriction was time-variant...etc.)


for (jj=1; jj<=(k_1-1); jj=jj+1) {
	//For-loop constructing (Nk_1 * T) x_mat(containing y_{i, t-1} and x_{i, t})

	x_temp = data_all[., 2+jj]
	x_temp = colshape(x_temp, T)
	x_mat = (x_mat \ x_temp)
}

if (flag_inst == 1) {
	data_IV = st_data(., inst_list, touse)
	data_IV = select(data_IV, ind_non_ms)
	
	k_inst = cols(data_IV)
	
	if (k_inst >= 2) {
	// When there exist 2 or more additional IV's
		
		IV_mat = colshape(data_IV[., 1], T)
		for (jj2=2; jj2<=k_inst; jj2=jj2+1) {
			//For-loop constructing (Nk_inst * T) IV_mat
			IV_temp = data_IV[., jj2]
			IV_temp = colshape(IV_temp, T)
			IV_mat = (IV_mat \ IV_temp)
		}
	}
	else {
	// When there exists only 1 IV.
		IV_mat = colshape(data_IV[., 1], T)
	}
}
else {
	flag_inst = flag_inst + 0
	IV_mat = .
	k_inst = 0
	//NO meaning. Just to fill out 'else' part
}


grid_th = .
grid_con(q_mat, trim_rate, grid_num, grid_th)
//grid_th = J(grid_num, 1, 0)

y_mat_fd = . // This will be (N*(T-1)) matrix
x_mat_fd = . // This will be (Nk_1*(T-1)) matrix
FD_con(y_mat, x_mat, T, y_mat_fd, x_mat_fd)


//#### 180719 updated part! (Start) ####//
k_endo_n = k_endo + k_exo 

moment_mat = .
moment_mat_con(T, k_1, k_endo_n, moment_mat)
mt_c = cols(moment_mat); mt_r = rows(moment_mat);
moment_mat = moment_mat + (J(mt_r, mt_c - k_exo, 0), J(mt_r, k_exo, 2))

Z_mat = .
Z_mat_con(y_mat, x_mat, x_mat_fd, moment_mat, k_endo_n, flag_static, Z_mat)

if (flag_inst == 1) {
	Z_mat_IV_add(Z_mat, moment_mat, IV_mat, k_inst, flag_static)
}
else{
	moment_mat = moment_mat :+ 0
	// NO meaning, just to fill out 'else' part
}
//#### 180719 updated part! (End) ####//


y_mat_rep = . 
x_temp1 = . 
x_temp2 = .
x_temp3 = .
q_temp1 = .
q_temp2 = .
//contains y_{it}, \delta x_{it}, x_{it}, x_{it-1}, q_{it}, q_{it-1} comparable with Z_mat
//Each should be (L*N), (L*Nk_1), (L*Nk_1), (L*Nk_1), (L*N), (L*N) matrix, respectively.

GMM_var_con(y_mat_fd, x_mat_fd, x_mat, q_mat, moment_mat, T, flag_static,	/// 
			y_mat_rep, x_temp1, x_temp2, x_temp3, q_temp1, q_temp2)


// ***(Step 2): 1st-step GMM estimate*** //
W_n = .

if (flag_static == 1) {
	GMM_W_n_static(Z_mat, moment_mat, T, W_n)
}
else {
	GMM_W_n_con(Z_mat, moment_mat, T, W_n)
}



g_1n_bar = rowsum(Z_mat:*y_mat_rep)/N 
//This is row-wise sum. Identical to g_1n_bar in Seo & Shin(2016)
 
J_result = J(grid_num, 1, 0)

if (flag_static == 1){
	end_p = cols(x_temp1);
	
	x_temp1_c = x_temp1[., N+1..end_p]
	x_temp2_c = x_temp2[., N+1..end_p]
	x_temp3_c = x_temp3[., N+1..end_p]
	// Copies of matrices, after removing columns that contain
	// y_{i, t-1} related terms
	
	// These matrices will be used if we assume static model
}


// Now we calculate 1st step GMM estimator w.r.t. 4 different cases.

if (flag_static == 1 & flag_kink == 0) {
// (Case 1): Static model, not kink design

	for (rr=1; rr<=grid_num; rr=rr+1) {

		r_th = grid_th[rr]
		
		g_2n_bar = .
		g_n_bar = .
		GMM_est = .

		GMM_cal(g_1n_bar, r_th, Z_mat, x_temp1_c, x_temp2_c, x_temp3_c, ///
				   q_temp1, q_temp2, W_n, g_2n_bar, g_n_bar, GMM_est)

		J_result[rr] = g_n_bar' * W_n * g_n_bar
		
	}
}

else if (flag_static == 0 & flag_kink == 1) {
// (Case 2): Dynamic model, with kink

	for (rr=1; rr<=grid_num; rr=rr+1) {
	
		r_th = grid_th[rr]
		
		g_2n_bar = .
		g_n_bar = .
		GMM_est = .
		
		GMM_cal_kink(g_1n_bar, r_th, Z_mat, x_temp1, x_temp2, x_temp3, ///
				q_temp1, q_temp2, W_n, g_2n_bar, g_n_bar, GMM_est)

		J_result[rr] = g_n_bar' * W_n * g_n_bar
		
	}
}

else if (flag_static == 1 & flag_kink == 1) {
// (Case 3): Static model, with kink 

	for (rr=1; rr<=grid_num; rr=rr+1) {

		r_th = grid_th[rr]
		
		g_2n_bar = .
		g_n_bar = .
		GMM_est = .
		GMM_cal_kink(g_1n_bar, r_th, Z_mat, x_temp1_c, x_temp2_c, x_temp3_c, ///
						q_temp1, q_temp2, W_n, g_2n_bar, g_n_bar, GMM_est)

		J_result[rr] = g_n_bar' * W_n * g_n_bar
		
	}
}

else {
// (Case 4): Dynamic model, not kink design (Default)

	for (rr=1; rr<=grid_num; rr=rr+1) {
	//This for-loop calculate 1st-step GMM w.r.t each point of the grid
	
		r_th = grid_th[rr]
		
		g_2n_bar = .
		g_n_bar = .
		GMM_est = .
		
		GMM_cal(g_1n_bar, r_th, Z_mat, x_temp1, x_temp2, x_temp3, ///
				q_temp1, q_temp2, W_n, g_2n_bar, g_n_bar, GMM_est)

		J_result[rr] = g_n_bar' * W_n * g_n_bar
		
	}
}

nouse_1 = . 
ind_min = .
minindex(J_result, 1, ind_min, nouse_1)
ind_min = ind_min[1] // If there are argmins more than two, we choose the first one.

r_hat = grid_th[ind_min]; // 1st-step threshold estimate


//Now with the estimated 1st-step threshold, we calculate 1st-step GMM estimator

if (flag_static == 1 & flag_kink == 0) {
	// (Case 1): Static model, not kink design
	GMM_cal(g_1n_bar, r_hat, Z_mat, x_temp1_c, x_temp2_c, x_temp3_c, ///
			   q_temp1, q_temp2, W_n, g_2n_bar, g_n_bar, GMM_est)
}
else if (flag_static == 0 & flag_kink == 1) {
	// (Case 2): Dynamic model, with kink
	GMM_cal_kink(g_1n_bar, r_hat, Z_mat, x_temp1, x_temp2, x_temp3, ///
			   q_temp1, q_temp2, W_n, g_2n_bar, g_n_bar, GMM_est)
}
else if (flag_static == 1 & flag_kink == 1) {
	// (Case 3): Static model, with kink 
	GMM_cal_kink(g_1n_bar, r_hat, Z_mat, x_temp1_c, x_temp2_c, x_temp3_c, ///
					q_temp1, q_temp2, W_n, g_2n_bar, g_n_bar, GMM_est)
}
else {
	// (Case 4): Dynamic model, not kink design (Default)
	GMM_cal(g_1n_bar, r_hat, Z_mat, x_temp1, x_temp2, x_temp3, ///
			q_temp1, q_temp2, W_n, g_2n_bar, g_n_bar, GMM_est)
}


//This is the 1st-step GMM estimate(GMM_est)


// ***(Step 3): 2nd-step GMM estimate*** //

// We collect estimated residuals for 4 difference cases
if (flag_static == 1 & flag_kink == 0) {
	ep_collect = ep_cal(r_hat, x_temp1_c, x_temp2_c, x_temp3_c, q_temp1, q_temp2, ///
						y_mat_rep, GMM_est)
	// Static, not kink - (Case 1)
}
else if (flag_static == 0 & flag_kink == 1) {
	ep_collect = ep_cal_kink(r_hat, x_temp1, x_temp2, x_temp3, q_temp1, q_temp2, ///
							 y_mat_rep, GMM_est)
	// Dynamic, kink - (Case 2)
}
else if (flag_static == 1 & flag_kink == 1) {
	ep_collect = ep_cal_kink(r_hat, x_temp1_c, x_temp2_c, x_temp3_c, q_temp1, q_temp2, ///
							 y_mat_rep, GMM_est)
	// Static, kink - (Case 3)
}
else {
	ep_collect = ep_cal(r_hat, x_temp1, x_temp2, x_temp3, q_temp1, q_temp2, ///
						y_mat_rep, GMM_est)
	// Dynamic, not kink (Case 4) (default)
}



W_n_2 = GMM_W_n_2_con(ep_collect, Z_mat)
J_result_2 = J(grid_num, 1, 0)

// Now we calculate 2nd step GMM estimator w.r.t. 4 different cases.

if (flag_static == 1 & flag_kink == 0) {
// (Case 1): Static model, not kink design

	for (rr2=1; rr2<=grid_num; rr2=rr2+1) {

		r_th = grid_th[rr2]
		
		g_2n_bar = .
		g_n_bar = .
		GMM_est_2 = .

		GMM_cal(g_1n_bar, r_th, Z_mat, x_temp1_c, x_temp2_c, x_temp3_c, ///
				   q_temp1, q_temp2, W_n_2, g_2n_bar, g_n_bar, GMM_est_2)

		J_result_2[rr2] = g_n_bar' * W_n_2 * g_n_bar
		
	}
}

else if (flag_static == 0 & flag_kink == 1) {
// (Case 2): Dynamic model, with kink

	for (rr2=1; rr2<=grid_num; rr2=rr2+1) {
	
		r_th = grid_th[rr2]
		
		g_2n_bar = .
		g_n_bar = .
		GMM_est_2 = .
		
		GMM_cal_kink(g_1n_bar, r_th, Z_mat, x_temp1, x_temp2, x_temp3, ///
					q_temp1, q_temp2, W_n_2, g_2n_bar, g_n_bar, GMM_est_2)

		J_result_2[rr2] = g_n_bar' * W_n_2 * g_n_bar
		
	}
}

else if (flag_static == 1 & flag_kink == 1) {
// (Case 3): Static model, with kink 

	for (rr2=1; rr2<=grid_num; rr2=rr2+1) {

		r_th = grid_th[rr2]
		
		g_2n_bar = .
		g_n_bar = .
		GMM_est_2 = .
		GMM_cal_kink(g_1n_bar, r_th, Z_mat, x_temp1_c, x_temp2_c, x_temp3_c, ///
						q_temp1, q_temp2, W_n_2, g_2n_bar, g_n_bar, GMM_est_2)

		J_result_2[rr2] = g_n_bar' * W_n_2 * g_n_bar
		
	}
}

else {
// (Case 4): Dynamic model, not kink design (Default)

	for (rr2=1; rr2<=grid_num; rr2=rr2+1) {
	//This for-loop calculate 1st-step GMM w.r.t each point of the grid
	
		r_th = grid_th[rr2]
		
		g_2n_bar = .
		g_n_bar = .
		GMM_est_2 = .
		
		GMM_cal(g_1n_bar, r_th, Z_mat, x_temp1, x_temp2, x_temp3, ///
				q_temp1, q_temp2, W_n_2, g_2n_bar, g_n_bar, GMM_est_2)

		J_result_2[rr2] = g_n_bar' * W_n_2 * g_n_bar
		
	}
}


nouse_2 = . 
ind_min_2 = .
minindex(J_result_2, 1, ind_min_2, nouse_2)

ind_min_2 = ind_min_2[1] // If there's more than two, we choose the first one.
r_hat_2 = grid_th[ind_min_2]; // 2nd-step threshold estimate

//Now with the estimated 2nd-step threshold, we calculate 2nd-step GMM estimator

if (flag_static == 1 & flag_kink == 0) {
	// (Case 1): Static model, not kink design
	GMM_cal(g_1n_bar, r_hat_2, Z_mat, x_temp1_c, x_temp2_c, x_temp3_c, ///
			   q_temp1, q_temp2, W_n_2, g_2n_bar, g_n_bar, GMM_est_2)
}
else if (flag_static == 0 & flag_kink == 1) {
	// (Case 2): Dynamic model, with kink 
	GMM_cal_kink(g_1n_bar, r_hat_2, Z_mat, x_temp1, x_temp2, x_temp3,  ///
				q_temp1, q_temp2, W_n_2, g_2n_bar, g_n_bar, GMM_est_2)
}
else if (flag_static == 1 & flag_kink == 1) {
	// (Case 3): Static model, with kink 
	GMM_cal_kink(g_1n_bar, r_hat_2, Z_mat, x_temp1_c, x_temp2_c, x_temp3_c, ///
					q_temp1, q_temp2, W_n_2, g_2n_bar, g_n_bar, GMM_est_2)
}
else {
	// (Case 4): Dynamic model, not kink design (Default)
	GMM_cal(g_1n_bar, r_hat_2, Z_mat, x_temp1, x_temp2, x_temp3, ///
			q_temp1, q_temp2, W_n_2, g_2n_bar, g_n_bar, GMM_est_2)
}

//This is the 2nd-step GMM estimate

// We collect 2nd-step estimated residuals for 4 difference cases
if (flag_static == 1 & flag_kink == 0) {
	ep_collect_2 = ep_cal(r_hat_2, x_temp1_c, x_temp2_c, x_temp3_c, q_temp1, q_temp2, ///
						  y_mat_rep, GMM_est_2)
	// Static, not kink - (Case 1)
}
else if (flag_static == 0 & flag_kink == 1) {
	ep_collect_2 = ep_cal_kink(r_hat_2, x_temp1, x_temp2, x_temp3, q_temp1, q_temp2, ///
							   y_mat_rep, GMM_est_2)
	// Dynamic, kink - (Case 2)
}
else if (flag_static == 1 & flag_kink == 1) {
	ep_collect_2 = ep_cal_kink(r_hat_2, x_temp1_c, x_temp2_c, x_temp3_c, q_temp1, q_temp2, ///
							   y_mat_rep, GMM_est_2)
	// Static, kink - (Case 3)
}
else {
	ep_collect_2 = ep_cal(r_hat_2, x_temp1, x_temp2, x_temp3, q_temp1, q_temp2, ///
						  y_mat_rep, GMM_est_2)
	// Dynamic, not kink (Case 4) (default)
}


omega_hat = GMM_W_n_2_con(ep_collect_2, Z_mat)
//This is actually an estimate of omega^(-1). Note the diffrence(with MATLAB ver.)


//Now we estimate covariance matrix for 4 cases

if (flag_static == 1 & flag_kink == 0) {
	// (Case 1): Static model, not kink design
	cov_mat = .
	cov_mat_con(GMM_est_2, r_hat_2, q_mat, Z_mat, x_temp1_c, x_temp2_c, ///
				x_temp3_c, q_temp1, q_temp2, h_0, omega_hat, cov_mat)
}
else if (flag_static == 0 & flag_kink == 1) {
	// (Case 2): Dynamic model, with kink 
	cov_mat = .
	cov_mat_con_kink(GMM_est_2, r_hat_2, q_mat, Z_mat, x_temp1, x_temp2, x_temp3, ///
					 q_temp1, q_temp2, h_0, omega_hat, cov_mat)
}
else if (flag_static == 1 & flag_kink == 1) {
	// (Case 3): Static model, with kink 
	cov_mat = .
	cov_mat_con_kink(GMM_est_2, r_hat_2, q_mat, Z_mat, x_temp1_c, x_temp2_c, ///
					 x_temp3_c, q_temp1, q_temp2, h_0, omega_hat, cov_mat)
}
else {
	// (Case 4): Dynamic model, not kink design (Default)
	cov_mat = .
	cov_mat_con(GMM_est_2, r_hat_2, q_mat, Z_mat, x_temp1, x_temp2, x_temp3, ///
				q_temp1, q_temp2, h_0, omega_hat, cov_mat)
}

GMM_estimates = (GMM_est_2 \ r_hat_2)


// ***(Step 4): Boostrap linearity test (optional)*** //

if (flag_b == 1) {
// (Case 1): Bootstrap needed
// First, with original data.
Z_mat_b = Z_mat
y_mat_rep_b = y_mat_rep
x_temp1_b = x_temp1
x_temp2_b = x_temp2
x_temp3_b = x_temp3
q_temp1_b = q_temp1
q_temp2_b = q_temp2

L = rows(Z_mat)

sup_wald_ori = . // This is the sup-Wald stat. for original sample
cov_wald_rec = .
sup_wald_cal(y_mat_rep_b, Z_mat_b, x_temp1_b, x_temp2_b, x_temp3_b, q_temp1_b, ///
			 q_temp2_b, moment_mat, grid_th, grid_num, flag_static, flag_kink, ///
			 N, T, sup_wald_ori, cov_wald_rec)

wald_rec = J(B, 1, 0) // To save sup-Wald*s


// Bootstrap sup-Wald stat.
for (bb2=1; bb2<=B; bb2=bb2+1) {

	sup_wald_b = .
	sup_wald_fb(y_mat_rep_b, Z_mat_b, x_temp1_b, x_temp2_b, x_temp3_b, q_temp1_b, ///
			 q_temp2_b, ep_collect_2, W_n_2, cov_wald_rec, moment_mat, grid_th, grid_num, ///
			 flag_static, flag_kink, N, T, sup_wald_b)

	wald_rec[bb2] = sup_wald_b
}

boots_p = sum(wald_rec :> sup_wald_ori)/B

}
else {
	// (Case 2): Bootstrap not needed
	boots_p = -1
}


// ** Display something, and calculate outputs (04/10), (06/25) ** //
if (flag_static == 1) {
	L = sum(moment_mat[., .]) - sum(moment_mat[., 2])
	}
else {
	L = sum(moment_mat[., .]) - sum(moment_mat[1, .])
	}
display(" ")
display("N = " + strofreal(N) + ", T = " + strofreal(T))
display("Panel Var. = " + pan_var)
display("Time Var. = " + time_var)
display("Number of moment conditions = " + strofreal(L))

if (flag_b == 1) {
	display("Bootstrap p-value for linearity test = " + strofreal(boots_p))
	}

z_975 = 1.959963984540054 // 0.975 quantile of Z
sds = (diagonal(cov_mat)):^(1/2) // S.E. of coefficients
CI_95 = (GMM_estimates - z_975*sds, GMM_estimates + z_975*sds)
// 95% asymptotic CI

}


/// **************************************************** ///
/// **************************************************** ///
/// ** The below lines are pre-defined mata functions ** ///
/// **************************************************** ///
/// **************************************************** ///


void grid_con(												///
				real matrix q_mat, real scalar trim_rate, 	///
				real scalar grid_num,						///
				grid_th)
{
	//This is a function to obtain the grid to use
	
	q_vec = vec(q_mat)
	grid_lower = mm_quantile(q_vec, 1, trim_rate/2)
	grid_upper = mm_quantile(q_vec, 1, 1-trim_rate/2)
	
	grid_th = rangen(grid_lower, grid_upper, grid_num) 
}


void FD_con(												///
				real matrix y_mat, real matrix x_mat, 		///
				real scalar T,								///
				y_mat_fd, x_mat_fd)
{
	//This is a function to obtain first-differenced variables
	y_mat_fd = y_mat[., 2..T] - y_mat[., 1..T-1]
	x_mat_fd = x_mat[., 2..T] - x_mat[., 1..T-1]
}


void Z_mat_con(												  ///
				real matrix y_mat, real matrix x_mat,		  ///
				real matrix x_mat_fd, real matrix moment_mat, ///
				real scalar k_endo, real scalar flag_static,  ///
				Z_mat)
{
	// Here we stack IV's for each t. Lagged y_{i, t}'s, differenced x_{i, t}, 
	// and 1's(constant) are included

	N = rows(y_mat)
	T = cols(y_mat)

	flag_z = 1
	//index indicator to fill out Z_mat

	k_1 = rows(x_mat) / N
	// # of dependent variables(excluding constant term)

	if (flag_static == 1) {
		//For static model 
		
		L = sum(moment_mat[., .]) - sum(moment_mat[., 2])
		//Here moment condition starts from t=2
		//However, the 2nd column is excluded since it is for lagged y_{i, t}
		//(The 1st row of moment_mat starts from t=2)
		
		Z_mat = J(L, N, 0) 
		
		for (kk1=1; kk1<=T-1; kk1=kk1+1) {
			m_temp = moment_mat[kk1, .]
			//Moment conditions for t=(kk1+1)
			
			for (ll1=1; ll1<=(k_1+1); ll1=ll1+1) {
				
				if(m_temp[ll1] == 0) {
					continue // No IV to stack.
				}
				else if(ll1 == 1){
					Z_mat[flag_z, .] = J(1, N, 1)
					
					flag_z = flag_z + 1
					// Stack 1's(constant)
				}
				else if(ll1 == 2){
					continue
					// Do not stack lagged y_{i, t}'s
				}
				else if (ll1 <= (k_1+1-k_endo)){
					Z_mat[flag_z, .] = (x_mat_fd[(ll1-2)*N+1..(ll1-1)*N, kk1])'
					
					flag_z = flag_z + 1
					// Stack \delta x_{i, t}'s for exogeneous x_{i, t}
				}
				else {
					Z_mat[flag_z..flag_z+m_temp[ll1]-1, .] = ///
					(x_mat[(ll1-2)*N+1..(ll1-1)*N, 1..m_temp[ll1]])'
					
					flag_z = flag_z + m_temp[ll1]
					// Stack x_{i, t}'s for endogenous x_{i, t}
				}
			}
		}
		
	}

	else {
		//For dynamic model (Default)
		
		L = sum(moment_mat[2..rows(moment_mat), .])
		//Here moment condition stars from t=3
		//(The 1st row of moment_mat starts from t=2)
		
		Z_mat = J(L, N, 0) 
		
		for (kk1=2; kk1<=T-1; kk1=kk1+1) {
			m_temp = moment_mat[kk1, .]
			//Moment conditions for t=(kk1+1)
			
			for (ll1=1; ll1<=(k_1+1); ll1=ll1+1) {
				
				if(m_temp[ll1] == 0) {
					continue // No IV to stack. But this will not happen normally.
				}
				else if(ll1 == 1){
					Z_mat[flag_z, .] = J(1, N, 1)
					
					flag_z = flag_z + 1
					// Stack 1's(constant)
				}
				else if(ll1 == 2){
					Z_mat[flag_z..flag_z+m_temp[ll1]-1, .] = ///
					(y_mat[., kk1-m_temp[ll1]..kk1-1])'
					
					flag_z = flag_z + m_temp[ll1]
					// Stack lagged y_{i, t}'s
				}
				else if (ll1 <= (k_1+1-k_endo)){
					Z_mat[flag_z, .] = (x_mat_fd[(ll1-2)*N+1..(ll1-1)*N, kk1])'
					
					flag_z = flag_z + 1
					// Stack \delta x_{i, t}'s for exogeneous x_{i, t}
				}
				else {
					Z_mat[flag_z..flag_z+m_temp[ll1]-1, .] = ///
					(x_mat[(ll1-2)*N+1..(ll1-1)*N, 1..m_temp[ll1]])'
					
					flag_z = flag_z + m_temp[ll1]
					// Stack x_{i, t}'s for endogenous x_{i, t}
				}
			}
		}
		
	}

}



void GMM_var_con(											///
				real matrix y_mat_fd, real matrix x_mat_fd, ///
				real matrix x_mat, real matrix q_mat,		///
				real matrix moment_mat, real scalar T,		///
				real scalar flag_static, 					///
				y_mat_rep, x_temp1, x_temp2, x_temp3, q_temp1, q_temp2)
{
	//This is a function to stack data for GMM calculation
	
	if (flag_static == 1) {
		L = sum(moment_mat) - sum(moment_mat[., 2])
		// In case of static model
	}
	else {
		L = sum(moment_mat[2..rows(moment_mat), .])
		//In case of dynamic model(default)
		//The first row was excluded since it is for t=2 
		//*For t=2 (LHS) is (y_{i, 2} - y_{i, 1}), so there's no proper (RHS)
	}



	N = rows(y_mat_fd)
	k_1 = rows(x_mat)/N

	y_mat_rep = J(L, N, 0)
	x_temp1 = J(L, N*k_1, 0)
	x_temp2 = J(L, N*k_1, 0)
	x_temp3 = J(L, N*k_1, 0)
	q_temp1 = J(L, N, 0)
	q_temp2 = J(L, N, 0)

	flag_r1 = 1

	if (flag_static == 1) {
		//For static model
		
		for (kk=1; kk<=T-1; kk=kk+1) {
			m_temp = sum(moment_mat[kk, .]) - sum(moment_mat[kk, 2])
			// # of stacked moment condition for t = kk

			if (m_temp == 0) {
				continue
			}

			else {

			y_mat_rep[flag_r1..flag_r1+m_temp-1, .] = (y_mat_fd[., kk] # J(1, m_temp, 1))'
			//Using # (kronecker product), repeat the desired column
		 
			x_temp1[flag_r1..flag_r1+m_temp-1, .] = (x_mat_fd[., kk] # J(1, m_temp, 1))'
			x_temp2[flag_r1..flag_r1+m_temp-1, .] = (x_mat[., kk+1] # J(1, m_temp, 1))'
			x_temp3[flag_r1..flag_r1+m_temp-1, .] = (x_mat[., kk] # J(1, m_temp, 1))'
			
			q_temp1[flag_r1..flag_r1+m_temp-1, .] = (q_mat[., kk+1] # J(1, m_temp, 1))'
			q_temp2[flag_r1..flag_r1+m_temp-1, .] = (q_mat[., kk] # J(1, m_temp, 1))'
		 
			flag_r1 = flag_r1 + m_temp
			} // end of if-else
			
		} // end of for-loop
		

	}
	else {
	//For dynamic model (Default)

		for (kk=2; kk<=T-1; kk=kk+1) {
			m_temp = sum(moment_mat[kk, .])  
			// # of stacked moment condition for t = kk

			if (m_temp == 0) {
				continue
			}

			else {

			y_mat_rep[flag_r1..flag_r1+m_temp-1, .] = (y_mat_fd[., kk] # J(1, m_temp, 1))'
			//Using # (kronecker product), repeat the desired column
		 
			x_temp1[flag_r1..flag_r1+m_temp-1, .] = (x_mat_fd[., kk] # J(1, m_temp, 1))'
			x_temp2[flag_r1..flag_r1+m_temp-1, .] = (x_mat[., kk+1] # J(1, m_temp, 1))'
			x_temp3[flag_r1..flag_r1+m_temp-1, .] = (x_mat[., kk] # J(1, m_temp, 1))'
			
			q_temp1[flag_r1..flag_r1+m_temp-1, .] = (q_mat[., kk+1] # J(1, m_temp, 1))'
			q_temp2[flag_r1..flag_r1+m_temp-1, .] = (q_mat[., kk] # J(1, m_temp, 1))'
		 
			flag_r1 = flag_r1 + m_temp
			} // end of if-else
			
		} // end of for-loop
	}

}


void GMM_W_n_con(											///
				real matrix Z_mat, real matrix moment_mat,	///
				real scalar T,								///
				W_n)
{
	//This is a function to construct weight matrix for 1st step
	//This is for dynamic model (default)
	
	L = rows(Z_mat)
	N = cols(Z_mat)
	
	W_n_temp1 = J(L, L, 0) 
	W_n_temp2 = J(L, L, 0)
	
	flag_r2 = 1 //To indicate the next column to fill up
	flag_z1 = 1
	flag_z2 = 1 //To indicate the index of Z_mat to copy
	
	moment_vec = rowsum(moment_mat)
	//for dynamic model(default)
		
	for (kk2=3; kk2<=T; kk2=kk2+1) {
			m_num = moment_vec[kk2-1] // # of IV's at period t (where LHS is y_{i, t} - y_{i, t-1})
			
			if (m_num == 0) {
				continue 
				// When there's no IV. But this will not happen normally.
				}
			
			else{
				if(flag_r2 == 1) {
					Z_temp1 = Z_mat[1..m_num, .]
					// (m_1 * N) matrix, where m_1 is # of moment condition at t_0
					
					W_n_temp2[1..m_num, 1..m_num] = (1/N)*(Z_temp1*(Z_temp1'))
					flag_r2 = flag_r2 + m_num
					flag_z2 = flag_z2 + m_num
					}
				
				else {
					m_num_l = moment_vec[kk2-2] // lagged m_num
					Z_temp1 = Z_mat[flag_z1..flag_z1+m_num_l-1, .] 
					Z_temp2 = Z_mat[flag_z2..flag_z2+m_num-1, .]
					//Both are building blocks of (10) in Seo & Shin(2016).
					
					flag_z1 = flag_z1 + m_num_l
					flag_z2 = flag_z2 + m_num	// updating
					
					W_n_temp1[flag_r2-(m_num_l)..flag_r2-1, flag_r2..flag_r2+m_num-1] ///
						= -(1/N)*(Z_temp1*Z_temp2')
					
					 W_n_temp2[flag_r2..flag_r2+m_num-1, flag_r2..flag_r2+m_num-1] ///
						= (1/N)*(Z_temp2 * Z_temp2')
					
					flag_r2 = flag_r2 + m_num
					}
				}
			
		}


	W_n = invsym(W_n_temp1 + W_n_temp1' + 2*W_n_temp2)
}


void GMM_cal(											///
			real matrix g_1n_bar, real scalar r_th,		///
			real matrix Z_mat, real matrix x_temp1,		///
			real matrix x_temp2, real matrix x_temp3,	///
			real matrix q_temp1, real matrix q_temp2,	///
			real matrix W_n,							///
			g_2n_bar, g_n_bar, GMM_est)
{
	//This is a function to calculate GMM estimator for given threshold r_th 
	
	L = rows(q_temp1)
	N = cols(q_temp1)
	k_1 = cols(x_temp2)/N
	
	ind_temp2 = (q_temp1 :> r_th)
	ind_temp3 = (q_temp2 :> r_th) 
	// This two (L*N) 0-1 matrices contain information whether q_{it} or q_{it-1} is 
	// larger than the threshold r_th
	
	g_2n_bar = J(L, 2*k_1 + 1, 0) //This is identical to g_2n_bar in Seo & Shin(2016)
	g_2n_bar[., k_1 + 1] = rowsum(Z_mat:*(ind_temp2 - ind_temp3))/N
	
	for (kk3=1; kk3<=k_1; kk3=kk3+1) {
		g_2n_bar[., kk3] = rowsum(Z_mat:*x_temp1[., (kk3-1)*N+1..kk3*N])/N
		g_2n_bar[., k_1 + 1 + kk3] = rowsum(Z_mat:*(x_temp2[., (kk3-1)*N+1..kk3*N]:*ind_temp2 ///
											- x_temp3[., (kk3-1)*N+1..kk3*N]:*ind_temp3))/N
	}
	
	GMM_est = invsym(g_2n_bar' * W_n * g_2n_bar) * (g_2n_bar' * W_n * g_1n_bar)
	g_n_bar = g_1n_bar - g_2n_bar*GMM_est;
}


function ep_cal(
				real scalar r_hat, real matrix x_temp1,		///
				real matrix x_temp2, real matrix x_temp3,	///
				real matrix q_temp1, real matrix q_temp2,	///
				real matrix y_mat_rep, real matrix GMM_est)
{
	// This is a function to stack estimated errors from the 1st step

	L = rows(q_temp1)
	N = cols(q_temp1)
	k_1 = cols(x_temp2)/N
	
	ep_collect = J(L, N, 0)
	//(L*N) matrix that will contain all collected residuals from the first GMM
	
	ind_temp2 = (q_temp1 :> r_hat)
	ind_temp3 = (q_temp2 :> r_hat) 
	
	ind_vec = rangen(0, (k_1-1)*N, k_1)
	// To pick general k_1 variables
	
	
	for (ii2=1; ii2<=L; ii2=ii2+1) {

		for (ii3=1; ii3<=N; ii3=ii3+1) {
					
					ii3_rep = J(k_1, 1, ii3)
					
					ep_collect[ii2, ii3] = y_mat_rep[ii2, ii3] - (GMM_est')* /// 
					(x_temp1[ii2, ind_vec+ii3_rep]' \ 						 ///
					ind_temp2[ii2, ii3] - ind_temp3[ii2, ii3] \ 			 ///
					(ind_temp2[ii2, ii3]:*x_temp2[ii2, ind_vec+ii3_rep] - 	 ///
					ind_temp3[ii2, ii3]:*x_temp3[ii2, ind_vec+ii3_rep])')
		
		}

	}
	
	return (ep_collect)
}


function ep_cal_kink(
					real scalar r_hat, real matrix x_temp1,		///
					real matrix x_temp2, real matrix x_temp3,	///
					real matrix q_temp1, real matrix q_temp2,	///
					real matrix y_mat_rep, real matrix GMM_est)
{
	// This is a function to stack estimated errors from the 1st step
	// Especially, for the kink model case - i.e. only the variable 
	// (q_{i, t} - r) * 1_{q_{i, t} > r} changes according to threshold

	L = rows(q_temp1)
	N = cols(q_temp1)
	k_1 = cols(x_temp2)/N
	
	ep_collect = J(L, N, 0)
	//(L*N) matrix that will contain all collected residuals from the first GMM
	
	ind_temp2 = (q_temp1 :> r_hat)
	ind_temp3 = (q_temp2 :> r_hat) 
	
	ind_vec = rangen(0, (k_1-1)*N, k_1)
	// To pick general k_1 variables
	
	
	for (ii2=1; ii2<=L; ii2=ii2+1) {

		for (ii3=1; ii3<=N; ii3=ii3+1) {
					
					ii3_rep = J(k_1, 1, ii3)
					
					ep_collect[ii2, ii3] = y_mat_rep[ii2, ii3] - (GMM_est')* /// 
					(x_temp1[ii2, ind_vec+ii3_rep]' \ 						 ///
					(ind_temp2[ii2, ii3]*(q_temp1[ii2, ii3] - r_hat) - 	 	 ///
					ind_temp3[ii2, ii3]*(q_temp2[ii2, ii3] - r_hat))')
		
		}

	}
	
	return (ep_collect)
}


function GMM_W_n_2_con(										///
						real matrix ep_collect, real matrix Z_mat)
{
	//This is a function to construct weight matrix for 2nd step
	
	L = rows(Z_mat)
	N = cols(Z_mat)
	
	g_temp_1 = J(L, L, 0)
	g_temp_2 = J(L, 1, 0) 
	//Matrices that are corresponding to the 1st, 2nd term of (11) in the paper. 

	for (ii3=1; ii3<=N; ii3=ii3+1) {
		g_i_hat = (ep_collect[., ii3]):*(Z_mat[., ii3]);

		g_temp_1 = g_temp_1 + g_i_hat*g_i_hat';
		g_temp_2 = g_temp_2 + g_i_hat;
	
	}

	W_n_2 = invsym((1/N)*g_temp_1 - (1/(N^2))*(g_temp_2*g_temp_2'))
	return(W_n_2);
}

void moment_mat_con(real scalar T, real scalar k_1,		///
					real scalar k_endo, moment_mat)
{

	if (k_endo == 0) {
		moment_mat = J(T-1, k_1 + 1, 1)
		moment_mat[., 2] = (0::T-2)
	/* 
	This matrix contains information about moment conditions. 
	(Here we assume there's no endogenous x_{i, t})
	The first row is for \delta x_{i2} = x_{i2} - x_{i1}
	The second row is for \delta x_{i3} = x_{i3} - x_{i2}
	... and the last row is for \delta x_{iT} = x_{iT} - x_{iT-1}

	The first column is for use of constant(1's).
	The second column is for use of lagged y_{i, t}'s. 
	(For different t's we use different lagged terms of y_{i, t})
	The other columns are for use of differenced x_{i, t}'s
	*/
	}
	else {
		moment_mat = J(T-1, k_1 + 1, 1)
		moment_mat[., 2] = (0::T-2)
		moment_mat[., (k_1 + 2 - k_endo) .. (k_1 + 1)] = J(1, k_endo, (0::T-2))
	/*
	Additionally, when there are (# = k_endo) endogenous x_{i, t},
	The last columns are for use of endogenous lagged x_{i, t}'s
	*/
	}

}

void GMM_W_n_static(											///
					real matrix Z_mat, real matrix moment_mat,	///
					real scalar T,								///
					W_n)
{
	//This is a function to construct weight matrix for 1st step
	//This is for static model (optional)
	
	L = rows(Z_mat)
	N = cols(Z_mat)
	
	W_n_temp1 = J(L, L, 0) 
	W_n_temp2 = J(L, L, 0)
	
	flag_r2 = 1 //To indicate the next column to fill up
	flag_z1 = 1
	flag_z2 = 1 //To indicate the index of Z_mat to copy
	
	moment_vec = rowsum(moment_mat) - moment_mat[., 2]
	//for static model
	
	for (kk2=2; kk2<=T; kk2=kk2+1) {
			m_num = moment_vec[kk2-1] 
			// # of IV's at period t (where LHS is y_{i, t} - y_{i, t-1})
			
			if (m_num == 0) {
				continue 
				// When there's no IV. But this will not happen normally.
				}
			
			else{
				if(flag_r2 == 1) {
					Z_temp1 = Z_mat[1..m_num, .]
					// (m_1 * N) matrix, where m_1 is # of moment condition at t_0
					
					W_n_temp2[1..m_num, 1..m_num] = (1/N)*(Z_temp1*(Z_temp1'))
					flag_r2 = flag_r2 + m_num
					flag_z2 = flag_z2 + m_num
					}
				
				else {
					m_num_l = moment_vec[kk2-2] // lagged m_num
					Z_temp1 = Z_mat[flag_z1..flag_z1+m_num_l-1, .] 
					Z_temp2 = Z_mat[flag_z2..flag_z2+m_num-1, .]
					//Both are building blocks of (10) in Seo & Shin(2016).
					
					flag_z1 = flag_z1 + m_num_l
					flag_z2 = flag_z2 + m_num	// updating
					
					W_n_temp1[flag_r2-(m_num_l)..flag_r2-1, flag_r2..flag_r2+m_num-1] ///
						= -(1/N)*(Z_temp1*Z_temp2')
					
					 W_n_temp2[flag_r2..flag_r2+m_num-1, flag_r2..flag_r2+m_num-1] ///
						= (1/N)*(Z_temp2 * Z_temp2')
					
					flag_r2 = flag_r2 + m_num
					}
				}
			
		}


	W_n = invsym(W_n_temp1 + W_n_temp1' + 2*W_n_temp2)
}


void GMM_cal_kink(													///
					real matrix g_1n_bar, real scalar r_th,			///
					real matrix Z_mat, real matrix x_temp1_c,		///
					real matrix x_temp2_c, real matrix x_temp3_c,	///
					real matrix q_temp1, real matrix q_temp2,		///
					real matrix W_n,								///
					g_2n_bar, g_n_bar, GMM_est)
{
	//This is a function to calculate GMM estimator for given threshold r_th 
	//Especially, this is for static model with kink - (Case 2, 3)
	
	L = rows(q_temp1)
	N = cols(q_temp1)
	k_1 = cols(x_temp1_c)/N
	
	ind_temp2 = (q_temp1 :> r_th)
	ind_temp3 = (q_temp2 :> r_th) 
	// This two (L*N) 0-1 matrices contain information whether q_{it} or q_{it-1} is 
	// larger than the threshold r_th
	
	g_2n_bar = J(L, k_1 + 1, 0) //This is identical to g_2n_bar in Seo & Shin(2016)
	
	
	for (kk3=1; kk3<=k_1; kk3=kk3+1) {
		g_2n_bar[., kk3] = rowsum(Z_mat:*x_temp1_c[., (kk3-1)*N+1..kk3*N])/N
	}
	
	kink_var = ((q_temp1 :- r_th):*ind_temp2 - (q_temp2 :- r_th):*ind_temp3)
	// This mat. contains the variable for kink model, which can be described as
	// (q_{i, t} - r) * 1_{q_{i, t} > r} - (q_{i, t-1} - r) * 1_{q_{i, t-1} > r}
	// Note that the above term is after first-differencing
	
	g_2n_bar[., k_1 + 1] = rowsum(Z_mat:*kink_var)/N
	
	GMM_est = invsym(g_2n_bar' * W_n * g_2n_bar) * (g_2n_bar' * W_n * g_1n_bar)
	g_n_bar = g_1n_bar - g_2n_bar*GMM_est;
}


void Z_mat_IV_add(												
				  real matrix Z_mat, real matrix moment_mat, ///
				  real matrix IV_mat, real scalar k_inst,	 ///
				  real scalar flag_static)
{
	//This is a function to extend Z_mat and moment_mat, 
	//so that they can respond to additional given IV's. 
	//Activated only if there are additional IV's (optional)
	
	L = rows(Z_mat)
	
	N = cols(Z_mat)
	T = rows(moment_mat) + 1
	
	moment_add = J(T-1, k_inst, 1)
	//Indicates IV's to add
	
	flag_z1 = 1 
	flag_z2 = 1
	//Index for Z_mat and new Z_mat, respectively
	
	Z_copy = Z_mat // copy of original Z_mat
	
	if (flag_static == 1){
	//For static model(special case)
	
	Z_mat = J(L + sum(moment_add), N, 0) 
	// This will be new Z_mat for static model
	
		for (kk4=2; kk4<=T; kk4=kk4+1) {
			m_i1 = sum(moment_mat[kk4-1, .]) - moment_mat[kk4-1, 2] 
			// # of original IV's at t = kk4
			
			IV_add = colshape(IV_mat[., kk4], N)
			//(k_inst * N) matrix, this will be added
			
			IV_ori = Z_copy[flag_z1..flag_z1+m_i1-1, .]
			// original IV's at t = kk4
			
			flag_z1 = flag_z1 + m_i1 // Index updating
			
			Z_mat[flag_z2..flag_z2+m_i1-1, .] = IV_ori
			Z_mat[flag_z2+m_i1..flag_z2+m_i1+k_inst-1, .] = IV_add
			
			flag_z2 = flag_z2 + m_i1 + k_inst
		}
	}
	else {
	//For dynamic model
	
	Z_mat = J(L + sum(moment_add) - k_inst, N, 0) 
	// This will be new Z_mat for dynamic model
	
		for (kk4=3; kk4<=T; kk4=kk4+1) {
			m_i1 = sum(moment_mat[kk4-1, .])
			// # of original IV's at t = kk4
			
			IV_add = colshape(IV_mat[., kk4], N) 
			//(k_inst * N) matrix, this will be added
			
			IV_ori = Z_copy[flag_z1..flag_z1+m_i1-1, .]
			// original IV's at t = kk4
			
			flag_z1 = flag_z1 + m_i1 // Index updating
			
			Z_mat[flag_z2..flag_z2+m_i1-1, .] = IV_ori
			Z_mat[flag_z2+m_i1..flag_z2+m_i1+k_inst-1, .] = IV_add
			
			flag_z2 = flag_z2 + m_i1 + k_inst
		}
	}
	
	moment_mat = (moment_mat, moment_add)
	//Add columns for additional IV
	
}


void cov_mat_con(												///
					real matrix GMM_est_2, real scalar r_hat_2,	///
					real matrix q_mat,							///
					real matrix Z_mat, real matrix x_temp1, 	///
					real matrix x_temp2, real matrix x_temp3,	///
					real matrix q_temp1, real matrix q_temp2,	///
					real scalar h_0, real matrix omega_hat,		///
					cov_mat)
{
	//This function estimates covariance matrix. First, we construct ingredients.
	L = rows(Z_mat)
	N = cols(Z_mat)
	h_band = h_0 * 1.06 * N^(-0.2) * sqrt(variance(colshape(q_mat, 1)))
	k_1 = cols(x_temp1)/N
	
	G_b = J(L, k_1, 0)
	G_d = J(L, k_1 + 1, 0)
	G_r = J(L, 1, 0)
	//Three matrices are \hat{G_b}, \hat{G_d}, \hat{G_r} on Seo & Shin(2016), 173p. 
	
	ind_temp2 = (q_temp1 :> r_hat_2)
	ind_temp3 = (q_temp2 :> r_hat_2) 
	// This two (L*N) 0-1 matrices contain information whether q_{it} or q_{it-1} is 
	// larger than the threshold r_th
	
	g_2n_bar = J(L, 2*k_1 + 1, 0) //This is identical to g_2n_bar in Seo & Shin(2016)
	g_2n_bar[., k_1 + 1] = rowsum(Z_mat:*(ind_temp2 - ind_temp3))/N
	

	for (ii6=1; ii6<=k_1; ii6=ii6+1) {
		temp_x = x_temp1[., (ii6-1)*N+1..ii6*N] //Extract (L*N) for one variable
		G_b[., ii6] = -rowsum(Z_mat:*x_temp1[., (ii6-1)*N+1 .. ii6*N])/N
		
		g_2n_bar[., ii6] = rowsum(Z_mat:*temp_x)/N
		g_2n_bar[., k_1 + 1 + ii6] = rowsum(Z_mat:*(x_temp2[., (ii6-1)*N+1..ii6*N]:*ind_temp2 ///
											- x_temp3[., (ii6-1)*N+1..ii6*N]:*ind_temp3))/N
	}
	
	G_d[., 1] = -g_2n_bar[., k_1 + 1]
	G_d[., 2..k_1 + 1] = -g_2n_bar[., k_1 + 2..2*k_1 + 1]
	
	
	for (ii7=1; ii7<=N; ii7=ii7+1){
		ind_vec = N*(0::k_1-1) :+ ii7
		
		temp_rc1 = (J(L, 1, 1), x_temp3[., ind_vec])
		temp_rc2 = normalden((r_hat_2:-q_temp2[., ii7])/h_band) # J(1, k_1+1, 1)
		temp_rc3 = (J(L, 1, 1), x_temp2[., ind_vec])
		temp_rc4 = normalden((r_hat_2:-q_temp1[., ii7])/h_band) # J(1, k_1+1, 1)
		
		temp_z = (temp_rc1:*temp_rc2 - temp_rc3:*temp_rc4) * GMM_est_2[k_1+1..2*k_1+1]
		G_r = G_r + Z_mat[., ii7]:*temp_z
		
	}
	
	G_r = G_r/(N*h_band)
	
	G = (G_b, G_d, G_r)
	cov_mat = cholinv(G'*omega_hat*G)/N
}


void cov_mat_con_kink(												///
						real matrix GMM_est_2, real scalar r_hat_2,	///
						real matrix q_mat,							///
						real matrix Z_mat, real matrix x_temp1, 	///
						real matrix x_temp2, real matrix x_temp3,	///
						real matrix q_temp1, real matrix q_temp2,	///
						real scalar h_0, real matrix omega_hat,		///
						cov_mat)
{
	//This function estimates covariance matrix. First, we construct ingredients.
	L = rows(Z_mat)
	N = cols(Z_mat)
	h_band = h_0 * 1.06 * N^(-0.2) * sqrt(variance(colshape(q_mat, 1)))
	k_1 = cols(x_temp1)/N
	
	G_b = J(L, k_1, 0)
	G_d = J(L, 1, 0)
	G_r = J(L, 1, 0)
	//Three matrices are \hat{G_b}, \hat{G_d}, \hat{G_r} on Seo & Shin(2016), 173p. 
	
	ind_temp2 = (q_temp1 :> r_hat_2)
	ind_temp3 = (q_temp2 :> r_hat_2) 
	// This two (L*N) 0-1 matrices contain information whether q_{it} or q_{it-1} is 
	// larger than the threshold r_th
	
	G_d[., 1] = rowsum(Z_mat:*((q_temp1:-r_hat_2):*ind_temp2 - ///
							   (q_temp2:-r_hat_2):*ind_temp3))/N
	// Adjusted for kink model
	
	for (ii6=1; ii6<=k_1; ii6=ii6+1) {
		temp_x = x_temp1[., (ii6-1)*N+1..ii6*N] //Extract (L*N) for one variable
		G_b[., ii6] = -rowsum(Z_mat:*temp_x)/N
	}
	
	d_0 = GMM_est_2[k_1 + 1]
	G_r = -d_0:*rowsum(Z_mat:*(ind_temp2 - ind_temp3))/N
	
	G = (G_b, G_d, G_r)
	cov_mat = cholinv(G'*omega_hat*G)/N
}


void index_non_ms( 
					real scalar N, real scalar T, real matrix PT_mat, 	///
					real scalar N_ori, real matrix ind_non_ms)
{
	// The derived ind_non_ms is column vector with binary elements,
	// 1 for samples without missing through all times, 0 for others.
	// Using ind_non_ms, I will filter observations with missing values (for any t)
	
	
	label_vec = uniqrows(PT_mat[., 1]) // Unique labels
	time_vec = uniqrows(PT_mat[., 2])

	N_obs = rows(PT_mat) // # of total observations

	N = length(label_vec) 
	T = length(time_vec) // N, T of panel data
	
	ind_non_ms = J(rows(PT_mat), 1, 1)
	
	flag_ms = 1 // For indexing
	
	for (ii_ms=1; ii_ms<=N; ii_ms = ii_ms + 1) {
		ilab = label_vec[ii_ms] // The label of i-th sample
		
		sub_ilab = select(PT_mat[., 2], (PT_mat[., 1] :== J(N_obs, 1, ilab)))
		// # of (time-varying) obs. for ilabel
		
		obs_ilab = length(sub_ilab)
		
		if (obs_ilab < T) {
			ind_non_ms[flag_ms..flag_ms+obs_ilab-1] = J(obs_ilab, 1, 0)
			flag_ms = flag_ms + obs_ilab
		}
		else {
			flag_ms = flag_ms + T
		}
		
	}
	
	N_ori = N
	N = sum(ind_non_ms) / T

}




void sup_wald_cal(													 	 ///
						real matrix y_mat_rep_b, real matrix Z_mat_b,	 ///
						real matrix x_temp1_b, real matrix x_temp2_b, 	 ///
						real matrix x_temp3_b,						 ///
						real matrix q_temp1_b, real matrix q_temp2_b, ///
						real matrix moment_mat,	real matrix grid_th, ///
						real scalar grid_num,						 ///
						real scalar flag_static,					 /// 
						real scalar flag_kink,						 ///
						real scalar N, real scalar T,				 ///
						sup_wald, cov_wald_rec)
{	
	k_1 = cols(x_temp1_b)/N
		
	// We separate delta_hat from the GMM_estimate
	// Note that delta_hat includes constant difference
	if (flag_static == 1 & flag_kink == 0) {
		// (Case 1): Static model, not kink design
		dlt_len = k_1
	}
	else if (flag_static == 0 & flag_kink == 1) {
		// (Case 2): Dynamic model, with kink 
		dlt_len = 2
	}
	else if (flag_static == 1 & flag_kink == 1) {
		// (Case 3): Static model, with kink 
		dlt_len = 2
	}
	else {
		// (Case 4): Dynamic model, not kink design (Default)
		dlt_len = k_1 + 1
	}
	
	cov_wald_rec = J(grid_num*dlt_len, dlt_len, 0)
	// To record sigma_hat for each r
	
	// ***(Step 2): 1st-step GMM estimate*** //
	W_n_b = .

	if (flag_static == 1) {
		GMM_W_n_static(Z_mat_b, moment_mat, T, W_n_b)
	}
	else {
		GMM_W_n_con(Z_mat_b, moment_mat, T, W_n_b)
	}


	g_1n_bar_b = rowsum(Z_mat_b:*y_mat_rep_b)/N 
	//This is row-wise sum. Identical to g_1n_bar in Seo & Shin(2016)
	 
	// J_result = J(grid_num, 1, 0) // We don't need this for bootstrapping. r is given.

	if (flag_static == 1){
		end_p = cols(x_temp1_b);
		
		x_temp1_c_b = x_temp1_b[., N+1..end_p]
		x_temp2_c_b = x_temp2_b[., N+1..end_p]
		x_temp3_c_b = x_temp3_b[., N+1..end_p]
		// Copies of matrices, after removing columns that contain
		// y_{i, t-1} related terms
		
		// These matrices will be used if we assume static model
	}


	//Now with the given threshold, we calculate 1st-step GMM estimator
	//In the following for-loop we will calculate sup Wald.

	wald_vec = J(grid_num, 1, 0)

	for (bb=1; bb<=grid_num; bb=bb+1) {
	// Specially I omit indentation for this loop (too long!)

	r_tt = grid_th[bb]

	g_2n_bar_b = .
	g_n_bar_b = .
	GMM_est_b = .
			
	if (flag_static == 1 & flag_kink == 0) {
		// (Case 1): Static model, not kink design
		GMM_cal(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_c_b, x_temp2_c_b, x_temp3_c_b, ///
				   q_temp1_b, q_temp2_b, W_n_b, g_2n_bar_b, g_n_bar_b, GMM_est_b)
	}
	else if (flag_static == 0 & flag_kink == 1) {
		// (Case 2): Dynamic model, with kink
		GMM_cal_kink(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_b, x_temp2_b, x_temp3_b, ///
				   q_temp1_b, q_temp2_b, W_n_b, g_2n_bar_b, g_n_bar_b, GMM_est_b)
	}
	else if (flag_static == 1 & flag_kink == 1) {
		// (Case 3): Static model, with kink 
		GMM_cal_kink(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_c_b, x_temp2_c_b, x_temp3_c_b, ///
						q_temp1_b, q_temp2_b, W_n_b, g_2n_bar_b, g_n_bar_b, GMM_est_b)
	}
	else {
		// (Case 4): Dynamic model, not kink design (Default)
		GMM_cal(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_b, x_temp2_b, x_temp3_b, ///
				q_temp1_b, q_temp2_b, W_n_b, g_2n_bar_b, g_n_bar_b, GMM_est_b)
	}

	//This is the 1st-step bootstrap GMM estimate(GMM_est_b)
	


	// ***(Step 3): 2nd-step GMM estimate*** //

	// We collect estimated residuals for 4 difference cases
	if (flag_static == 1 & flag_kink == 0) {
		ep_collect_b = ep_cal(r_tt, x_temp1_c_b, x_temp2_c_b, x_temp3_c_b, q_temp1_b, q_temp2_b, ///
							y_mat_rep_b, GMM_est_b)
		// Static, not kink - (Case 1)
	}
	else if (flag_static == 0 & flag_kink == 1) {
		ep_collect_b = ep_cal_kink(r_tt, x_temp1_b, x_temp2_b, x_temp3_b, q_temp1_b, q_temp2_b, ///
								 y_mat_rep_b, GMM_est_b)
		// Dynamic, kink - (Case 2)
	}
	else if (flag_static == 1 & flag_kink == 1) {
		ep_collect_b = ep_cal_kink(r_tt, x_temp1_c_b, x_temp2_c_b, x_temp3_c_b, q_temp1_b, q_temp2_b, ///
								 y_mat_rep_b, GMM_est_b)
		// Static, kink - (Case 3)
	}
	else {
		ep_collect_b = ep_cal(r_tt, x_temp1_b, x_temp2_b, x_temp3_b, q_temp1_b, q_temp2_b, ///
							y_mat_rep_b, GMM_est_b)
		// Dynamic, not kink (Case 4) (default)
	}

	W_n_2_b = GMM_W_n_2_con(ep_collect_b, Z_mat_b)


	//Now with the given threshold r_tt, we calculate 2nd-step GMM estimator

	g_2n_bar_b = .
	g_n_bar_b = .
	GMM_est_2_b = .
			
	if (flag_static == 1 & flag_kink == 0) {
		// (Case 1): Static model, not kink design
		GMM_cal(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_c_b, x_temp2_c_b, x_temp3_c_b, ///
				   q_temp1_b, q_temp2_b, W_n_2_b, g_2n_bar_b, g_n_bar_b, GMM_est_2_b)
	}
	else if (flag_static == 0 & flag_kink == 1) {
		// (Case 2): Dynamic model, with kink 
		GMM_cal_kink(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_b, x_temp2_b, x_temp3_b,  ///
					q_temp1_b, q_temp2_b, W_n_2_b, g_2n_bar_b, g_n_bar_b, GMM_est_2_b)
	}
	else if (flag_static == 1 & flag_kink == 1) {
		// (Case 3): Static model, with kink 
		GMM_cal_kink(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_c_b, x_temp2_c_b, x_temp3_c_b, ///
						q_temp1_b, q_temp2_b, W_n_2_b, g_2n_bar_b, g_n_bar_b, GMM_est_2_b)
	}
	else {
		// (Case 4): Dynamic model, not kink design (Default)
		GMM_cal(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_b, x_temp2_b, x_temp3_b, ///
				q_temp1_b, q_temp2_b, W_n_2_b, g_2n_bar_b, g_n_bar_b, GMM_est_2_b)
	}

	//This is the 2nd-step GMM estimate

	// We collect 2nd-step estimated residuals for 4 difference cases
	if (flag_static == 1 & flag_kink == 0) {
		ep_collect_2_b = ep_cal(r_tt, x_temp1_c_b, x_temp2_c_b, x_temp3_c_b, q_temp1_b, q_temp2_b, ///
							  y_mat_rep_b, GMM_est_2_b)
		// Static, not kink - (Case 1)
	}
	else if (flag_static == 0 & flag_kink == 1) {
		ep_collect_2_b = ep_cal_kink(r_tt, x_temp1_b, x_temp2_b, x_temp3_b, q_temp1_b, q_temp2_b, ///
								   y_mat_rep_b, GMM_est_2_b)
		// Dynamic, kink - (Case 2)
	}
	else if (flag_static == 1 & flag_kink == 1) {
		ep_collect_2_b = ep_cal_kink(r_tt, x_temp1_c_b, x_temp2_c_b, x_temp3_c_b, q_temp1_b, q_temp2_b, ///
								   y_mat_rep_b, GMM_est_2_b)
		// Static, kink - (Case 3)
	}
	else {
		ep_collect_2_b = ep_cal(r_tt, x_temp1_b, x_temp2_b, x_temp3_b, q_temp1_b, q_temp2_b, ///
							  y_mat_rep_b, GMM_est_2_b)
		// Dynamic, not kink (Case 4) (default)
	}



	omega_hat_b = GMM_W_n_2_con(ep_collect_2_b, Z_mat_b)
	//This is actually an estimate of omega^(-1). Note the diffrence(with MATLAB ver.)

	//Now we estimate covariance matrix for 4 cases

	if (flag_static == 1 & flag_kink == 0) {
		// (Case 1): Static model, not kink design
		V_s = .
		cov_con_b(GMM_est_2_b, r_tt, Z_mat_b, x_temp1_c_b, x_temp2_c_b, ///
					x_temp3_c_b, q_temp1_b, q_temp2_b, omega_hat_b, V_s)
	}
	else if (flag_static == 0 & flag_kink == 1) {
		// (Case 2): Dynamic model, with kink 
		V_s = .
		cov_con_kink_b(GMM_est_2_b, r_tt, Z_mat_b, x_temp1_b, x_temp2_b, x_temp3_b, ///
						 q_temp1_b, q_temp2_b, omega_hat_b, V_s)
	}
	else if (flag_static == 1 & flag_kink == 1) {
		// (Case 3): Static model, with kink 
		V_s = .
		cov_con_kink_b(GMM_est_2_b, r_tt, Z_mat_b, x_temp1_c_b, x_temp2_c_b, ///
						 x_temp3_c_b, q_temp1_b, q_temp2_b, omega_hat_b, V_s)
	}
	else {
		// (Case 4): Dynamic model, not kink design (Default)
		V_s = .
		cov_con_b(GMM_est_2_b, r_tt, Z_mat_b, x_temp1_b, x_temp2_b, x_temp3_b, ///
					q_temp1_b, q_temp2_b, omega_hat_b, V_s)
	}
	

	size_c = cols(V_s) // size of (V_s' * V_s)
	dlt_hat = GMM_est_2_b[size_c-dlt_len+1 .. size_c]

	sel_m = J(dlt_len, size_c, 0) 
	sel_m[., size_c-dlt_len+1 .. size_c] = I(dlt_len)
	// To select delta_cov_mat from the whole (V_s' * V_s)
	// This is R mat. in the paper (172p).

	cov_wald = sel_m * cholinv(V_s' * V_s) * sel_m'
	wald_b =  N * dlt_hat' * cholinv(cov_wald) * dlt_hat

	wald_vec[bb] = wald_b
	
	cov_wald_rec[(bb-1)*dlt_len+1 .. bb*dlt_len, .] = cov_wald

	}

	sup_wald = max(wald_vec)
	// this is the desired sup wald stat.
}


void sup_wald_fb(													 	 ///
						real matrix y_mat_rep_b, real matrix Z_mat_b,	 ///
						real matrix x_temp1_b, real matrix x_temp2_b, 	 ///
						real matrix x_temp3_b,						 ///
						real matrix q_temp1_b, real matrix q_temp2_b, ///
						real matrix ep_collect_2, real matrix W_n_2, ///
						real matrix cov_wald_rec, 					 ///
						real matrix moment_mat,	real matrix grid_th, ///
						real scalar grid_num,						 ///
						real scalar flag_static,					 /// 
						real scalar flag_kink,						 ///
						real scalar N, real scalar T,				 ///
						sup_wald_b)
{	
	// ***(Step 2): 1st-step GMM estimate*** //
	k_1 = cols(x_temp1_b)/N
	L = rows(W_n_2)
	
	// We separate delta_hat from the GMM_estimate
	// Note that delta_hat includes constant difference
	if (flag_static == 1 & flag_kink == 0) {
		// (Case 1): Static model, not kink design
		dlt_len = k_1
	}
	else if (flag_static == 0 & flag_kink == 1) {
		// (Case 2): Dynamic model, with kink 
		dlt_len = 2
	}
	else if (flag_static == 1 & flag_kink == 1) {
		// (Case 3): Static model, with kink 
		dlt_len = 2
	}
	else {
		// (Case 4): Dynamic model, not kink design (Default)
		dlt_len = k_1 + 1
	}
	
	
	norm_vec = invnormal(uniform(1, N))
	b_eta = J(L, 1, norm_vec) // For bootstrap weighting
	g_1n_bar_b = rowsum(Z_mat_b:*ep_collect_2:*b_eta)/N 
	//This part makes difference for Fast Bootstrap (vs. original)
	

	if (flag_static == 1){
		end_p = cols(x_temp1_b);
		
		x_temp1_c_b = x_temp1_b[., N+1..end_p]
		x_temp2_c_b = x_temp2_b[., N+1..end_p]
		x_temp3_c_b = x_temp3_b[., N+1..end_p]
		// Copies of matrices, after removing columns that contain
		// y_{i, t-1} related terms
		
		// These matrices will be used if we assume static model
	}


	//Now with the given threshold, we calculate 1st-step GMM estimator
	//In the following for-loop we will calculate sup Wald.

	wald_vec = J(grid_num, 1, 0)

	for (bb=1; bb<=grid_num; bb=bb+1) {
	// Specially I omit indentation for this loop (too long!)

	r_tt = grid_th[bb]

	g_2n_bar_b = .
	g_n_bar_b = .
	GMM_est_b = .
			
	if (flag_static == 1 & flag_kink == 0) {
		// (Case 1): Static model, not kink design
		GMM_cal(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_c_b, x_temp2_c_b, x_temp3_c_b, ///
				   q_temp1_b, q_temp2_b, W_n_2, g_2n_bar_b, g_n_bar_b, GMM_est_b)
	}
	else if (flag_static == 0 & flag_kink == 1) {
		// (Case 2): Dynamic model, with kink
		GMM_cal_kink(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_b, x_temp2_b, x_temp3_b, ///
				   q_temp1_b, q_temp2_b, W_n_2, g_2n_bar_b, g_n_bar_b, GMM_est_b)
	}
	else if (flag_static == 1 & flag_kink == 1) {
		// (Case 3): Static model, with kink 
		GMM_cal_kink(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_c_b, x_temp2_c_b, x_temp3_c_b, ///
						q_temp1_b, q_temp2_b, W_n_2, g_2n_bar_b, g_n_bar_b, GMM_est_b)
	}
	else {
		// (Case 4): Dynamic model, not kink design (Default)
		GMM_cal(g_1n_bar_b, r_tt, Z_mat_b, x_temp1_b, x_temp2_b, x_temp3_b, ///
				q_temp1_b, q_temp2_b, W_n_2, g_2n_bar_b, g_n_bar_b, GMM_est_b)
	}

	//This is the 1st-step bootstrap GMM estimate(GMM_est_b)
	//But in Fast bootstrap scheme, 2nd step is NOT NEEDED (W_n_2 is given)
	//So we omit the 2nd-step GMM
	//Moreover we don't need to calculate cov_wald (it is also given)

	
	size_c = rows(GMM_est_b)
	dlt_hat = GMM_est_b[size_c-dlt_len+1 .. size_c]
	
	cov_wald = cov_wald_rec[(bb-1)*dlt_len+1 .. bb*dlt_len, .]
	// This is the same as the cov_wald's for original wald stat.
	
	wald_b =  N * dlt_hat' * cholinv(cov_wald) * dlt_hat

	wald_vec[bb] = wald_b
	}

	sup_wald_b = max(wald_vec) // this is the desired sup wald stat.
	
}


void cov_con_b(													///
					real matrix GMM_est_2, real scalar r_hat_2,	///
					real matrix Z_mat, real matrix x_temp1, 	///
					real matrix x_temp2, real matrix x_temp3,	///
					real matrix q_temp1, real matrix q_temp2,	///
					real matrix omega_hat,		///
					V_s)
{
	//This function estimates covariance matrix for bootstrapping section.
	//Refer to (173p.) for \hat{V}_s
	
	L = rows(Z_mat)
	N = cols(Z_mat)
	//h_band = h_0 * 1.06 * N^(-0.2) * sqrt(variance(colshape(q_mat, 1))) // Not needed here
	k_1 = cols(x_temp1)/N
	
	G_b = J(L, k_1, 0)
	G_d = J(L, k_1 + 1, 0)
	//G_r = J(L, 1, 0) // Not needed here
	//Three matrices are \hat{G_b}, \hat{G_d}, \hat{G_r} on Seo & Shin(2016), 173p. 
	
	ind_temp2 = (q_temp1 :> r_hat_2)
	ind_temp3 = (q_temp2 :> r_hat_2) 
	// This two (L*N) 0-1 matrices contain information whether q_{it} or q_{it-1} is 
	// larger than the threshold r_th
	
	g_2n_bar = J(L, 2*k_1 + 1, 0) //This is identical to g_2n_bar in Seo & Shin(2016)
	g_2n_bar[., k_1 + 1] = rowsum(Z_mat:*(ind_temp2 - ind_temp3))/N
	

	for (ii6=1; ii6<=k_1; ii6=ii6+1) {
		temp_x = x_temp1[., (ii6-1)*N+1..ii6*N] //Extract (L*N) for one variable
		G_b[., ii6] = -rowsum(Z_mat:*x_temp1[., (ii6-1)*N+1 .. ii6*N])/N
		
		g_2n_bar[., ii6] = rowsum(Z_mat:*temp_x)/N
		g_2n_bar[., k_1 + 1 + ii6] = rowsum(Z_mat:*(x_temp2[., (ii6-1)*N+1..ii6*N]:*ind_temp2 ///
											- x_temp3[., (ii6-1)*N+1..ii6*N]:*ind_temp3))/N
	}
	
	G_d[., 1] = -g_2n_bar[., k_1 + 1]
	G_d[., 2..k_1 + 1] = -g_2n_bar[., k_1 + 2..2*k_1 + 1]
	
	V_s = cholesky(omega_hat)' * (G_b, G_d)
}


void cov_con_kink_b(													///
						real matrix GMM_est_2, real scalar r_hat_2,	///
						real matrix Z_mat, real matrix x_temp1, 	///
						real matrix x_temp2, real matrix x_temp3,	///
						real matrix q_temp1, real matrix q_temp2,	///
						real matrix omega_hat,		///
						V_s)
{
	//This function estimates covariance matrix for bootstrapping section. (w.r.t kink model)
	//Refer to (173p.) for \hat{V}_s
	
	L = rows(Z_mat)
	N = cols(Z_mat)
	//h_band = h_0 * 1.06 * N^(-0.2) * sqrt(variance(colshape(q_mat, 1))) // Not needed here
	k_1 = cols(x_temp1)/N
	
	G_b = J(L, k_1, 0)
	G_d = J(L, 1, 0)
	//G_r = J(L, 1, 0) // Not needed here
	//Three matrices are \hat{G_b}, \hat{G_d}, \hat{G_r} on Seo & Shin(2016), 173p. 
	
	ind_temp2 = (q_temp1 :> r_hat_2)
	ind_temp3 = (q_temp2 :> r_hat_2) 
	// This two (L*N) 0-1 matrices contain information whether q_{it} or q_{it-1} is 
	// larger than the threshold r_th
	
	G_d[., 1] = rowsum(Z_mat:*((q_temp1:-r_hat_2):*ind_temp2 - ///
							   (q_temp2:-r_hat_2):*ind_temp3))/N
	// Adjusted for kink model
	
	for (ii6=1; ii6<=k_1; ii6=ii6+1) {
		temp_x = x_temp1[., (ii6-1)*N+1..ii6*N] //Extract (L*N) for one variable
		G_b[., ii6] = -rowsum(Z_mat:*temp_x)/N
	}
	
	V_s = cholesky(omega_hat)' * (G_b, G_d)
}
end
