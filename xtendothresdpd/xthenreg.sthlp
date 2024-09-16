{smcl}
{* April 05, 2019}{...}
{cmd:help xthenreg}{right: ({browse "https://doi.org/10.1177/1536867X19874243":SJ19-3: st0573})}
{hline}

{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{cmd:xthenreg} {hline 2}}Dynamic panel-data model allowing threshold
and endogeneity (regression)


{title:Syntax}

{p 8 16 2}
{cmd:xthenreg} {depvar} {it:indepvars} {ifin} [{cmd:,} {it:options}]

{synoptset 25}{...}
{synopthdr :options}
{synoptline}
{synopt:{opt static}}specify a static model; default model is dynamic; 
unlike dynamic model, static model does not automatically include
{cmd:L.y} as independent variable{p_end}
{synopt:{opt kink}}specify a kink model{p_end}
{synopt:{opt endo:genous(varlist)}}specify endogenous independent variables; it must be excluded from the list of independent
variables before comma{p_end}
{synopt:{opt inst(varlist)}}specify the list of additional instrumental variables{p_end}
{synopt:{opt grid_num(integer)}}determine the number of grid points to estimate the
threshold r; default is {cmd:grid_num(100)}{p_end}
{synopt:{opt trim_rate(real)}}determine the trim rate when constructing a grid
for estimating r; default is {cmd:trim_rate(0.4)}{p_end}
{synopt:{opt h_0(real)}}determine a parameter for Silverman's rule of thumb
used for kernel estimation; default is {cmd:h_0(1.5)}{p_end}
{synopt:{opt boost(integer)}}number of bootstrapping for linearity test; default
is {cmd:boost(0)}{p_end}
{synoptline}
{p2colreset}{...}
{phang}
Pressing {cmd:Break} in the middle of the execution may not restore the original dataset.{p_end}
{phang}
{cmd:xtset} must be executed to declare data to be panel data before using this command.
{p_end}
{phang}
The {cmd:moremata} package is required to run this command.  Moreover, data
have to be strongly balanced panel data.  Both requirements are automatically checked.{p_end} 


{title:Description}

{p 4 4 2}
{cmd:xthenreg} estimates coefficients of a dynamic panel model with threshold
and endogeneity.  The default model is a dynamic threshold model, but options
can change the model form.  The estimator proposed by Seo and Shin (2016) is
implemented.  The result includes point estimates and asymptotic standard
errors.  The threshold r in the model is estimated by grid search to minimize
the objective function of the generalized method of moments.

{title:Examples}

{p 4 8 2}{stata "xthenreg lggdppccstd lginvestgdpr, endo(lginvestgdpr) inst(lginvestgdpr) kink static"}{p_end}

{p 4 8 2}{stata "xthenreg lggdppccstd lginvestgdpr, endo(lginvestgdpr) exo(period) inst(lginvestgdpr) kink static"}{p_end}

{p 4 8 2}{stata "xthenreg lggdppccstd lginvestgdpr, endo(lginvestgdpr) exo(period) inst(lginvestgdpr) static"}{p_end}

{title:Options}

{phang}
{opt static} specifies a static model; the default model is dynamic.  Unlike
the dynamic model, the static model does not automatically include {cmd:L.y}
as an independent variable.

{phang}
{opt kink} specifies a kink model.

{phang}
{opt endogenous(varlist)} specifies endogenous independent variables.  The
endogenous variables must be excluded from the list of independent variables
before the comma.  For example, if {cmd:x1} and {cmd:x2} are independent
variables and {cmd:x2} is endogenous, you must use a command like
{cmd:xthenreg y q x1, endogenous(x2)}, not 
{cmd:xthenreg y q x1 x2, endogenous(x2)}.

{phang}
{opt inst(varlist)} specifies the list of additional instrumental variables.

{phang}
{opt grid_num(integer)} determines the number of grid points to estimate the
threshold r.  The default is {cmd:grid_num(100)}.

{phang}
{opt trim_rate(real)} determines the trim rate when constructing a grid for
estimating r.  It must be a positive real number smaller than 1.  For example,
if the trim rate is 0.2, the start and end points of the grid are the 0.1
quantile and 0.9 quantile of the observed variable q, respectively.  The
default is {cmd:trim_rate(0.4)}.

{phang}
{opt h_0(real)} determines a parameter for Silverman's rule of thumb, which is
used to determine the bandwidth of the Nadaraya-Watson kernel for covariance
matrix estimation.  It must be a positive real number.  The default is
{cmd:h_0(1.5)}.

{phang}
{opt boost(integer)} is the number of bootstrapping iterations for the
linearity test.  The default is {cmd:boost(0)}.  If you give a positive
integer for this option, the bootstrapping linearity test in Seo and Shin
(2016) is executed.  However, this test requires a lot of computation time.


{title:Stored results}

{pstd}
{cmd:xthenreg} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of units of panel data{p_end}
{synopt:{cmd:e(T)}}time length of panel data{p_end}
{synopt:{cmd:e(boots_p)}}p-value for bootstrap linearity test; {cmd:-1} if the test
is not used{p_end}
{synopt:{cmd:e(grid)}}number of grid points used{p_end}
{synopt:{cmd:e(trim)}}trim rate for grid search{p_end}
{synopt:{cmd:e(bs)}}number of bootstrap iterations{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvar)}}name of independent variable or variables{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(zx)}}name of instrumental variables{p_end}
{synopt:{cmd:e(qx)}}name of the threshold variable{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}estimates of coefficients{p_end}
{synopt:{cmd:e(V)}}estimate of the covariance matrix{p_end}
{synopt:{cmd:e(CI)}}95% asymptotic confidence interval for {cmd:b}{p_end}


{title:Package requirements}

{p 4 4 2}
The {cmd:moremata} package is required to calculate the quantile of data.


{title:Methods and formulas}

{pstd}
See Seo and Shin (2016).


{title:Reference}

{phang}
Seo, M. H., and Y. Shin. 2016. Dynamic panels with threshold effect and
endogeneity. {it:Journal of Econometrics} 195: 169-186.


{title:Disclaimer}

{p 4 4 2}
THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
PROGRAM IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST
OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

{p 4 4 2}
IN NO EVENT WILL THE COPYRIGHT HOLDERS OR THEIR EMPLOYERS, OR ANY OTHER PARTY
WHO MAY MODIFY AND/OR REDISTRIBUTE THIS SOFTWARE, BE LIABLE TO YOU FOR
DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM.


{title:Authors}

{pstd}
Myung Hwan Seo{break}
Department of Economics & Institute of Economic Research{break}
Seoul National University{break}
Seoul, Korea{break}
myunghseo@snu.ac.kr

{pstd}
Sueyoul Kim{break}
Department of Economics{break}
University of Maryland{break}
College Park, MD{break}
shelkim@econ.umd.edu

{pstd}
Young-Joo Kim{break}
Department of Economics{break}
Hongik University{break}
Seoul, Korea{break}
y.j.kim@hongik.ac.kr


{marker alsosee}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 19, number 3: {browse "https://doi.org/10.1177/1536867X19874243":st0573}{p_end}
