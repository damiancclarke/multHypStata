/* singleExample.do              damiancclarke             yyyy-mm-dd:2022-04-26
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

  This code implements a single example of a table exporting regression results
with multiple hypothesis testing corrections based on a simulated example, for
a talk at University of Wisconsin-Madison:
     "Correcting for Multiple Hypothesis Testing"
       by: Damian Clarke

  This code requires a number of user-written ados which must be installed from
the SSC or Stata Journal's repository.  If these are not available, they will
be installed on lines 37-47.  If net aware Stata is not used and these ados are
not installed, the code will fail to run.

  Contact email: dclarke@fen.uchile.cl

*/

*-------------------------------------------------------------------------------
*--- (0) General set-up
*-------------------------------------------------------------------------------
vers 14
clear all
set more off
cap log close

**Check for all necesary user-written ados, install if not available
local ados swindex rwolf2 estout
foreach ado of local ados  {
    cap which `ado'
    if _rc!=0 ssc install `ado'               
}


*-------------------------------------------------------------------------------
*--- (1) Simulate some data: 6 outcomes, 3 significant, 3 insignificant
*-------------------------------------------------------------------------------
set seed 156
set obs 1000
local r = 0.3
#delimit ;
mat corr = (1,`r',`r',`r',`r',`r' \
            `r',1,`r',`r',`r',`r' \
            `r',`r',1,`r',`r',`r' \
            `r',`r',`r',1,`r',`r' \
            `r',`r',`r',`r',1,`r' \
            `r',`r',`r',`r',`r',1 );
#delimit cr
mat corsim = cholesky(corr)
    
foreach num of numlist 1(1)6 {
    gen c`num'     = rnormal()
}
gen treat = runiform()<0.5
foreach num of numlist 1(1)6 {
    matrix eps`num' = corsim[`num',1..6]
    matrix score epsilon`num' = eps`num'

    if `num'<=3 local tau=0.12
    if `num'>3  local tau=0
    
    gen y`num'  = 1 + `tau'*treat   + epsilon`num'
    dis "Seed is `seed'"
}
drop c* epsilon*

gen control = treat==0
swindex y1 y2 y3 y4 y5 y6, normby(control) generate(AndersonIndex)
    
*-------------------------------------------------------------------------------
*--- (2) Conduct regressions (uncorrected)
*-------------------------------------------------------------------------------
eststo: reg y1 treat, nohead
eststo: reg y2 treat, nohead
eststo: reg y3 treat, nohead
eststo: reg y4 treat, nohead
eststo: reg y5 treat, nohead
eststo: reg y6 treat, nohead
eststo: reg AndersonIndex treat, nohead

#delimit ;
rwolf2 (reg y1 treat) (reg y2 treat) (reg y3 treat)
       (reg y4 treat) (reg y5 treat) (reg y6 treat),
indepvars(treat, treat, treat, treat, treat, treat) reps(1000) holm; 
#delimit cr

foreach num of numlist 1(1)6 {
    local pv`num' = e(rw_y`num'_treat)
}
foreach num of numlist 1(1)6 {
    estadd scalar pRW = `pv`num'' : est`num'
}

lab var treat "Treatment"
#delimit ;
esttab est1 est2 est3 est4 est5 est6 est7 using correctedRegression.tex,
replace booktabs cells(b(fmt(%-9.3f)) se(fmt(%-9.3f) par([ ]) )) label
stats(pRW N r2, fmt(%05.3f %9.0gc %5.3f)
      label("Romano Wolf p-value" Observations R-Squared))
mlabels("y1" "y2" "y3" "y4" "y5" "y6" "Index") style(tex)
collabels(none);

esttab est1 est2 est3 est4 est5 est6 est7 using correctedRegression_version2.tex,
replace booktabs cells(b(fmt(%-9.3f) star) se(fmt(%-9.3f) par([ ]) )) label
stats(pRW N r2, fmt(%05.3f %9.0gc %5.3f)
      label("Romano Wolf p-value" Observations R-Squared))
mlabels("y1" "y2" "y3" "y4" "y5" "y6" "Index") style(tex)
collabels(none) starlevel("*" 0.1 "**" 0.05 "***" 0.01);
#delimit cr

