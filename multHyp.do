/* multHyp.do                    damiancclarke             yyyy-mm-dd:2021-11-10
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

  This code runs some simulations and conducts multiple hypothesis test correct-
ions for a talk at the Stata Economics Symposium:
     "Multiple Hypothesis Testing in Stata"
       by: Damian Clarke

  This code replicates results shown in the slides using a range of procedures
correcting for multiple testing in a number of ways.  For full details please
refer to the slides available with this code.  Many of these procedures will
take a considerable amount of time to run given that they are based on many
simulations with bootstrap iterations.

  This code requires a number of user-written ados which must be installed from
the SSC or Stata Journal's repository.  If these are not available, they will
be installed on lines 37-47.  If net aware Stata is not used and these ados are
not installed, the code will fail to run.

  If you wish to explore faster (though less precise) versions of these simulat-
ions, the locals defined as S below could be set at lower values. Throughout the
code, the local S sets the number of simulations, while the local B sets the nu-
mber of boottsrap replicates.

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
local ados graph3d swindex rwolf2 rwolf 
foreach ado of local ados  {
    cap which `ado'
    if _rc!=0 ssc install `ado'               
}
cap which multproc
if _rc!=0 {
    net from http://www.stata-journal.com/software/sj10-4
    net install st0035_1
}
cap which plotplainblind
if _rc!=0 ssc install blindschemes

cap mkdir results
set scheme plotplainblind

*-------------------------------------------------------------------------------
*--- (1) Demonstration of over rejection (class wise) when testing multiple hyps
*-------------------------------------------------------------------------------
set seed 1213
local S = 5000

**Consider up to K=10 depdendent variables (Table on slide 2)
foreach K of numlist 1(1)10 {
    local reject`K' = 0
    local prop`K'   = 0
    local prop2`K'  = 0
    forvalues s = 1/`S' {
        set obs 100
        gen treat = runiform()>0.5
        local rejectAny = 0
        foreach k of numlist 1(1)`K' {
            gen y`k' = 2 + 0*treat + rnormal()
            reg y`k' treat        
            if abs(_b[treat]/_se[treat])>invnormal(0.975) {
                local ++reject`K'
                local ++rejectAny
            }
        }
        dis "Reject any is `rejectAny'"
        if `rejectAny'>0 local ++prop`K'
        if `rejectAny'>1 local ++prop2`K'        
        clear
    }
}
local TTR "Total Tests Rejected "
local MTR "Mean Tests Rejected"
local PTR "Proportion $\geq$ 1 Rejection"
local P2R "Proportion $\geq$ 2 Rejection"

**Export output for Table on slide 2
foreach num of numlist 1(1)10 {
    local ttr = `reject`num''
    local mtr = string(`reject`num''/(`num'*`S'), "%04.3f")
    local ptr = string(`prop`num''/`S', "%04.3f")
    local p2r = string(`prop2`num''/`S', "%04.3f")
    
    dis "Total rejected (`num') is " `reject`num''
    dis "Rejected (`num') is " `reject`num''/(`num'*`S')
    dis "Proportion rejected is " `prop`num''/`S'

    local TTR "`TTR' & `ttr'"
    local MTR "`MTR' & `mtr'"
    local PTR "`PTR' & `ptr'"
    local P2R "`P2R' & `p2r'"
}
dis "`TTR'"
dis "`MTR'"
dis "`PTR'"
dis "`P2R'"


*-------------------------------------------------------------------------------
*--- (2) Demonstration of mechanics which will be used below (MVN simulation)
*-------------------------------------------------------------------------------
foreach x of numlist 0 50 99 {
    clear
    set obs 1000
    local r=`x'/100
    gen treat  = runiform()>0.5
    mat corr = (1, `r', `r' \ `r', 1, `r' \ `r', `r', 1)

    **use Cholesky and matrix score for this (eg Gould's FAQ)
    mat corsim = cholesky(corr)    
    foreach num of numlist 1(1)3 {
        gen c`num'     = rnormal()
    }
    foreach num of numlist 1(1)3 {
        matrix eps`num' = corsim[`num',1..3]
        matrix score epsilon`num' = eps`num'
        gen y`num'  = 1 + 0*treat   + epsilon`num'
    }
    
    #delimit ;
    graph3d y1 y2 y3, ycam(-4) zcam(-18) cuboid innergrid blv
    perspective colorscheme(bcgyr) xangle(80) yangle(179) zangle(45);
    #delimit cr
    graph export "./results/correlationY_`x'.eps", replace
}
clear

*-------------------------------------------------------------------------------
*--- (3) Indexes: Anderson, 1st Principal Component, Summary Index
*-------------------------------------------------------------------------------
local S = 1000
set obs `S'
foreach index in Anderson PC Sum {
    gen `index'_t  = .
    gen `index'B_t = .
    gen `index'C_t = .
    gen `index'2_t = .
}

** Run S simulation testing for real impacts depending on correlational structure
** Slides 7-9
forvalues s=1/`S' {
    dis "Simulation Number `s'"
    preserve
    clear
    set obs 1000
    local r = 0.9
    #delimit ;
    mat corr =
       (1,`r',`r',`r',`r',`r',`r',`r',`r',`r' \
        `r',1,`r',`r',`r',`r',`r',`r',`r',`r' \
        `r',`r',1,`r',`r',`r',`r',`r',`r',`r' \
        `r',`r',`r',1,`r',`r',`r',`r',`r',`r' \
        `r',`r',`r',`r',1,`r',`r',`r',`r',`r' \
        `r',`r',`r',`r',`r',1,`r',`r',`r',`r' \
        `r',`r',`r',`r',`r',`r',1,`r',`r',`r' \
        `r',`r',`r',`r',`r',`r',`r',1,`r',`r' \
        `r',`r',`r',`r',`r',`r',`r',`r',1,`r' \
        `r',`r',`r',`r',`r',`r',`r',`r',`r',1 );
    #delimit cr
    mat corsim = cholesky(corr)
    
    foreach num of numlist 1(1)10 {
        gen c`num'     = rnormal()
    }
    gen treat = runiform()<0.5
    foreach num of numlist 1(1)10 {
        matrix eps`num' = corsim[`num',1..10]
        matrix score epsilon`num' = eps`num'
        gen y`num'  = 1 + 0*treat   + epsilon`num'
        gen yr`num' = 1 + 0.5*treat + epsilon`num'
    }
    gen yr11 = 1+0.5*treat + rnormal()

    gen control = treat==0

    **Anderson Index
    swindex y1 y2 y3 y4 y5 y6 y7 y8 y9 y10, normby(control) generate(AndersonIndex)
    swindex y1 y2 y3 y4 y5 y6 y7 y8 y9 yr10, normby(control) generate(AndersonIndexB)
    swindex y1 y2 y3 y4 y5 y6 y7 y8 y9 yr11, normby(control) generate(AndersonIndexC)
    swindex yr1 yr2 yr3 yr4 yr5 yr6 yr7 yr8 yr9 yr10, normby(control) generate(AndersonIndex2)

    **1st Principal Component
    pca y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 
    predict PCIndex, score
    pca y1 y2 y3 y4 y5 y6 y7 y8 y9 yr10
    predict PCIndexB, score
    pca y1 y2 y3 y4 y5 y6 y7 y8 y9 yr11
    predict PCIndexC, score
    pca yr1 yr2 yr3 yr4 yr5 yr6 yr7 yr8 yr9 yr10
    predict PCIndex2, score

    **Summary Index
    egen SumIndex  = rowtotal(y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 )
    egen SumIndexB = rowtotal(y1 y2 y3 y4 y5 y6 y7 y8 y9 yr10)
    egen SumIndexC = rowtotal(y1 y2 y3 y4 y5 y6 y7 y8 y9 yr11)
    egen SumIndex2 = rowtotal(yr1 yr2 yr3 yr4 yr5 yr6 yr7 yr8 yr9 yr10)

    foreach index in Anderson PC Sum {
        reg `index'Index treat
        local `index'I = _b[treat]/_se[treat]
        reg `index'Index2 treat
        local `index'I2 = _b[treat]/_se[treat]    
        reg `index'IndexB treat
        local `index'IB = _b[treat]/_se[treat]
        reg `index'IndexC treat
        local `index'IC = _b[treat]/_se[treat]
    }
    restore
    foreach index in Anderson PC Sum {
        replace `index'_t  = ``index'I'  in `s'
        replace `index'B_t = ``index'IB' in `s'
        replace `index'C_t = ``index'IC' in `s'
        replace `index'2_t = ``index'I2' in `s'
    }
}
**Output results
foreach index in Anderson PC Sum {
    foreach var of varlist `index'_t `index'_t `index'C_t `index'2_t {
        #delimit ;
        twoway hist `var', xtitle("{it:t}-statistic") ylabel(, format(%3.1f))
          || function y=normalden(x), range(-3 3) lcolor(red) lpattern(solid)
        legend(off) xlabel(-3(1)11);
        #delimit cr
        graph export "./results/Index_`var'.eps", replace
    }
}
clear

*-------------------------------------------------------------------------------
*--- (4) FWER tests (Naive, Bonferroni, Holm, Romano-Wolf)
*-------------------------------------------------------------------------------
**These simulations take a long time to get smoothish power graphs, but can be
**seen iteratively...
cap log close
log using "./results/FWERtests.txt", replace
foreach rho in 25 75 {
    set seed 12011303
    local S=50
    local B=250


    set obs 20
    gen beta  = (_n/100)*5-0.04
    gen probNA = .
    gen probRW = .
    gen probHO = .
    gen probBO = .
    local k = 1
    foreach beta of numlist 0.01(0.05)1 {
        dis "\beta is: `beta'"
        local propNA = 0
        local propHO = 0
        local propBO = 0
        local propRW = 0
    
        preserve
        forval s=1/`S' {
            dis "SIM NUMBER `s'"
            clear all
            set obs 100
            gen treat  = runiform()>0.5
            foreach num of numlist 1(1)10 {
                gen c`num'     = rnormal()
            }
            local r=`rho'/100
            #delimit ;
            mat corr = (1,`r',`r',`r',`r',`r',`r',`r',`r',`r' \
                `r',1,`r',`r',`r',`r',`r',`r',`r',`r' \
                `r',`r',1,`r',`r',`r',`r',`r',`r',`r' \
                `r',`r',`r',1,`r',`r',`r',`r',`r',`r' \
                `r',`r',`r',`r',1,`r',`r',`r',`r',`r' \
                `r',`r',`r',`r',`r',1,`r',`r',`r',`r' \
                `r',`r',`r',`r',`r',`r',1,`r',`r',`r' \
                `r',`r',`r',`r',`r',`r',`r',1,`r',`r' \
                `r',`r',`r',`r',`r',`r',`r',`r',1,`r' \
                `r',`r',`r',`r',`r',`r',`r',`r',`r',1);
            #delimit cr
            mat corsim = cholesky(corr)
    
    
            foreach num of numlist 1(1)10 {
                matrix eps`num' = corsim[`num',1..10]
                matrix score epsilon`num' = eps`num'
                gen y`num'  = 1 + `beta'*treat   + epsilon`num'
            }
            rwolf y1 y2 y3 y4 y5 y6 y7 y8 y9 y10, indepvar(treat) reps(`B') holm
            matrix pvals = e(RW)
            mat list pvals
            foreach num of numlist 1(1)10 {
                local pna = pvals[`num',2]
                local bon = `pna'*10
                local prw = pvals[`num',3]
                local pho = pvals[`num',4]

                dis `pna'
                
                if `pna'<=0.05 local ++propNA
                if `bon'<=0.05 local ++propBO
                if `prw'<=0.05 local ++propRW
                if `pho'<=0.05 local ++propHO
            }
        }
        restore
        dis "when beta is `beta', proportion Naive rejected is `=`propNA'/`S''"
        dis "when beta is `beta', proportion Bonferroni rejected is `=`propBO'/`S''"
        dis "when beta is `beta', proportion Holm rejected is `=`propHO'/`S''"
        dis "when beta is `beta', proportion R-W rejected is `=`propRW'/`S''"
    
        replace probNA = `propNA'/(`S'*10) in `k'
        replace probBO = `propBO'/(`S'*10) in `k'
        replace probHO = `propHO'/(`S'*10) in `k'
        replace probRW = `propRW'/(`S'*10) in `k'
        sum prob*
        local ++k
        
        #delimit ;
        twoway connected probNA beta
        ||  connected probBO beta
        ||  connected probHO beta
        ||  connected probRW beta,
        ytitle("Proportion of Nulls Correctly Rejected")
        xtitle("Value for {&beta}")
        legend(lab(1 "Uncorrected") lab(2 "Bonferroni")
               lab(3 "Holm") lab(4 "Romano-Wolf")
               rows(1) position(12));
        graph export "./results/powerExample_`rho'.eps", replace;
        #delimit cr
        save "./results/simulationResults_`rho'.dta", replace
    }
    clear
}
log close

*-------------------------------------------------------------------------------
*--- (5) Comparison of FDR with FWER
*-------------------------------------------------------------------------------
foreach rho of numlist 0 33 67 {
    clear
    log using ./results/allTest_`rho'.txt, replace text
    set obs 1000
    foreach correc in Naive BH BY Bonf Holm WY RW {
        gen FP_`correc' = .
        gen TP_`correc' = .
    }

    local S=500
    forvalues s=1/`S' {
        dis "Simulation Number `s'"
        preserve
        clear
        set obs 100
        local r = `rho'/100
        #delimit ;
        mat corr =
           (1,`r',`r',`r',`r',`r',`r',`r',`r',`r' \
            `r',1,`r',`r',`r',`r',`r',`r',`r',`r' \
            `r',`r',1,`r',`r',`r',`r',`r',`r',`r' \
            `r',`r',`r',1,`r',`r',`r',`r',`r',`r' \
            `r',`r',`r',`r',1,`r',`r',`r',`r',`r' \
            `r',`r',`r',`r',`r',1,`r',`r',`r',`r' \
            `r',`r',`r',`r',`r',`r',1,`r',`r',`r' \
            `r',`r',`r',`r',`r',`r',`r',1,`r',`r' \
            `r',`r',`r',`r',`r',`r',`r',`r',1,`r' \
            `r',`r',`r',`r',`r',`r',`r',`r',`r',1 );
        #delimit cr
        mat corsim = cholesky(corr)
    
        foreach num of numlist 1(1)10 {
            gen c`num'     = rnormal()
        }
        gen treat = runiform()<0.5
        foreach num of numlist 1(1)10 {
            matrix eps`num' = corsim[`num',1..10]
            matrix score epsilon`num' = eps`num'
            gen y`num'  = 1 + 0*treat   + epsilon`num'
            gen yr`num' = 1 + 0.5*treat + epsilon`num'
        }
    
        **Outcomes: y1, y2, y3, y4, y5, yr6, yr7, yr8, yr9, yr10
        gen p = .
        local i = 1
        foreach y of varlist y1 y2 y3 y4 y5 y6 y7 yr8 yr9 yr10 {
        qui reg `y' treat
        qui test treat
        qui replace p=r(p) in `i'
        local ++i
        } 
        **Naive
        count if p<0.05 in 1/7
        local FP_Naive = r(N)
        count if p<0.05 in 8/10
        local TP_Naive = r(N)
        **Benjamini-Hochberg
        multproc, pvalue(p) reject(pBH) method(simes)
        count if pBH==1 in 1/7
        local FP_BH = r(N)
        count if pBH==1 in 8/10
        local TP_BH = r(N)
        **Benjamini-Yekuetieli
        multproc, pvalue(p) reject(pBY) method(yekutieli)
        count if pBY==1 in 1/7
        local FP_BY = r(N)
        count if pBY==1 in 8/10
        local TP_BY = r(N)
        **Bonferroni
        multproc, pvalue(p) reject(pBonf) method(bonferroni)
        count if pBonf==1 in 1/7
        local FP_Bonf = r(N)
        count if pBonf==1 in 8/10
        local TP_Bonf = r(N)
        **Holm
        multproc, pvalue(p) reject(pHolm) method(holm)
        count if pHolm==1 in 1/7
        local FP_Holm = r(N)
        count if pHolm==1 in 8/10
        local TP_Holm = r(N)
        #delimit ;
        wyoung y1 y2 y3 y4 y5 y6 y7 yr8 yr9 yr10, cmd(regress OUTCOMEVAR treat)
        bootstraps(56) familyp(treat);
        #delimit cr
        mat WY=r(table)
        svmat WY
        count if WY4<0.05 in 1/7
        local FP_WY = r(N)
        count if WY4<0.05 in 8/10
        local TP_WY = r(N)
        #delimit ;
        rwolf y1 y2 y3 y4 y5 y6 y7 yr8 yr9 yr10, indepvar(treat) reps(56)
        noplusone nodots;
        #delimit cr
        matrix RW=e(RW)	
        svmat RW
        count if RW3<0.05 in 1/7
        local FP_RW = r(N)
        count if RW3<0.05 in 8/10
        local TP_RW = r(N)
    
        restore
        foreach correc in Naive BH BY Bonf Holm WY RW {
            replace FP_`correc' = `FP_`correc'' in `s'
            replace TP_`correc' = `TP_`correc'' in `s'       
        }
    }
    sum FP*
    sum TP*
    foreach correc in Naive BH BY Bonf Holm WY RW {
        qui count if FP_`correc'>0&FP_`correc'<=8
        local PrA = string(r(N)/`S', "%05.3f")
        gen FDR_`correc' = FP_`correc'/(FP_`correc'+TP_`correc')
        qui sum FDR_`correc'
        local PrB = string(r(mean), "%05.3f")
        qui sum TP_`correc'
        local PrC = string(r(mean)/8, "%05.3f")

        dis "`correc' & `PrA' & `PrB' & `PrC' "
    }
    log close    
}


*-------------------------------------------------------------------------------
*--- (6) Some other things
*-------------------------------------------------------------------------------
foreach rho of numlist 0 50 90 {
    clear
    set obs 1000
    local r = `rho'/100
    #delimit ;
    mat corr =
       (1,`r',`r',`r',`r' \
        `r',1,`r',`r',`r' \
        `r',`r',1,`r',`r' \
        `r',`r',`r',1,`r' \
        `r',`r',`r',`r',1 );
    #delimit cr
    mat corsim = cholesky(corr)
    foreach num of numlist 1(1)5 {
        gen c`num'     = rnormal()
    }
    gen treat = runiform()<0.5
    foreach num of numlist 1(1)5 {
        matrix eps`num' = corsim[`num',1..5]
        matrix score epsilon`num' = eps`num'
        gen y`num'  = 1 + ((`num'-1)/10)*treat   + epsilon`num'
    }
    #delimit ;
    rwolf2 (reg y1 treat) (reg y2 treat) (reg y3 treat) (reg y4 treat) (reg y5 treat),
    indepvars(treat, treat, treat, treat, treat) reps(500) nodots graph;
    graph export "./results/RW_rho`rho'.eps", replace;
    #delimit cr
}



**Examining FDR with many more variables
**This can be set at maximum 1000 variables given observation numbers 
local numTests = 200
forvalues num = 6/`numTests' {
    gen y`num'= rnormal()
}

**FDR rates
gen p = .
forvalues num = 1/`numTests' {
    reg y`num' treat
    test treat
    replace p=r(p) in `num'
}
andersonFDR p
multproc, pvalue(p) reject(pBH) method(simes)
multproc, pvalue(p) reject(pBY) method(yekutieli)
rwolf y*, indepvar(treat)




