* This code generates BKY (2006) sharpened two-stage q-values as described in Anderson (2008), "Multiple Inference and Gender Differences in the Effects of Early Intervention: A Reevaluation of the Abecedarian, Perry Preschool, and Early Training Projects", Journal of the American Statistical Association, 103(484), 1481-1495

* BKY (2006) sharpened two-stage q-values are introduced in Benjamini, Krieger, and Yekutieli (2006), "Adaptive Linear Step-up Procedures that Control the False Discovery Rate", Biometrika, 93(3), 491-507

* Last modified: M. Anderson, 11/20/07
* Test Platform: Stata/MP 10.0 for Macintosh (Intel 32-bit), Mac OS X 10.5.1
* Should be compatible with Stata 10 or greater on all platforms
* Likely compatible with with Stata 9 or earlier on all platforms (remove "version 10" line below)

**This has been slightly modified 10/11/2021 to accept a variable containing
**  p-values as arguments.  This is for illustrative purposes only for a talk
**  at Stata's economics symposium (Damian Clarke), and maintains all the
**  original elements of Anderson's code.  This has not otherwise been altered
**  and generated variables may clash with those in the user's memory.

version 10
cap program drop andersonFDR
program andersonFDR
syntax varname

gen pval = `1'

quietly sum pval
local totalpvals = r(N)

* Sort the p-values in ascending order and generate a variable that codes each p-value's rank

quietly gen int original_sorting_order = _n
quietly sort pval
quietly gen int rank = _n if pval~=.

* Set the initial counter to 1 

local qval = 1

* Generate the variable that will contain the BKY (2006) sharpened q-values

gen bky06_qval = 1 if pval~=.

* Set up a loop that begins by checking which hypotheses are rejected at q = 1.000, then checks which hypotheses are rejected at q = 0.999, then checks which hypotheses are rejected at q = 0.998, etc.  The loop ends by checking which hypotheses are rejected at q = 0.001.


while `qval' > 0 {
	* First Stage
	* Generate the adjusted first stage q level we are testing: q' = q/1+q
	local qval_adj = `qval'/(1+`qval')
	* Generate value q'*r/M
	gen fdr_temp1 = `qval_adj'*rank/`totalpvals'
	* Generate binary variable checking condition p(r) <= q'*r/M
	gen reject_temp1 = (fdr_temp1>=pval) if pval~=.
	* Generate variable containing p-value ranks for all p-values that meet above condition
	gen reject_rank1 = reject_temp1*rank
	* Record the rank of the largest p-value that meets above condition
	egen total_rejected1 = max(reject_rank1)

	* Second Stage
	* Generate the second stage q level that accounts for hypotheses rejected in first stage: q_2st = q'*(M/m0)
	local qval_2st = `qval_adj'*(`totalpvals'/(`totalpvals'-total_rejected1[1]))
	* Generate value q_2st*r/M
	gen fdr_temp2 = `qval_2st'*rank/`totalpvals'
	* Generate binary variable checking condition p(r) <= q_2st*r/M
	gen reject_temp2 = (fdr_temp2>=pval) if pval~=.
	* Generate variable containing p-value ranks for all p-values that meet above condition
	gen reject_rank2 = reject_temp2*rank
	* Record the rank of the largest p-value that meets above condition
	egen total_rejected2 = max(reject_rank2)

	* A p-value has been rejected at level q if its rank is less than or equal to the rank of the max p-value that meets the above condition
	replace bky06_qval = `qval' if rank <= total_rejected2 & rank~=.
	* Reduce q by 0.001 and repeat loop
	drop fdr_temp* reject_temp* reject_rank* total_rejected*
	local qval = `qval' - .001
}
	

quietly sort original_sorting_order
pause off
set more on

display "Code has completed."
display "Benjamini Krieger Yekutieli (2006) sharpened q-vals are in variable 'bky06_qval'"
display	"Sorting order is the same as the original vector of p-values"


* Note: Sharpened FDR q-vals can be LESS than unadjusted p-vals when many hypotheses are rejected, because if you have many true rejections, then you can tolerate several false rejections too (this effectively just happens for p-vals that are so large that you are not going to reject them regardless).
end
