 /*_______________________________________________________________________
	|                                                                      |
	| PROJECT: Event Study Template (PSM + Dynamic DiD)                    |
	| Author: Jesus Gutierrez Andrade                                      |
	| Version: 25/04/2026                                                  |
	|______________________________________________________________________|	*/
	
	clear all
	
	set more off
	set varabbrev off
	set type double
	set seed 94783

   /*_______________________________________________________________________
	|                                                                      |
	|                  	    Set working directory                  		   |
	|_____________________________________________________________________*/

*	cd "your/working/directory"


/*	_______________________________________________________________________
	|                                                                      |
	|     		      EVENT STUDY ESTIMATION (one outcome)		           |
	|______________________________________________________________________|	*/

*	Required panel structure:
*		firm	   : panel identifier (cross-section unit)
*		year	   : time variable
*		outcome	   : variable of interest (e.g. sales, profits, etc.)
*		treated	   : 1 if unit is in the treatment group, 0 otherwise
*		period     : 2017 - 2022
*		event year : t = 2020 (replace with your own event date)

	u temp/your_panel_data, clear
		
*	Outcome label (used for graph titles and exports):

	gl yt_outcome	"Outcome (units)"
		
*	Log-transform outcome:
	
	cap drop ln_outcome
		 gen ln_outcome = ln(outcome)			   // In case of huge quantities
	
	sort firm year
	
*	Pre-event mean of the outcome (matching variable):
	
	cap drop matvar
	by firm: egen matvar = mean(ln_outcome) if inrange(year,2017,2019)
	
	
   /*_______________________________________________________________________
	|                                                                      |
	|     		    1. Propensity Score Matching (pre-event)		       |
	|______________________________________________________________________|	*/
	
	psmatch2 treated matvar if year <= 2019, cal(0.01)
	
	cap drop new_T
		 gen new_T = _treated
	
	cap drop weight
		 gen weight = _weight
	
*	Balance test:
	
	pstest matvar if year <= 2019, sum rubin
	psgraph, bin(30)
	
	drop _*
	
*	Keep only matched observations from here on:
	
	cap drop _aux
	bys firm: egen _aux = max(weight)
	
	keep if _aux == 1
	
	drop _aux
	
	
/*	_______________________________________________________________________
	|                                                                      |
	|     		    2. Outcome Evolution (matched sample)		           |
	|______________________________________________________________________|	*/
	
	set scheme s2gcolor, permanently
	
	xtset firm year
	
	cap drop ave_outcome
	bys year new_T: egen ave_outcome = mean(outcome)
	
	graph two (con ave_outcome year if new_T == 1, msymbol(T) mcolor(blue)) ///
			  (con ave_outcome year if new_T == 0, msymbol(O) mcolor(red)), ///
				xline(2020) ytitle("${yt_outcome}") xtitle("Year")		   ///
				legend(order(1 "Treated" 2 "Control"))					   ///
				graphregion(color(gs16))
	
	graph export "output/evolution_${yt_outcome}.png", replace

	
   /*_______________________________________________________________________
	|                                                                      |
	|     		     3. Event-time dummies (leads & lags)		    	   |
	|______________________________________________________________________|	*/
	
	xtset firm year
	
	cap drop t_*
	tab year, gen(t_)				// t_4 corresponds to event year (2020)
	
	forval i = 1/6		{
		
		cap drop T`i'
			 gen T`i' = new_T * t_`i'
			 
	}
	
	cap drop N
		 gen N = _n
	 replace N = N - 3				// relative time (-3, -2, -1, 0, 1, 2)
	
	gl pre 		= "T1 T2"			// pre-event leads (placebo test)
	gl post	 	= "T4 T5 T6"		// post-event lags  (treatment effect)
	gl controls = "$pre $post"		// T3 = omitted base category
	
	gl Z = "year firm"				// fixed effects (if apply)
	gl y = "ln_outcome"
	
	
   /*_______________________________________________________________________
	|                                                                      |
	|     		    			 4. Estimation		    				   |
	|______________________________________________________________________|	*/
	
	reghdfe $y $controls, a(${Z})
	est store reg1
	
	mat A = r(table)'
	
*	Build a matrix with point estimates and 95% CI bounds:
	
	mat X = J(5,3,.)
	
	forval i = 1/5 {
	
		matrix X[`i', 1] = A[`i', 1]		// coefficient
		matrix X[`i', 2] = A[`i', 5]		// 95% CI lower bound
		matrix X[`i', 3] = A[`i', 6]		// 95% CI upper bound
	
	}
	
	outreg2 [reg1] using "output/r_${yt_outcome}.xls", bd(2) sd(3) replace
	
	cap drop X*
	svmat X
	
	
   /*_______________________________________________________________________
	|                                                                      |
	|     		    4.1 Event Study graph (coefficients + 95% CI)		   |
	|______________________________________________________________________|	*/
	
	twoway (rcap X2 X3 N, lcolor(gs8))									///
		   (scatter X1 N, mcolor(black) msymbol(O)),					///
				yline(0, lcolor(blue) lwidth(vthin))					///
				xline(0, lpattern(dash))								///
				ytitle("Estimated coefficient") legend(off)				///
				xtitle("Years pre and post event")						///
				graphregion(color(gs16))
	
	graph export "output/${yt_outcome}_T.png", replace
