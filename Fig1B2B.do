clear
capture log close
cd "YOUR_PATH_HERE"

log using "Fig1B2B.log", replace

/* X missing dependent on Y, Y MCAR */
local sampsize = 1000

set seed 123456789



set obs `sampsize'
gen id=_n



gen truecoefyx=.
gen ccacoefyx=. 

gen micoefyxce=.

gen micoefyxcesuby=.
gen micoefyxcesubx=.


gen setruecoefyx=.
gen seccacoefyx=. 
gen semicoefyxce=.

gen semicoefyxcesuby=.
gen semicoefyxcesubx=.



gen pmxmy=.
gen pxmy=.
gen pmxy=.
gen pxy=.


local nsims=500
local xcoef=0.15
local b=sqrt(0.5^2-(`xcoef'^2*0.5^2))

local sim=1
foreach sim of numlist 1/`nsims' {
	di `sim' 
	quietly{ 
	
gen x=invnorm(uniform())*0.5+1
gen y=x*`xcoef'+invnorm(uniform())*`b'+0.85



/* Scenario 1B, 2B*/
centile y, centile(50)
gen try1=r(c_1)

gen pselx=0.9 if y<=try1
replace pselx=0.1 if y>try1


gen psely=0.5

gen sely=psely>=uniform()
gen selx=pselx>=uniform()
gen sel=selx==1&sely==1


gen cons=1
summ cons if sely==0&selx==0
replace pmxmy=r(N) if id==`sim'
summ cons if sely==0&selx==1
replace pxmy=r(N) if id==`sim'

summ cons if sely==1&selx==0
replace pmxy=r(N) if id==`sim'

summ cons if sely==1&selx==1
replace pxy=r(N) if id==`sim'


reg y x
replace truecoefyx=_coef[x] if id==`sim' 
replace setruecoefyx=_se[x] if id==`sim' 
reg y x if sel==1
replace ccacoefyx=_coef[x] if id==`sim' 
replace seccacoefyx=_se[x] if id==`sim' 



replace y=. if sely~=1
replace x=. if selx~=1



preserve


mi set wide
mi register imputed x y
mi impute chained (regress)  y x, add(20) 

mi estimate: regress y x
matrix define C=r(table)


restore
preserve 

keep if x~=.
mi set wide
mi register imputed y
mi impute regress  y x, add(20) 


mi estimate: regress y x
matrix define M=r(table)

restore
preserve 

keep if y~=.
mi set wide
mi register imputed x
mi impute regress  x y, add(20) 


mi estimate: regress y x
matrix define K=r(table)


restore

gen trythis3=C[1,1]
replace micoefyxce=trythis3 if id==`sim'

gen trythis9=K[1,1]
replace micoefyxcesuby=trythis9 if id==`sim'

gen trythis10=M[1,1]
replace micoefyxcesubx=trythis10 if id==`sim'
drop trythis*

gen trythis3=C[2,1]
replace semicoefyxce=trythis3 if id==`sim'


gen trythis9=K[2,1]
replace semicoefyxcesuby=trythis9 if id==`sim'
drop trythis*


gen trythis10=M[2,1]
replace semicoefyxcesubx=trythis10 if id==`sim'
drop trythis*


matrix drop C  K M


keep id true* cca* micoef* setrue* secca* semicoef* pmxmy pxmy pmxy pxy
		

	}
		local sim = `sim' + 1
		}
		
keep if micoefyxce~=.		
summ

gen biasyxce=micoefyxce-`xcoef'
gen biasyxcca=ccacoefyx-`xcoef'
gen biasyxtrue=truecoefyx-`xcoef'
gen biasyxcesuby=micoefyxcesuby-`xcoef'
gen biasyxcesubx=micoefyxcesubx-`xcoef'


gen inyxtrue=((truecoefyx-1.96*setruecoefyx)<`xcoef')&(`xcoef'<(truecoefyx+1.96*setruecoefyx))
gen inyxcca=((ccacoefyx-1.96*seccacoefyx)<`xcoef')&(`xcoef'<(ccacoefyx+1.96*seccacoefyx))
gen inyxce=((micoefyxce-1.96*semicoefyxce)<`xcoef')&(`xcoef'<(micoefyxce+1.96*semicoefyxce))
gen inyxcesuby=(micoefyxcesuby-1.96*semicoefyxcesuby)<`xcoef'&`xcoef'<(micoefyxcesuby+1.96*semicoefyxcesuby)
gen inyxcesubx=(micoefyxcesubx-1.96*semicoefyxcesubx)<`xcoef'&`xcoef'<(micoefyxcesubx+1.96*semicoefyxcesubx)

summ bias*
summ in*
summ p*

ttest biasyxce==0
ttest biasyxcca==0
ttest biasyxtrue==0
ttest biasyxcesuby==0
ttest biasyxcesubx==0


exit

