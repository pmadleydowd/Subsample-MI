clear
capture log close
cd "YOUR_PATH_HERE"

log using "Fig4B.log", replace

* Y missing dependent on X and C, C missing dependent on X and U, X complete

local sampsize = 1000

set seed 123456789



set obs `sampsize'
gen id=_n



gen truecoefyx=.
gen ccacoefyx=. 

gen micoefyxce=.

gen micoefyxcesubc=.
gen micoefyxcesuby=.



gen setruecoefyx=.
gen seccacoefyx=. 
gen semicoefyxce=.

gen semicoefyxcesubc=.
gen semicoefyxcesuby=.





gen pmxmymc=.
gen pxmymc=.
gen pmxymc=.
gen pxymc=.
gen pmxmyc=.
gen pxmyc=.
gen pmxyc=.
gen pxyc=.


local nsims=500
local xcoef=0.15
local yerr1=0.15^2/4+(0.5^2)/4-2*0.15*0.5*(sqrt(0.5))/4

local yerr2=(0.25-`yerr1')
di `yerr2'

local yerr=sqrt(`yerr2'*2)
di `yerr1'
di `yerr2'
di `yerr'

local sim=1
foreach sim of numlist 1/`nsims' {
	di `sim' 
	quietly{ 
gen u=invnorm(uniform())*0.5+1
gen c=invnorm(uniform())*0.5+1
gen x=sqrt(0.5)*c+sqrt(0.5)*invnorm(uniform())*0.5+(1-sqrt(0.5))
gen y=x*`xcoef'+invnorm(uniform())*0.5*`yerr'+(1-0.15+1+0.5-0.35)-0.5*c-u*`yerr'



/* Scenario RY dep on C and X, Rc dep on U and X*/
centile c, centile(50)
gen try1=r(c_1)
centile x, centile(50)
gen try2=r(c_1)
centile u, centile(50)
gen try3=r(c_1)

gen pselc=0.5
replace pselc=0.9 if x<=try2&u<=try3
replace pselc=0.1 if x>try2&u>try3

gen psely=0.5
replace psely=0.9 if x<=try2&c<=try1
replace psely=0.1 if x>try2&c>try1


gen selx=1

gen selc=pselc>=uniform()
gen sely=psely>=uniform()
gen sel=selx==1&sely==1&selc==1


gen cons=1
summ cons if sely==0&selx==0&selc==0
replace pmxmymc=r(N) if id==`sim'
summ cons if sely==0&selx==1&selc==0
replace pxmymc=r(N) if id==`sim'

summ cons if sely==1&selx==0&selc==0
replace pmxymc=r(N) if id==`sim'

summ cons if sely==1&selx==1&selc==0
replace pxymc=r(N) if id==`sim'
summ cons if sely==0&selx==0&selc==1
replace pmxmyc=r(N) if id==`sim'
summ cons if sely==0&selx==1&selc==1
replace pxmyc=r(N) if id==`sim'

summ cons if sely==1&selx==0&selc==1
replace pmxyc=r(N) if id==`sim'

summ cons if sely==1&selx==1&selc==1
replace pxyc=r(N) if id==`sim'

reg y x c
replace truecoefyx=_coef[x] if id==`sim' 
replace setruecoefyx=_se[x] if id==`sim' 
reg y x c if sel==1
replace ccacoefyx=_coef[x] if id==`sim' 
replace seccacoefyx=_se[x] if id==`sim' 

replace y=. if sely~=1
replace x=. if selx~=1
replace c=. if selc~=1


preserve

mi set wide
mi register imputed x y c
mi impute chained (regress)  y x c, add(20) 

mi estimate: regress y x c
matrix define C=r(table)

restore
preserve 

keep if y~=.
mi set wide
mi register imputed  c
mi impute chained (regress)   c = x y, add(20) 


mi estimate: regress y x c
matrix define K=r(table)

restore
preserve 

keep if c~=.
mi set wide
mi register imputed  y 
mi impute chained (regress)   y =  x c, add(20) 


mi estimate: regress y x c
matrix define L=r(table)



restore


gen trythis3=C[1,1]
replace micoefyxce=trythis3 if id==`sim'

gen trythis9=K[1,1]
replace micoefyxcesuby=trythis9 if id==`sim'

gen trythis10=L[1,1]
replace micoefyxcesubc=trythis10 if id==`sim'



drop trythis*

gen trythis3=C[2,1]
replace semicoefyxce=trythis3 if id==`sim'


gen trythis9=K[2,1]
replace semicoefyxcesuby=trythis9 if id==`sim'

gen trythis10=L[2,1]
replace semicoefyxcesubc=trythis10 if id==`sim'


drop trythis*


matrix drop C  K L 


keep id true* cca* micoef* setrue* secca* semicoef* pmxmy* pxmy* pmxy* pxy*
		

	}
		local sim = `sim' + 1
		}
		
keep if micoefyxce~=.		
summ

gen biasyxce=micoefyxce-`xcoef'
gen biasyxcca=ccacoefyx-`xcoef'
gen biasyxtrue=truecoefyx-`xcoef'
gen biasyxcesubc=micoefyxcesubc-`xcoef'
gen biasyxcesuby=micoefyxcesuby-`xcoef'



gen inyxtrue=((truecoefyx-1.96*setruecoefyx)<`xcoef')&(`xcoef'<(truecoefyx+1.96*setruecoefyx))
gen inyxcca=((ccacoefyx-1.96*seccacoefyx)<`xcoef')&(`xcoef'<(ccacoefyx+1.96*seccacoefyx))
gen inyxce=((micoefyxce-1.96*semicoefyxce)<`xcoef')&(`xcoef'<(micoefyxce+1.96*semicoefyxce))
gen inyxcesubc=(micoefyxcesubc-1.96*semicoefyxcesubc)<`xcoef'&`xcoef'<(micoefyxcesubc+1.96*semicoefyxcesubc)
gen inyxcesuby=(micoefyxcesuby-1.96*semicoefyxcesuby)<`xcoef'&`xcoef'<(micoefyxcesuby+1.96*semicoefyxcesuby)



summ bias*
summ in*
summ p*

ttest biasyxce==0
ttest biasyxcca==0
ttest biasyxtrue==0
ttest biasyxcesuby==0
ttest biasyxcesubc==0



exit

