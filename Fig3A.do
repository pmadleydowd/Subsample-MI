clear
capture log close
cd "YOUR_PATH_HERE"

log using "Fig3A.log", replace

* X missing dependent on C, C missing dependent on X, Y MCAR

local sampsize = 1000

set seed 12345678



set obs `sampsize'
gen id=_n



gen truecoefyx=.
gen ccacoefyx=. 

gen micoefyxce=.

gen micoefyxcesubc=.
gen micoefyxcesubx=.
gen micoefyxcesubxc=.


gen setruecoefyx=.
gen seccacoefyx=. 
gen semicoefyxce=.

gen semicoefyxcesubc=.
gen semicoefyxcesubx=.
gen semicoefyxcesubxc=.




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


local yerr=sqrt(`yerr2')




local sim=1
foreach sim of numlist 1/`nsims' {
	di `sim' 
	quietly{ 
gen c=invnorm(uniform())*0.5+1
gen x=sqrt(0.5)*c+sqrt(0.5)*invnorm(uniform())*0.5+(1-sqrt(0.5))

gen y=x*`xcoef'+invnorm(uniform())*`yerr'+(1-0.15+0.5)-0.5*c



/* Scenario Rx dep on C, Rc dep on X*/
centile c, centile(50)
gen try1=r(c_1)

gen pselx=0.9 if c<=try1
replace pselx=0.1 if c>try1

centile x, centile(50)
gen try2=r(c_1)

gen pselc=0.9 if x<=try2
replace pselc=0.1 if x>try2

gen sely=uniform()>=0.5

gen selc=pselc>=uniform()
gen selx=pselx>=uniform()
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

keep if x~=.
mi set wide
mi register imputed y c
mi impute chained (regress)  y c =  x, add(20) 


mi estimate: regress y x c
matrix define K=r(table)

restore
preserve 

keep if c~=.
mi set wide
mi register imputed  y x
mi impute chained (regress)   y x=  c, add(20) 


mi estimate: regress y x c
matrix define L=r(table)

restore
preserve 

keep if x~=.&c~=.
mi set wide
mi register imputed y
mi impute regress  y=c x, add(20) 


mi estimate: regress y x c
matrix define M=r(table)

restore


gen trythis3=C[1,1]
replace micoefyxce=trythis3 if id==`sim'

gen trythis9=K[1,1]
replace micoefyxcesubx=trythis9 if id==`sim'

gen trythis10=L[1,1]
replace micoefyxcesubc=trythis10 if id==`sim'

gen trythis11=M[1,1]
replace micoefyxcesubxc=trythis11 if id==`sim'

drop trythis*

gen trythis3=C[2,1]
replace semicoefyxce=trythis3 if id==`sim'


gen trythis9=K[2,1]
replace semicoefyxcesubx=trythis9 if id==`sim'

gen trythis10=L[2,1]
replace semicoefyxcesubc=trythis10 if id==`sim'

gen trythis11=M[2,1]
replace semicoefyxcesubxc=trythis11 if id==`sim'
drop trythis*



matrix drop C  K L M


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
gen biasyxcesubx=micoefyxcesubx-`xcoef'
gen biasyxcesubxc=micoefyxcesubxc-`xcoef'


gen inyxtrue=((truecoefyx-1.96*setruecoefyx)<`xcoef')&(`xcoef'<(truecoefyx+1.96*setruecoefyx))
gen inyxcca=((ccacoefyx-1.96*seccacoefyx)<`xcoef')&(`xcoef'<(ccacoefyx+1.96*seccacoefyx))
gen inyxce=((micoefyxce-1.96*semicoefyxce)<`xcoef')&(`xcoef'<(micoefyxce+1.96*semicoefyxce))
gen inyxcesubc=(micoefyxcesubc-1.96*semicoefyxcesubc)<`xcoef'&`xcoef'<(micoefyxcesubc+1.96*semicoefyxcesubc)
gen inyxcesubx=(micoefyxcesubx-1.96*semicoefyxcesubx)<`xcoef'&`xcoef'<(micoefyxcesubx+1.96*semicoefyxcesubx)
gen inyxcesubxc=(micoefyxcesubxc-1.96*semicoefyxcesubxc)<`xcoef'&`xcoef'<(micoefyxcesubxc+1.96*semicoefyxcesubxc)


summ bias*
summ in*
summ p*

ttest biasyxce==0
ttest biasyxcca==0
ttest biasyxtrue==0
ttest biasyxcesubc==0
ttest biasyxcesubx==0
ttest biasyxcesubxc==0


exit

