;define x = t/12, F(x)-C = B*x+k11*sin(2Pi*x)+k12*sin(4Pi*x)+k13*sin(6Pi*x)+k14*sin(8Pi*x)+k21*cos(2Pi*x)+k22*cos(4Pi*x)+k23*cos(6Pi*x)+k24*cos(8Pi*x)
;define a procedure to return F(x) and the partial derivatives,given x. Note that A is an array containing the values B,k11,k12,k13,k14,k21,k22,k23,k24,C
;the partial derivatives are computed:
;dF/dB = x
;dF/dk11 = sin(!DPI*2*x)
;dF/dk12 = sin(!DPI*4*x)
;dF/dk13 = sin(!DPI*6*x)
;dF/dk14 = sin(!DPI*8*x)
;dF/dk21 = cos(!DPI*2*x)
;dF/dk22 = cos(!DPI*4*x)
;dF/dk23 = cos(!DPI*6*x)
;dF/dk24 = cos(!DPI*8*x)
;dF/dC = 1
pro gfunct,x,A,F,pder
        k11x = A[1]*sin(!DPI/6*x)
        k12x = A[2]*sin(!DPI/3*x)
        k13x = A[3]*sin(!DPI/2*x)
        k14x = A[4]*sin(!DPI*(2.0/3.0)*x)
        k21x = A[5]*cos(!DPI/6*x)
        k22x = A[6]*cos(!DPI/3*x)
        k23x = A[7]*cos(!DPI/2*x)
        k24x = A[8]*cos(!DPI*(2.0/3.0)*x)
        F =(1.0/12.0)* A[0]*x+k11x+k12x+k13x+k14x+k21x+k22x+k23x+k24x+A[9]
; If the procedure is called with four parameters, calculate the partial derivatives.
        if N_PARAMS() ge 4 then $
                pder = [[(1.0/12.0)*x],[sin(!DPI/6*x)],[sin(!DPI/3*x)],[ sin(!DPI/2*x)],[sin(!DPI*(2.0/3.0)*x)],[cos(!DPI/6*x)],[cos(!DPI/3*x)],[cos(!DPI/2*x)],[cos(!DPI*(2.0/3.0)*x)],[replicate(1.0,N_ELEMENTS(x))]]
end



pro trend_fit_omi_125x125
;input data
nlon_g=2880
nlat_g=1440
;global density
global = fltarr(nlon_g,nlat_g)

;china density,limit=[15,72,55,136]
;nlon = 256*2
;nlat = 160*2
;bay density,limit=[37,117,41,121]
start_lon=117.25
end_lon=123
start_lat=37
end_lat=41
nlon =(end_lon-start_lon)/0.125
nlat =(end_lat-start_lat)/0.125
no2 = fltarr(nlon,nlat)
average = fltarr(nlon,nlat)
number = fltarr(nlon,nlat)
;star point of grid for China in global map
;slon = (72+180)/0.125
;slat = (90-55)/0.125
;star point of grid for bay in global map
slon = (start_lon+180)/0.125
slat = (90-end_lat)/0.125

flag = 0U
flag2 = 0U
For year = 2004, 2011 do begin
For month = 1,12 do begin

	Yr4  = string(year,format='(i4.4)')
	Mon2 = string(month,format='(i2.2)')
	nymd = year * 10000L + month * 100L + 1 * 1L
        
	if nymd eq 20040101 then continue
        if nymd eq 20040201 then continue
        if nymd eq 20040301 then continue
        if nymd eq 20040401 then continue
        if nymd eq 20040501 then continue
        if nymd eq 20040601 then continue
        if nymd eq 20040701 then continue
        if nymd eq 20040801 then continue
        if nymd eq 20040901 then continue


	header = strarr(7,1)
	filename='/z6/satellite/OMI/no2/KNMI_L3_v2/no2_'+Yr4+Mon2+'.grd'
	openr,lun,filename,/get_lun
	readf,lun,header,global
	no2 = global[slon:slon+nlon-1,slat:slat+nlat-1]/100

	;calculate long-term avearage
	For J=0, nlat-1 do begin
	For I = 0,nlon-1 do begin
		if no2[I,J] ge 0.0 then begin
			average[I,J]+= no2[I,J]
			number[I,J] +=1
		endif
	endfor
	endfor


	if flag  then begin
	;m is used for counting the total number of months
	m = m+1
	;define no2_month to save no2 for all months
	no2_month =[ [no2_month],[no2] ]
	endif else begin
	m = 1U
	no2_month = no2
	flag = 1U
	endelse
	;print,Yr4,Mon2
	close,/all

endfor
endfor
print,'m',m
print,'size of no2_month',size(no2_month)

For J=0, nlat-1 do begin
        For I = 0,nlon-1 do begin
                if number[I,J] gt 0.0 then begin
                        average[I,J] = average[I,J]/number[I,J]   
;filter the data 9-year average lt 1*10^15 molecules/cm2
			if average[I,J] lt 1.0 then begin
				average[I,J] = -999.0
			endif
                endif else begin
			average[I,J] = -999.0
		endelse
        endfor
endfor


;convert no2_month to 3-D array no2_data
no2_data =fltarr (nlon,nlat,m)
For num = 0,m-1 do begin
	no2_data[*,*,num] = no2_month[0:nlon-1,nlat*num:nlat*(num+1)-1]
endfor
;find no2_data with reasonable value( filter value < 0)
loc = where(no2_data lt 0)
;print,'loc',loc
if not array_equal(loc,[-1]) then begin
	no2_data[loc] = -999.0
endif
undefine,loc
print,'size of no2_data',size(no2_data)
print,'MAX OF no2_data',MAX(no2_data[*,*,*]),'MIN OF no2_data',MIN(no2_data[*,*,*])

;Compute the fit to the function we have just defined. First, define the independent and dependent variables
;define para to save fit parameters
para = fltarr(nlon,nlat)
;deltab=deltan/n^(3/2)*((1+cor)/(1-cor))^1/2
;deltan is the variance in the remainder
cor=fltarr(nlon,nlat)
deltan=fltarr(nlon,nlat)
deltab=fltarr(nlon,nlat)

nodata_delta_min = 0
nodata_para=0
For J = 0,nlat-1 do begin
FOR I = 0,nlon-1 do begin
    flag2 = 0U
    nodata_no2 =0
    For num = 0,m-1 do begin
	 if no2_data[I,J,num] gt -999.0 then begin
		if flag2 then begin
			Y = [Y,no2_data[I,J,num] ]
			x = [x,num+8]
		endif else begin
			Y = no2_data[I,J,num] 
			x = [num+8]
			flag2 = 1U
		endelse
          endif else begin
		nodata_no2 = nodata_no2+1
 	  endelse
     endfor
	
   if nodata_no2 eq m then begin
		para[I,J]=-999.0
		nodata_para=nodata_para+1
   endif else begin	

	if N_elements(x) le 10 then begin
		para[I,J]=-999.0
		deltab[I,J]=1.0
		print,'Elements <10 -999.0',I,J
	endif else begin
		weights = make_array(N_elements(x),value=1.0) 
		;Provide an initial guess of the function's parameters.
		A = [0.02,2.666,4.475,2.725,3.785,2.003,4.498,6.22,5.0345,5]
		yfit = CURVEFIT(x,y,weights,A,SIGMA,FUNCTION_NAME='gfunct',ITMAX=100)
		para[I,J]= A[0]
		;if para[I,J] le 0.0 then begin
			;print,'I',I,'J',J,'PARA <0',para[I,J]
		;endif

		k11x = A[1]*sin(!DPI/6*x)
	        k12x = A[2]*sin(!DPI/3*x)
	       	k13x = A[3]*sin(!DPI/2*x)
	       	k14x = A[4]*sin(!DPI*(2.0/3.0)*x)
	        k21x = A[5]*cos(!DPI/6*x)
	        k22x = A[6]*cos(!DPI/3*x)
	        k23x = A[7]*cos(!DPI/2*x)
	        k24x = A[8]*cos(!DPI*(2.0/3.0)*x)
	        F = A[0]*(1.0/12.0)*x+k11x+k12x+k13x+k14x+k21x+k22x+k23x+k24x+A[9]
		remain = F-Y
		;print,'remain',remain
		lag_cor = [1]
		cor = A_CORRELATE(remain,lag_cor,/double)
		;print,'cor',cor
		deltan[I,J]=STDDEV(remain)
		;n is total year
		;n=9
		n=N_elements(x)/12
		deltab[I,J] = deltan[I,J]/(n^(1.5))*(((1+cor[0])/(1-cor[0]))^(0.5))
		;filter deltab lt delta_min
;		delta_min = 0.65 - A[9] + 0.3*mean(y)
;		if  deltab[I,J] lt delta_min then begin
;			deltab[I,J] = delta_min
;			nodata_delta_min +=1
;		endif
	
	endelse

    endelse

endfor
endfor

print,'nodata_para',nodata_para
print,'nodata_delta_min',nodata_delta_min
;filter the data 9-year average lt 1*10^15 molecules/cm2
loc = where(average  eq -999.0)
print,'number of background',n_elements(loc)
if not array_equal(loc,[-1]) then begin
para[loc]=-999.0
endif
undefine,loc
;tw is the value of  the student's t-distribution for a significance level of 0.05 and the degrees of free given for the time series
;select abs(B/deltab)>tw
tw = T_CVF(0.05,N_elements(x)-1)
para2=para
loc = where(abs(para/deltab) le tw)
if not array_equal(loc,[-1]) then begin
para2[loc]= -999.0
print,'nodata_para2',N_elements(loc)
endif
undefine,loc

loc = where( (para2 gt -999.0) and (average gt 0.0) )
para3 = fltarr(nlon,nlat)
For J=0, nlat-1 do begin
para3[*,J] = replicate(-999.0,nlon)
endfor
para3[loc]= para2[loc]/average[loc]
undefine,loc

loc = where( (para gt -999.0) and (average gt 0.0) )
para4 = fltarr(nlon,nlat)
For J=0, nlat-1 do begin
para4[*,J] = replicate(-999.0,nlon)
endfor
para4[loc]= para[loc]/average[loc]


;header_output = [['ncols 512'],['nrows 320'],['xllcorner 72'],['yllcorner 15'],['cellsize 0.125'],['nodata_value -999.0']]
header_output = [['ncols '+string(nlon)],['nrows '+string(nlat)],['xllcorner '+string(start_lon)],['yllcorner '+string(start_lat)],['cellsize 0.125'],['nodata_value -999.0']]
outfile = '/home/liufei/Data/Satellite/NO2/trend/trend_fit_omi_0.125X0.125.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,para
close,lun

outfile = '/home/liufei/Data/Satellite/NO2/trend/trend_fit2_omi_0.125X0.125.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,para2
close,lun

outfile = '/home/liufei/Data/Satellite/NO2/trend/trend_fit3_omi_0.125X0.125.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,para3
close,lun

outfile = '/home/liufei/Data/Satellite/NO2/trend/trend_fit4_omi_0.125X0.125.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,para4
close,lun
end
