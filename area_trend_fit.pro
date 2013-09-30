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


pro area_trend_fit
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


;convert no2_month to 3-D array no2_data
no2_data =fltarr (nlon,nlat,m)
For num = 0,m-1 do begin
	no2_data[*,*,num] = no2_month[0:nlon-1,nlat*num:nlat*(num+1)-1]
endfor
;find no2_data with reasonable value( filter value < 0)
loc = where(no2_data lt 0)
if not array_equal(loc,[-1]) then begin
	no2_data[loc] = -999.0
endif
undefine,loc
print,'size of no2_data',size(no2_data)
print,'MAX OF no2_data',MAX(no2_data[*,*,*]),'MIN OF no2_data',MIN(no2_data[*,*,*])

;mask the sea boundary
mask=fltarr (nlon,nlat)
header = strarr(6,1)
filename='/home/liufei/r5/bay/sea_boundary.asc'
openr,lun,filename,/get_lun
readf,lun,header,mask
mask_3d=fltarr (nlon,nlat,m)
For i=0,m-1 do begin
	mask_3d[*,*,i]=mask
endfor
no2_data[where(mask_3d ne -9999)]=-999.0
undefine,mask_3d

For area_case=0,4 do begin
;area_case 4 stands for the whole area
area_data=fltarr (nlon,nlat,m)
area_data=no2_data

if area_case ne 4 then begin
;mask the area boundary
area_mask=fltarr (nlon,nlat)
header2 = strarr(6,1)
filename2='/home/liufei/r5/bay/area_boundary.asc'
openr,lun,filename2,/get_lun
readf,lun,header2,area_mask
mask_3d=fltarr (nlon,nlat,m)
For i=0,m-1 do begin
	mask_3d[*,*,i]=area_mask
endfor
area_data[where(mask_3d ne area_case)]=-999.0
undefine,mask_3d
endif
;calculate long-term avearage
average=0.d
average_data=0.d
;average valid data
Y=fltarr(m,1)
valid_num=fltarr(m,1)
For num = 0,m-1 do begin
	For J = 0,nlat-1 do begin
		FOR I = 0,nlon-1 do begin
			if area_data[I,J,num] gt -999.0 then begin
				average+=area_data[I,J,num]
				average_data+=1
				Y[num]+=area_data[I,J,num]
				valid_num[num]+=1
			endif
		endfor
	endfor
endfor	
average=average/average_data
print,'average',average
Y=Y/valid_num
x=indgen(m)+8
;print,'y',y

;Compute the fit to the function we have just defined. First, define the independent and dependent variables
weights = make_array(N_elements(x),value=1.0)
;Provide an initial guess of the function's parameters.
A = [0.02,2.666,4.475,2.725,3.785,2.003,4.498,6.22,5.0345,5]
yfit = CURVEFIT(x,y,weights,A,SIGMA,FUNCTION_NAME='gfunct',ITMAX=100)
para= A[0]
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
lag_cor = [1]
cor = A_CORRELATE(remain,lag_cor,/double)
deltan=STDDEV(remain)
n=N_elements(x)/12
deltab = deltan/(n^(1.5))*(((1+cor[0])/(1-cor[0]))^(0.5))
;filter deltab lt delta_min
;               delta_min = 0.65 - A[9] + 0.3*mean(y)
;               if  deltab lt delta_min then begin
;                       deltab = delta_min
;                       nodata_delta_min +=1
;               endif

;filter the data 9-year average lt 1*10^15 molecules/cm2
;loc = where(average  eq -999.0)
;print,'number of background',n_elements(loc)

;tw is the value of  the student's t-distribution for a significance level of 0.05 and the degrees of free given for the time series
;select abs(B/deltab)>tw
;tw = T_CVF(0.05,N_elements(x)-1)
;para2=para
;loc = where(abs(para/deltab) le tw)

plot,x,y,psym=2,$
	yrange=[min(y),max(y)],$
        title='case:'+string(area_case),xtitle='number of months',ytitle='no2/(10^15moles/cm2)'
oplot,x,F,psym=0
xyouts,10,max(y)-3,$
'absolute trend(molec/(cm2*yr)):'+string(para)+'!Crelative trend(%/yr):'+string(100*para/average),$
charsize=1.5, charthick=1.8
image = tvrd(true =1)
write_jpeg,'trend_analysis_case'+string(area_case)+'.jpg',image,true=1
endfor
end

