pro  direction_column

file_in_lat=string('/home/liufei/Data/Large_point/Steffen/latitude.nc')
fid=NCDF_OPEN(file_in_lat)
varid=NCDF_VARID(fid,'map_TVCD')
NCDF_VARGET,fid,varid,lat
NCDF_CLOSE, fid

file_in_lon=string('/home/liufei/Data/Large_point/Steffen/longitude.nc')
fid=NCDF_OPEN(file_in_lon)
varid=NCDF_VARID(fid,'map_TVCD')
NCDF_VARGET,fid,varid,lon
NCDF_CLOSE, fid

file_in=string('/home/liufei/Data/Large_point/Steffen/map_TVCD.nc')
fid=NCDF_OPEN(file_in)
varid=NCDF_VARID(fid,'map_TVCD')
NCDF_VARGET,fid,varid,DATA
NCDF_CLOSE, fid
;DATA(lat,lon,sea,winddirection,windspeed)
close,/all
;bay density,limit=[37,117.25,41,123]
start_lon=117.25
end_lon=123
start_lat=37
end_lat=41
lon_array=where((lon le end_lon) and (lon ge start_lon) )
lat_array=where((lat le end_lat) and (lat ge start_lat) )
nlon =N_elements(lon_array)
nlat =N_elements(lat_array)
;clip no2 column
no2=DATA[lat_array,lon_array,*,*,*]
print,'nlon',nlon,'nlat',nlat

;mask the sea boundary
nlon_mask=(end_lon-start_lon)/0.125
nlat_mask=(end_lat-start_lat)/0.125
print,'nlon_mask',nlon_mask,'nlat_mask',nlat_mask
mask=fltarr (nlon_mask,nlat_mask)
header = strarr(6,1)
filename='/home/liufei/r5/bay/sea_boundary.asc'
openr,lun,filename,/get_lun
readf,lun,header,mask
;location of mask
;mask_I:start_lon+mask_I*0.125,mask_J:end_lat-mask_J*0.125

;build an array to save paramaters corresponding to the location after mask
para=fltarr(nlon,nlat)
For J= 0 , nlat-1 do begin
	For I= 0 , nlon-1 do begin
		For mask_J=0,nlat_mask-1 do begin
		    if (lat[ lat_array[J] ] ge end_lat-mask_J*0.125-0.125/2) and $
		       (lat[ lat_array[J] ] le end_lat-mask_J*0.125+0.125/2) then begin
			For mask_I=0,nlon_mask-1 do begin
				if (lon[ lon_array[I] ] ge start_lon+mask_I*0.125-0.125/2) and $
				   (lon[ lon_array[I] ] le start_lon+mask_I*0.125+0.125/2) then begin
					para[I,J]=mask[mask_I,mask_J]
				endif
			endfor
		    endif
		endfor
	endfor
endfor
					
header_output = [['ncols '+string(nlon)],['nrows '+string(nlat)],['xllcorner '+string(start_lon)],$
		['yllcorner '+string(start_lat)],['cellsize 0.176'],['nodata_value -9999.0']]
outfile = '/home/liufei/Data/Satellite/NO2/trend/para.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,para
close,lun

average=fltarr(9,1)
num=fltarr(9,1)
For dir = 0,8 do begin
        TVCD=fltarr(nlon,nlat)
        For J= 0 , nlat-1 do begin
                For I= 0 , nlon-1 do begin
			;For season=0,3 do begin
			season=2
                        if (finite(no2[J,I,season,dir,2]) eq 0) or (no2[J,I,season,dir,2] lt 0) then begin
                                TVCD[I,J]=-999
                        endif else begin
                                TVCD[I,J]= no2[J,I,season,dir,2]
				if para[I,J] eq -9999 then begin
					average[dir]+=no2[J,I,season,dir,2]
					num[dir]+=1
				endif
                        endelse
			;endfor
		endfor
	endfor
	header_output = [['ncols '+string(nlon)],['nrows '+string(nlat)],['xllcorner '+string(start_lon)],['yllcorner '+string(start_lat)],['cellsize 0.176'],['nodata_value -999.0']]
	outfile = '/home/liufei/Data/Satellite/NO2/trend/VCD_dir'+string(dir)+'.asc'
	openw,lun,outfile,/get_lun
	printf,lun,header_output,TVCD
	close,lun
endfor
print,'average',average/num

loadct,5
Y=average/num
bar_plot,y,background=!d.n_colors-1,barnames=indgen(9)+1,barwidth=0.6,barspace=0.3,colors=replicate(0.25*!d.n_colors,9),$
	baselines=replicate(3,9),baserange=1,$
	title='average NO2 column by wind direction',xtitle='direction',$
	ytitle='NO2/(10^15 molec/cm2)'
image = tvrd(true =1)
write_jpeg,'average column.jpg',image,true=1
end
