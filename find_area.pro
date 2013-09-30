pro  find_area

nlon = 640
nlat = 400
;China,limit=[15,72,55,136]
no2 = fltarr(nlon,nlat)
lon = fltarr(nlon)
lat = fltarr(nlat)
lon = 72+indgen(nlon)*0.1
lat = 55-indgen(nlat)*0.1

;ID is used for mark the study area
;ID=4 out of Bohai; ID =1 land; ID=2 coastline; ID=3 near coast; ID=0 study area
ID = fltarr(nlon,nlat)
;size is used for defining how many km to be exclude from the coast
size = 50

header = strarr(6,1)
filename = '/home/liufei/Data/Bay/bohai_land.asc'
openr,lun,filename,/get_lun
readf,lun,header,no2
;exclued the small island [120.9,38.4]
no2[488,166]=-9999

Loc_lon=[-9999]
Loc_lat=[-9999]
For I=0, nlon-1 do begin
	For J=0,nlat-1 do begin
		if no2[I,J] gt 0 then begin
			ID[I,J]=1
			if (no2[I-1,J] eq -9999) or (no2[I+1,J] eq -9999) $
			or (no2[I,J-1] eq -9999) or (no2[I,J+1] eq -9999) then begin
				Loc_lon=[Loc_lon,I]
				Loc_lat=[Loc_lat,J]
				ID[I,J]=2
			endif
		endif
	endfor
endfor

print,'num of coastline',n_elements(Loc_lat)

For I=0, nlon-1 do begin
        For J=0,nlat-1 do begin
		if no2[I,J] eq -9999 then begin
			For k=1, n_elements(Loc_lat)-1 do begin
			distance= MAP_2POINTS(lon[I],lat[J],lon[Loc_lon[k]],lat[Loc_lat[k]],/meters)/1000
			if distance le size then begin
				ID[I,J]=3
			endif
			;****************Bohai Boundary***********************
			if (lon[I] gt 121.1) or (lon[I] lt 117.35) or (lat[J] gt 41) or (lat[J] lt 37.15) then begin
                                ID[I,J]=4
                        endif
		endfor
		endif
	endfor
endfor

print,'num of area for study', n_elements(where(ID lt 1))

location=[[-9999,-9999]]
For I=0, nlon-1 do begin
        For J=0,nlat-1 do begin
		if ID[I,J] lt 1 then begin
			location=[location,[lat[J]-0.05,lon[I]+0.05]]
		endif
	endfor
endfor

location2= reform(location,2,n_elements(location)/2)	
header_output = [['ncols 640'],['nrows 400'],['xllcorner 72'],['yllcorner 15'],['cellsize 0.1'],['nodata_value -999.0']]
outfile = '/home/liufei/Data/Bay/mask_area_0.1.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,ID
close,lun			

outfile = '/home/liufei/Data/Bay/loc_area_0.1.asc'
openw,lun,outfile,/get_lun
printf,lun,location2
close,lun
Free_LUN,lun

end
