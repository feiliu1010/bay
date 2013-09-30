pro line_density

pp_loc_lon=[120.121]
pp_loc_lat=[39.016]
file_in_lat=string('/home/liufei/Data/Large_point/Steffen/latitude.nc')
fid=NCDF_OPEN(file_in_lat)
varid=NCDF_VARID(fid,'map_TVCD')
NCDF_VARGET,fid,varid,lat
NCDF_CLOSE, fid
;reverse the latitude
lat=reverse(lat)

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
close,/all

For pp=0,N_elements(pp_loc_lon)-1 do begin
lat_point=pp_loc_lat[pp]
lon_point=pp_loc_lon[pp]

;transform TVCD coordinates from longitude and latitude to Cartesian (x, y) coordinates
xloc=fltarr(n_elements(lon),n_elements(lat))
yloc=fltarr(n_elements(lon),n_elements(lat))
For num=0,n_elements(lat)-1 do begin
	xloc[*,num]=lon
endfor
For num=0,n_elements(lon)-1 do begin
	yloc[num,*]=lat
endfor


mapStruct = MAP_PROJ_INIT( 'Lambert Conic', $
                            STANDARD_PAR1=25, STANDARD_PAR2=40, $
                            CENTER_LONGITUDE=lon_point, CENTER_LATITUDE=lat_point )
coordinates= MAP_PROJ_FORWARD( xloc, yloc, MAP_STRUCTURE = mapStruct)/1000
lon_new=reform(coordinates[0,*],[N_Elements(lon),N_Elements(lat)])
lat_new=reform(coordinates[1,*],[N_Elements(lon),N_Elements(lat)])

;map_TVCD(lat,lon,sea,winddirection,windspeed)
For dir = 0,8 do begin
	TVCD=fltarr(N_Elements(lon),N_Elements(lat))
        For J= 0 , N_Elements(lat)-1 do begin
                For I= 0 , N_Elements(lon)-1 do begin
		;************season***************
			season=2
			if (finite(DATA[J,I,season,dir,2]) eq 0) or (DATA[J,I,season,dir,2] lt 0) then begin
			;if (finite(DATA[J,I,season,dir,0]) eq 0) then begin
				TVCD[I,J]=-999
			endif else begin
				TVCD[I,J]= DATA[J,I,season,dir,2]
			endelse

		endfor
	endfor
	;reverse TVCD at latitude
	TVCD=reverse(TVCD,2)
	;rotate the coordinates of lon_new, lat_new
	case dir+1 of 
	    1:theta=135
	    2:theta=90
	    3:theta=45
            4:theta=0
	    5:theta=0
	    6:theta=0
	    7:theta=45
	    8:theta=90
	    9:theta=135
	endcase
        theta=theta/!radeg
        c=cos(theta)
        s=sin(theta)
	lon_new_rot=fltarr(N_Elements(lon),N_Elements(lat))
	lat_new_rot=fltarr(N_Elements(lon),N_Elements(lat))
        For J= 0 , N_Elements(lat)-1 do begin
                For I= 0 , N_Elements(lon)-1 do begin
			lon_new_rot[I,J]=c*lon_new[I,J]+s*lat_new[I,J]
			lat_new_rot[I,J]=c*lat_new[I,J]-s*lon_new[I,J]
		endfor
        endfor
	;x direction interval/km
	a_up=-250
	a_down=300
	;y direction interval/km
	b=150
	;find the point which is nearest to the origin point
	vector=fltarr(N_Elements(lon),N_Elements(lat))
	For J= 0 , N_Elements(lat)-1 do begin
        	For I= 0 , N_Elements(lon)-1 do begin
			if lon_new_rot[I,J] ge 0 then begin
	                	vector[I,J]= distance_measure([[0,0],[lon_new_rot[I,J],lat_new_rot[I,J]]])
			endif else begin
				vector[I,J]=-distance_measure([[0,0],[lon_new_rot[I,J],lat_new_rot[I,J]]])	
			endelse
	        endfor
	endfor
	print,'distance vector_zero',where (vector eq 0)
	print,'min of distance vector',min(abs(vector))
	nearest=where( (vector eq min(abs(vector))) or (vector eq -min(abs(vector)))  )
	nearest_new=nearest[0]
	near_row= Fix( (nearest_new+1)/N_Elements(lon) )
	near_col=( (nearest_new+1) mod N_Elements(lon))-1
	print,'near_row',near_row,'near_col',near_col
	;find boundary of interval
	bou_y=where( ((vector le b) and (vector ge 0)) or ((vector ge -b) and (vector lt 0)) )
	bou_x_up=where( (vector ge a_up) and (vector lt 0) )
	bou_x_down=where( (vector le a_down) and (vector ge 0) )
	y_temp_row= Fix( (bou_y+1)/N_Elements(lon) )
	y_temp_col=( (bou_y+1) mod N_Elements(lon))-1
        x_temp_row_up= Fix( (bou_x_up+1)/N_Elements(lon) )
        x_temp_col_up=( (bou_x_up+1) mod N_Elements(lon))-1 
	x_temp_row_down= Fix( (bou_x_down+1)/N_Elements(lon) )
        x_temp_col_down=( (bou_x_down+1) mod N_Elements(lon))-1
	;east-west
	if (dir eq 3) or (dir eq 5) then begin
		max_y1= max( y_temp_row )
		max_y2= min( y_temp_row )
		max_x1= max( x_temp_col_down )
		max_x2= min( x_temp_col_up )
		num_y= max_y1 - max_y2 +1
		num_x= max_x1 - max_x2 +1
	
		;Line Density:LD
		LD=fltarr(num_x)
		For x=max_x2 , max_x1 do begin
			For y=max_y2 , max_y1 do begin
			    if y eq 0 then begin
				dy =abs(lat_new[x,y+1]-lat_new[x,y])
			    endif else if y eq N_elements(lat)-1 then begin
				dy=abs(lat_new[x,y]-lat_new[x,y-1])
			    endif else begin
				dy=abs( (lat_new[x,y+1]-lat_new[x,y-1]) )/2
			    endelse
				if TVCD[x,y] ne -999 then begin
                                        LD[x-max_x2]+=TVCD[x,y]*dy/10000
                                endif 
			endfor
			x_axle_r= lon_new[max_x2:max_x1,near_row]
		endfor
		;convert 2D x_axle_r to 1D x_axle
                x_axle= x_axle_r[0:n_elements(x_axle_r)-1]
	
	;south-north
	endif else if (dir eq 1) or (dir eq 7) then begin
                max_y1= max( y_temp_col )
                max_y2= min( y_temp_col )
                max_x1= max( x_temp_row_down )
                max_x2= min( x_temp_row_up )
                num_y= max_y1 - max_y2 +1
                num_x= max_x1 - max_x2 +1
                ;Line Density:LD
                LD=fltarr(num_x)
                For x=max_x2 , max_x1 do begin
                        For y=max_y2 , max_y1 do begin
			    if y eq 0 then begin
				dy=abs( lon_new[y+1,x]-lon_new[y,x] )
                            endif else if y eq N_elements(lon)-1 then begin
                                dy=abs( lon_new[y,x]-lon_new[y-1,x])
			    endif else begin
                                dy=abs( (lon_new[y+1,x]-lon_new[y-1,x]) )/2
			    endelse
				if TVCD[y,x] ne -999 then begin
	                                LD[x-max_x2]+=TVCD[y,x]*dy/10000
                                endif
                        endfor
                        x_axle_r= lat_new[near_col,max_x2:max_x1]
                endfor
		;convert 2D x_axle_r to 1D x_axle
		x_axle= x_axle_r[0:n_elements(x_axle_r)-1]
	
	;southwest-northeast
	endif else if (dir eq 2) or (dir eq 6) then begin
		num_y_positive=0
	        num_y_negative=0
        	x_positive=[near_col,near_row]
	        x_negative=[near_col,near_row]
		For count=1,min([N_Elements(lon),N_Elements(lat)]) do begin
			if (not array_equal(where(y_temp_col eq (near_col-count)),-1))$
			   and(not array_equal(where(y_temp_row eq (near_row+count)),-1)) then begin
				num_y_positive+=1
			endif
			if (not array_equal(where(y_temp_col eq (near_col+count)),-1))$
                           and(not array_equal(where( y_temp_row eq (near_row-count)),-1)) then begin
                                num_y_negative+=1
			endif
		

			if (not array_equal(where(x_temp_col_down eq near_col+count-1),-1))$
                           and(not array_equal(where(x_temp_row_down eq (near_row+count)),-1)) then begin
                                x_positive=[x_positive,[near_col+count-1,near_row+count]]
			endif
			if (not array_equal(where(x_temp_col_down eq (near_col+count)),-1))$
                           and(not array_equal(where( x_temp_row_down eq (near_row+count)),-1)) then begin
                                x_positive=[x_positive,[near_col+count,near_row+count]]
                        endif
                        if (not array_equal(where(x_temp_col_up eq (near_col-count)),-1))$
                           and(not array_equal(where(x_temp_row_up eq near_row-count+1),-1)) then begin
                               x_negative=[x_negative,[near_col-count,near_row-count-1]]
			endif
			if (not array_equal(where( x_temp_col_up eq (near_col-count)),-1))$
                           and(not array_equal(where(x_temp_row_up eq (near_row-count)),-1)) then begin
                               x_negative=[x_negative,[near_col-count,near_row-count]]
                        endif
		endfor
		;print,'num_y_positive',num_y_positive,'num_y_negative',num_y_negative
		x_positive=reform(x_positive,[2,N_Elements(x_positive)/2])
		x_negative=reform(x_negative,[2,N_Elements(x_negative)/2])
		x_negative_r=reverse(x_negative,2)
		;exclude last parameter of x_negative_r in case double-counting
		x_negative_new=x_negative_r[*,0:N_elements(x_negative)/2-2]
		x_sequence=[[x_negative_new],[x_positive]]
		;print,'size of x_sequence',size(x_sequence)
		;print,'x_sequence',x_sequence
		LD=fltarr(N_elements(x_sequence)/2)
		x_axle=fltarr(N_elements(x_sequence)/2)
		For x=0,N_elements(x_sequence)/2-1 do begin
			I=x_sequence[0,x]
			J=x_sequence[1,x]
			For y=0, num_y_positive do begin
			  ;pay attention to the database boundary
			  if ( I-y ge 0) and ( J+y le  N_elements(lat)-1) then begin
			    if (I-y eq 0) and (J+y ne 0) then begin
				dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y+1,J+y-1])
			    endif else if (I-y eq 0) and (J+y eq 0) then begin
				dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y+1,J+y+1])
			    endif else if (I-y ne 0) and (J+y eq 0) then begin
				dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y-1,J+y+1])
			    endif else if (I-y eq N_elements(lon)-1) and (J+y ne N_elements(lat)-1) then begin
				dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y-1,J+y+1])
                            endif else if (I-y eq N_elements(lon)-1) and (J+y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y-1,J+y-1])
                            endif else if (I-y ne N_elements(lon)-1) and (J+y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y+1,J+y-1])
			    endif else begin				
				dy=abs( (lat_new_rot[I-y-1,J+y+1]-lat_new_rot[I-y+1,J+y-1]) )/2
			    endelse
				if TVCD[I-y,J+y] ne -999 then begin
					LD[x]+= TVCD[I-y,J+y]*dy/10000
                                endif   
			  endif
			endfor
			For y=1,num_y_negative do begin
                          ;pay attention to the database boundary
                          if (J-y ge 0) and (I+y le  N_elements(lon)-1) then begin
                            if (I+y eq 0) and (J-y ne 0) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y+1,J-y-1])
                            endif else if (I+y eq 0) and (J-y eq 0) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y+1,J-y+1])
                            endif else if (I+y ne 0) and (J-y eq 0) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y-1,J-y+1])
                            endif else if (I+y eq N_elements(lon)-1) and (J-y ne N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y-1,J-y+1])
                            endif else if (I+y eq N_elements(lon)-1) and (J-y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y-1,J-y-1])
                            endif else if (I+y ne N_elements(lon)-1) and (J-y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y+1,J-y-1])
                            endif else begin
				dy=abs( (lat_new_rot[I+y-1,J-y+1]-lat_new_rot[I+y+1,J-y-1]) )/2
                            endelse
				if TVCD[I+y,J-y] ne -999 then begin
					LD[x]+= TVCD[I+y,J-y]*dy/10000
                                endif
  
			  endif				
			endfor
			x_axle[x]=lon_new_rot[x_sequence[0,x],x_sequence[1,x]]
		endfor
        ;southeast-northwest
        endif else if (dir eq 0) or (dir eq 8) then begin
                num_y_positive=0
                num_y_negative=0
                x_positive=[near_col,near_row]
                x_negative=[near_col,near_row]
                For count=1,min([N_Elements(lon),N_Elements(lat)]) do begin
                        if (not array_equal(where(y_temp_col eq (near_col-count)),-1))$
                           and(not array_equal(where(y_temp_row eq (near_row-count)),-1)) then begin
                                num_y_positive+=1
                        endif
                        if (not array_equal(where(y_temp_col eq (near_col+count)),-1))$
                           and(not array_equal(where( y_temp_row eq (near_row+count)),-1)) then begin
                                num_y_negative+=1
                        endif


                        if (not array_equal(where(x_temp_col_down eq near_col-count),-1))$
                           and(not array_equal(where(x_temp_row_down eq (near_row+count-1)),-1)) then begin
                                x_positive=[x_positive,[near_col-count,near_row+count-1]]
                        endif
                        if (not array_equal(where(x_temp_col_down eq (near_col-count)),-1))$
                           and(not array_equal(where( x_temp_row_down eq (near_row+count)),-1)) then begin
                                x_positive=[x_positive,[near_col-count,near_row+count]]
                        endif
                        if (not array_equal(where(x_temp_col_up eq (near_col+count-1)),-1))$
                           and(not array_equal(where(x_temp_row_up eq near_row-count),-1)) then begin
                               x_negative=[x_negative,[near_col+count-1,near_row-count]]
                        endif
                        if (not array_equal(where( x_temp_col_up eq (near_col+count)),-1))$
                           and(not array_equal(where(x_temp_row_up eq (near_row-count)),-1)) then begin
                               x_negative=[x_negative,[near_col+count,near_row-count]]
                        endif
                endfor
                ;print,'num_y_positive',num_y_positive,'num_y_negative',num_y_negative
                x_positive=reform(x_positive,[2,N_Elements(x_positive)/2])
                x_negative=reform(x_negative,[2,N_Elements(x_negative)/2])
                x_negative_r=reverse(x_negative,2)
                ;exclude last parameter of x_negative_r in case double-counting
                x_negative_new=x_negative_r[*,0:N_elements(x_negative)/2-2]
                x_sequence=[[x_negative_new],[x_positive]]
                ;print,'size of x_sequence',size(x_sequence)
                ;print,'x_sequence',x_sequence
                LD=fltarr(N_elements(x_sequence)/2)
                x_axle=fltarr(N_elements(x_sequence)/2)
                For x=0,N_elements(x_sequence)/2-1 do begin
                        I=x_sequence[0,x]
                        J=x_sequence[1,x]
                        For y=0, num_y_positive do begin
                          ;pay attention to the database boundary
                          if ( I-y ge 0) and ( J-y ge 0) then begin
                            if (I-y eq 0) and (J-y ne  N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I-y,J-y]-lat_new_rot[I-y+1,J-y+1])
                            endif else if (I-y eq 0) and (J-y eq  N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I-y,J-y]-lat_new_rot[I-y+1,J-y])
                            endif else if (I-y ne N_elements(lon)-1) and (J-y eq 0) then begin
                                dy=abs(lat_new_rot[I-y,J-y]-lat_new_rot[I-y+1,J-y+1])
                            endif else if (I-y eq N_elements(lon)-1) and (J-y eq 0) then begin
                                dy=abs(lat_new_rot[I-y,J-y]-lat_new_rot[I-y,J-y+1])
                            endif else if (I-y eq N_elements(lon)-1) or (J+y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I-y,J-y]-lat_new_rot[I-y-1,J-y-1])
                            endif else begin
                                dy=abs( (lat_new_rot[I-y-1,J-y-1]-lat_new_rot[I-y+1,J-y+1]) )/2
                            endelse
			    if TVCD[I-y,J-y] ne -999 then begin
	                            LD[x]+= TVCD[I-y,J-y]*dy/10000
  			    endif
			  endif	
                        endfor

                        For y=1,num_y_negative do begin
                          ;pay attention to the database boundary
                          if (J+y le  N_elements(lat)-1) and (I+y le  N_elements(lon)-1) then begin
                            if (I+y eq 0) and (J+y ne  N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J+y]-lat_new_rot[I+y+1,J+y+1])
                            endif else if (I+y eq 0) and (J+y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J+y]-lat_new_rot[I+y+1,J+y])
                            endif else if (I+y ne  N_elements(lon)-1) and (J+y eq 0) then begin
                                dy=abs(lat_new_rot[I+y,J+y]-lat_new_rot[I+y+1,J+y+1])
                            endif else if (I+y eq N_elements(lon)-1) and (J+y eq 0) then begin
                                dy=abs(lat_new_rot[I+y,J+y]-lat_new_rot[I+y,J+y+1])
                            endif else if (I+y eq N_elements(lon)-1) or (J+y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J+y]-lat_new_rot[I+y-1,J+y-1])
                            endif else begin
                                dy=abs( (lat_new_rot[I+y+1,J+y+1]-lat_new_rot[I+y-1,J+y-1]) )/2
                            endelse
			    if TVCD[I+y,J+y] ne -999 then begin
                                LD[x]+= TVCD[I+y,J+y]*dy/10000
                            endif
 
                          endif
                        endfor
                        x_axle[x]=lon_new_rot[x_sequence[0,x],x_sequence[1,x]]
                endfor
		;print,'southeast-northwest',LD	
	endif else begin
		LD=0
		x_axle=0
	endelse
	
        if (dir eq 3) or (dir eq 6) or (dir eq 7) or (dir eq 8) then begin
                x_axle=-x_axle
        endif


	LD=LD*10.0D
	;sort x_axle ascend
	LD=LD[sort(x_axle)]
	x_axle=x_axle[sort(x_axle)]

	namein=string(lat_point)+'N'+string(lon_point)+'E'+'season'+string(season+1)+'dir'+string(dir+1)
	name=strcompress(namein,/REMOVE)
	PLOT,x_axle,LD,psym=2,$
	        title='wind direction'+string(dir+1),xtitle='x(km)',ytitle='NO2 Line Density(10^23 molec/cm)',$
		xrange=[-250,300]
	image = tvrd(true =1)
	write_jpeg,name+'LD.jpg',image,true=1

  endfor
endfor	

end
