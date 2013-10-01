pro read_no2,year,month,day,I,J,density 
;read NO2 profile from Loc_ data
Yr4 = string(year,format='(i4.4)')
mon2 = string(month, format='(i2.2)')
day2 = string(day, format='(i2.2)')
filename = file_search('/z6/satellite/OMI/no2/OMI_Lok/'+Yr4+'/'+Yr4+'_'+Mon2+'_'+Day2+'_NO2TropCS30.hdf5')

file_id=H5F_OPEN(filename)
dataset_id_ColumnAmountNO2Trop=H5D_OPEN(file_id,'/NO2.COLUMN.VERTICAL.TROPOSPHERIC.CS30_BACKSCATTER.SOLAR')
ColumnAmountNO2Trop=H5D_READ(dataset_id_ColumnAmountNO2Trop)
H5D_CLOSE,dataset_id_ColumnAmountNO2Trop
if ColumnAmountNO2Trop[I,J] gt 0U  then begin
	density=ColumnAmountNO2Trop[I,J]
endif else begin
	density=-999.0
endelse

end

pro loc_density_v1
;this program is used for fitler NO2 based on back trajector data
;in v1, the area is defined as pathway only inside stude area, sea and background (NO2 <1*10^15)

;the number of target point is 214, but the header is -9999,-9999
;&&&&&&&&&&&&&&input&&&&&&&&&&&&&&&&&
point_num=214

loc=fltarr(2,point_num)
loc_I=intarr(point_num)
loc_J=intarr(point_num)




;************************************************************************
;build up look-up table for location of target point 
header=fltarr(2,1)
filename = '/home/liufei/Data/Bay/loc_area_0.1.asc'
openr,lun,filename,/get_lun
readf,lun,header,loc
close,lun
Free_LUN,lun

nlon = 3600
nlat = 1800
lon = fltarr(nlon)
lat = fltarr(nlat)
lon = -179.95+indgen(nlon)*0.1
lat = -89.95 +indgen(nlat)*0.1

For i=0,point_num-1 do begin
	loc_I[i]= where(abs(lon-loc[1,i]) lt 0.0001)
	loc_J[i]= where(abs(lat-loc[0,i]) lt 0.001)
endfor

;print,lat
;print,loc[1,0]-lon[2982]
;print,loc_I
;print,loc_J




;************************************************************************
;read shapefile of bohai see
shapefile ='/home/liufei/Data/Bay/Bohai_2/bohai_see.shp'
oshp_see=OBJ_NEW('IDLffShape',shapefile)
oshp_see->getproperty,n_entities=n_ent_see,Attribute_info=attr_info_see,n_attributes=n_attr_see,Entity_type=ent_type_see
print,'n_entities',n_ent_see
print,'n_attributes',n_attr_see
print, attr_info_see[*].name

;;see is code_see =1
Code_see = 1

flag =1
For i=0,n_ent_see-1 Do Begin
        attr_see=oshp_see->GetAttributes(i)
        ;index=3 is GRID
        index=3
        IF total(attr_see.(index) eq Code_see) eq 1.0 then begin
            ent_see=oshp_see->getEntity(i)
            vert_see=*(ent_see.vertices)
            if flag then begin
                    Obj_see = OBJ_NEW('IDLanROI',vert_see)
                    flag =0
            endif else begin
                    Obj_see->IDLanROI::AppendData,vert_see
            endelse
        endif
endfor


;read shapefile of target area
;shapefile ='/home/liufei/Data/Bay/Prov_Boundary/Prov_Boundary.shp'
;shapefile ='/home/liufei/Data/Bay/Prov_Boundary_Reg/Prov_Boundary.shp'

shapefile ='/home/liufei/Data/Bay/Boundary/Hebei/Hebei_Add.shp'
;shapefile ='/home/liufei/Data/Bay/Boundary/HebeiN/HebeiN_Add.shp'

oshp=OBJ_NEW('IDLffShape',shapefile)
oshp->getproperty,n_entities=n_ent,Attribute_info=attr_info,n_attributes=n_attr,Entity_type=ent_type
print,'n_entities',n_ent
print,'n_attributes',n_attr
print, attr_info[*].name
;name = ['Beijing','Tianjin','Hebei','Shanxi','Liaoning','Shandong','Jiangsu','Shanghai','Zhejiang','Henan','Anhui']
;ProvCode=[11,12,13,140000,21,37,320000,310000,330000,410000,340000]
;Prov_Boundary/Prov_Boundary.shp should be changed to Prov_Boundary_Regi/Prov_Boundary.shp
;name = ['Beijing','Tianjin','HebeiN','HebeiS','Liaoning','Shandong']
;RegCode=[1,2,3,33,6,15]

name = ['Hebei','Liaoning','Shandong']
CatCode=[1,6,15]
;name = ['HebeiN','HebeiS','Liaoning','Shandong']
;CatCode=[1,2,6,15]

;&&&&&&&&&&&&&&input&&&&&&&&&&&&&&&&&
;cycle for area
;For area_ID=0,0 do begin
;case area_ID of
;0: select=[0,1,2]
;1: select=[4]
;2: select=[5]
;endcase

;For area_ID=0,3 do begin
;case area_ID of
;0: select=[0,1,2]
;1: select=[3]
;2: select=[4]
;3: select=[5]
;endcase



For area_ID=0,0 do begin
case area_ID of
0: select=[0]
1: select=[1]
2: select=[2]
endcase

;For area_ID=0,3 do begin
;case area_ID of
;0: select=[0]
;1: select=[1]
;2: select=[2]
;3: select=[3]
;endcase



name_select=name[select[0]]
;Code=ProvCode[select]
;Code=RegCode[select]
Code=CatCode[select]

flag = 1
For i=0,n_ent-1 Do Begin
        attr=oshp->GetAttributes(i)
        ;index=0 is ProvCode
        ;index=0
	;index=2 is RegCode
	;index=2

	;index=0 is CatCode
	index=0

        IF  where(attr.(index) eq Code) ge 0U then begin
            ent0=oshp->getEntity(i)
            vert0=*(ent0.vertices)
            ;help,vert0
            if flag then begin
                    Obj = OBJ_NEW('IDLanROI',vert0)
                    flag =0
            endif else begin
                    Obj->IDLanROI::AppendData,vert0
            endelse
	endif

endfor



;&&&&&&&&&&&&&&input&&&&&&&&&&&&&&&&&
For year =2005,2005 do begin


;************************************************************************
;read shapefile of background area
shapefile ='/home/liufei/Data/Bay/Background_'+string(year,format='(i4.4)')+'/background_'+string(year,format='(i4.4)')+'.shp'
oshp_b=OBJ_NEW('IDLffShape',shapefile)
oshp_b->getproperty,n_entities=n_ent_b,Attribute_info=attr_info_b,n_attributes=n_attr_b,Entity_type=ent_type_b
print,'n_entities',n_ent_b
print,'n_attributes',n_attr_b
print, attr_info_b[*].name

;;Code_b =1 background
Code_b =1
flag =1
For i=0,n_ent_b-1 Do Begin
        attr_b=oshp_b->GetAttributes(i)
        ;index=1 is GRIDCODE
        index=1
        IF total(attr_b.(index) eq Code_b) eq 1.0 then begin
            ent_b=oshp_b->getEntity(i)
            vert_b=*(ent_b.vertices)
            ;help,vert_b
	    if flag then begin
                    Obj_b = OBJ_NEW('IDLanROI',vert_b)
                    flag =0
            endif else begin
                    Obj_b->IDLanROI::AppendData,vert_b
            endelse

        endif
endfor


;************************************************************************
;&&&&&&&&&&&&&&input&&&&&&&&&&&&&&&&&
;hour_num is used for display the run time: 24h? 48h? 72h?
hour_num = 24

;density is used to yearly mean NO2
density= fltarr(point_num)
;sample is used to sample number of satellite data in one year
sample = fltarr(point_num)
;sample_traj is used to sample number of trajector data in one year
sample_traj = fltarr(point_num)

;&&&&&&&&&&&&&&input&&&&&&&&&&&&&&&&&
  month_start=1
  n_month=12
  ;summary mark the data [trajector number, satellite number, NO2 Average]
  summary = fltarr(3,n_month)
  For month = month_start,month_start+n_month-1 do begin

    ;density_m is used to monthly mean NO2
    density_m= fltarr(point_num)
    ;sample_m is used to sample number of satellite data in one month
    sample_m = fltarr(point_num)
    ;sample_traj_m is used to sample number of trajector data in one month
    sample_traj_m = fltarr(point_num)

;&&&&&&&&&&&&&&input&&&&&&&&&&&&&&&&&
    daymonth=[31,28,31,30,31,30,31,31,30,31,30,31]
    For day = 1,daymonth[month-1] do begin
	Yr2 = strmid(string(year,format='(i4.4)'),2,2)
	mon2 = string(month, format='(i2.2)')
	day2 = string(day, format='(i2.2)')	
	print,'process',Yr2+mon2+day2
	;read data from back trajecto
	data=fltarr(14,point_num*(hour_num+1))
	;first is used to read the line of GDAS
	first=intarr(2)
	filename = '/home/liufei/Data/Bay/Traj_data/'+string(hour_num,format='(i2.2)')+'h/'+string(year,format='(i4.4)')+'/'+Yr2+mon2+day2+'14'
	openr,lun, filename,/GET_LUN
	readf,lun,first
	header=strarr(first[0]+2+point_num)
        header[*]=''
	readf,lun,header
	readf,lun,data
	close,lun
	free_lun,lun
	;print,data
	
	;get the trajector of target point [lon,lat]
	traj=fltarr(2,hour_num,point_num)
	For k=0,point_num-1 do begin
		For n=0, hour_num-1 do begin
		    traj[1,n,k]= data[9,k+(n+1)*point_num ]
		    traj[0,n,k]= data[10,k+(n+1)*point_num]
		endfor
	endfor

	;find out the trajector is inside or outside target area+background area
	;result is used to mark the trajector is inside  or outside the area
	;0 = Exterior. The point lies strictly out of bounds of the ROI.
	;1 = Interior. The point lies strictly inside the bounds of the ROI.
	;2 = On edge. The point lies on an edge of the ROI boundary.
	;3 = On vertex. The point matches a vertex of the ROI.
	result0 = fltarr(hour_num,point_num)
	result1 = fltarr(hour_num,point_num)
	result2 = fltarr(hour_num,point_num)
	For k=0,point_num-1 do begin
		result0[*,k] = Obj->ContainsPoints(traj[*,*,k])
		result1[*,k] = Obj_see->ContainsPoints(traj[*,*,k])
		result2[*,k] = Obj_b->ContainsPoints(traj[*,*,k])
	endfor
	result_final = fltarr(hour_num,point_num)
	For k=0,point_num-1 do begin
		if not array_equal(where(result0[*,k] ge 1),[-1]) then begin
			temp= where( (result0[*,k] ge 1) or (result1[*,k] ge 1) or (result2[*,k] ge 1) )
      		        if not array_equal(temp,[-1]) then begin
                	        result_final[temp,k] = 1
          	        endif
		endif
	endfor
				
			
	;find out the point with trajector inside area,then record their no2
	For k=0,point_num-1 do begin
;&&&&&&&&&&&&&&input&&&&&&&&&&&&&&&&&
                if (hour_num - total(result_final[*,k])) le ceil(hour_num*0.2) then begin
			sample_traj_m[k]+=1
			I=loc_I[k]
			J=loc_J[k]
			;PRINT,I,J
			temp=0
			read_no2,year,month,day,I,J,temp
			if temp gt 0U then begin
				density_m[k]+= temp		
				sample_m[k]+=1
			endif
		endif
	endfor



    endfor

    For k=0,point_num-1 do begin
	if sample_m[k] gt 0U then begin
		density_m[k]=density_m[k]/sample_m[k]
	endif
    endfor

    ;print,'traj:', sample_traj_m
    ;print,'satellite:',sample_m
    print,name_select,':number of trajector sample in ', Yr2+Mon2, total(sample_traj_m)
    print,name_select,':number of satellite sample in ', Yr2+Mon2, total(sample_m)
    ;print,density_m
    summary[*,month-month_start] =[total(sample_traj_m),total(sample_m),mean(density_m[where(density_m gt 0U)]) ]
    
    ;output sample_traj_m, sample_m, density_m and convert them to global map
    month1=fltarr(nlon,nlat)
    month2=fltarr(nlon,nlat)
    month3=fltarr(nlon,nlat)

    For k=0,point_num-1 do begin
        I=loc_I[k]
        J=loc_J[k]
        month1[I,J]=sample_traj_m[k]
        month2[I,J]=sample_m[k]
        month3[I,J]=density_m[k]
    endfor
    ;output month1,month2,month3 into nc
    ID_list=lonarr(n_elements(name))
    Out_nc_file='/home/liufei/Data/Bay/Result/'+string(hour_num,format='(i2.2)')+'h/v1/'+name_select+Yr2+mon2+'.nc'
    FileId = NCDF_Create( Out_nc_file, /Clobber )
    NCDF_Control, FileID, /NoFill
    xID   = NCDF_DimDef( FileID, 'X', nlon )
    yID   = NCDF_DimDef( FileID, 'Y', nlat )
    LonID = NCDF_VarDef( FileID, 'x',    [xID],     /Float )
    LatID = NCDF_VarDef( FileID, 'y',    [yID],     /Float )
    num_trajID= NCDF_VarDef( FileID,'Num_Traj',  [xID,yID], /Long )
    num_satID   = NCDF_VarDef( FileID,'Num_Sat',   [xID,yID], /Long )
    NO2ID= NCDF_VarDef( FileID,'NO2',  [xID,yID], /Double )
    
    NCDF_Attput, FileID, /Global, 'Title', 'name_select+Yr2+mon2'
    NCDF_Control, FileID, /EnDef
    NCDF_VarPut, FileID, LonID, lon ,   Count=[ nlon ]
    NCDF_VarPut, FileID, LatID, lat ,   Count=[ nlat ]
    NCDF_VarPut, FileID, num_trajID , month1, Count=[ nlon,nlat ]
    NCDF_VarPut, FileID, num_satID , month2, Count=[ nlon,nlat ]
    NCDF_VarPut, FileID, NO2ID , month3, Count=[ nlon,nlat ]
    NCDF_Close, FileID



  endfor;month
  OBJ_DESTROY,oshp_b
  OBJ_DESTROY,Obj_b

  header_output=['trajector_number','satellite_number','NO2_average']
  outfile ='/home/liufei/Data/Bay/Result/'+string(hour_num,format='(i2.2)')+'h/v1/'+name_select+Yr2+'_summary.asc'
  openw,lun,outfile,/get_lun
  printf,lun,header_output,summary
  close,lun
  free_lun,lun

endfor;year


endfor ;area_ID

OBJ_DESTROY,oshp
OBJ_DESTROY,Obj

OBJ_DESTROY,oshp_see
OBJ_DESTROY,Obj_see


end     

