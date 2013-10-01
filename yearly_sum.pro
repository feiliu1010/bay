pro yearly_sum

hour_num=24

nlon = 3600
nlat = 1800
lon = fltarr(nlon)
lat = fltarr(nlat)
lon = -179.95+indgen(nlon)*0.1
lat = -89.95 +indgen(nlat)*0.1

name = ['HebeiN','HebeiS','Liaoning','Shandong','Hebei']
For area_ID=2,4 do begin
case area_ID of
0: select=[0]
1: select=[1]
2: select=[2]
3: select=[3]
4: select=[4]
endcase
name_select=name[select[0]]

For year =2009,2009 do begin

Yr2 = strmid(string(year,format='(i4.4)'),2,2)
;density is used to yearly mean NO2
density= fltarr(nlon,nlat)
;month_num is number of month where NO2>0
month_num= fltarr(nlon,nlat)
;sample is used to sample number of satellite data in one year
sample = fltarr(nlon,nlat)
;sample_traj is used to sample number of trajector data in one year
sample_traj = fltarr(nlon,nlat)

;&&&&&&&&&&&&&&input&&&&&&&&&&&&&&&&&
  month_start=1
  n_month=12
  For month = month_start,month_start+n_month-1 do begin

    mon2 = string(month, format='(i2.2)')
    ;density_m is used to monthly mean NO2
    density_m= fltarr(nlon,nlat)
    ;sample_m is used to sample number of satellite data in one month
    sample_m = fltarr(nlon,nlat)
    ;sample_traj_m is used to sample number of trajector data in one month
    sample_traj_m = fltarr(nlon,nlat)

    filename = '/home/liufei/Data/Bay/Result/'+string(hour_num,format='(i2.2)')+'h/v2/'+name_select+Yr2+mon2+'.nc'
    fid=NCDF_OPEN(filename)
    varid1=NCDF_VARID(fid,'Num_Traj')
    NCDF_VARGET, fid, varid1, sample_traj_m
    varid2=NCDF_VARID(fid,'Num_Sat')
    NCDF_VARGET, fid, varid2,sample_m
    varid3=NCDF_VARID(fid,'NO2')
    NCDF_VARGET, fid, varid3,density_m
    NCDF_CLOSE, fid

    sample_traj+=sample_traj_m
    sample+=sample_m
    density+=density_m
    month_num[where(density_m gt 0U)]+=1

  endfor
  
  density[where(month_num gt 0U)]=density[where(month_num gt 0U)]/month_num[where(month_num gt 0U)]
  density[where(month_num le 0U)]=-999

  Out_nc_file='/home/liufei/Data/Bay/Result/'+string(hour_num,format='(i2.2)')+'h/v2/'+name_select+Yr2+'.nc'
  FileId = NCDF_Create( Out_nc_file, /Clobber )
  NCDF_Control, FileID, /NoFill
  xID   = NCDF_DimDef( FileID, 'X', nlon )
  yID   = NCDF_DimDef( FileID, 'Y', nlat )
  LonID = NCDF_VarDef( FileID, 'x',    [xID],     /Float )
  LatID = NCDF_VarDef( FileID, 'y',    [yID],     /Float )
  num_trajID= NCDF_VarDef( FileID,'Num_Traj',  [xID,yID], /Long )
  num_satID   = NCDF_VarDef( FileID,'Num_Sat',   [xID,yID], /Long )
  NO2ID= NCDF_VarDef( FileID,'NO2',  [xID,yID], /Double )

  NCDF_Attput, FileID, /Global, 'Title', 'name_select+Yr2'
  NCDF_Control, FileID, /EnDef
  NCDF_VarPut, FileID, LonID, lon ,   Count=[ nlon ]
  NCDF_VarPut, FileID, LatID, lat ,   Count=[ nlat ]
  NCDF_VarPut, FileID, num_trajID , sample_traj, Count=[ nlon,nlat ]
  NCDF_VarPut, FileID, num_satID , sample, Count=[ nlon,nlat ]
  NCDF_VarPut, FileID, NO2ID , density, Count=[ nlon,nlat ]
  NCDF_Close, FileID    

endfor;year
endfor;area_ID

end
