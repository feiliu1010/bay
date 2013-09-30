pro MEIC_average

;read shapefile of target area

shapefile ='/home/liufei/Data/Bay/Boundary/HebeiN/HebeiN_Add.shp'

oshp=OBJ_NEW('IDLffShape',shapefile)
oshp->getproperty,n_entities=n_ent,Attribute_info=attr_info,n_attributes=n_attr,Entity_type=ent_type
print,'n_entities',n_ent
print,'n_attributes',n_attr
print, attr_info[*].name

name = ['HebeiN','HebeiS','Liaoning','Shandong']
CatCode=[1,2,6,15]

;&&&&&&&&&&&&&&input&&&&&&&&&&&&&&&&&

For area_ID=0,3 do begin
case area_ID of
0: select=[0]
1: select=[1]
2: select=[2]
3: select=[3]
endcase


name_select=name[select[0]]
Code=CatCode[select]

flag=1
For i=0,n_ent-1 Do Begin
        attr=oshp->GetAttributes(i)
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

;read NOx Emissions from MEIC
nlon = 800
nlat = 500
lon = fltarr(nlon)
lat = fltarr(nlat)
lon = 70+indgen(nlon)*0.1
lat = 10+indgen(nlat)*0.1

mark=fltarr(nlon,nlat)
For I=0,nlon-1 do begin
        For J=0,nlat-1 do begin
                mark[I,J] = Obj->ContainsPoints(lon[I],lat[J])
	endfor
endfor    

mark2=fltarr(nlon,nlat)
mark2=mark

for j=1, nlat/2 do begin
        tmp = mark2[*,j-1]
        mark2[*,j-1] = mark2[*,nlat-j]
        mark2[*,nlat-j] = tmp
endfor
undefine, tmp


header_output = [['ncols 800'],['nrows 500'],['xllcorner 70'],['yllcorner 10'],['cellsize 0.1'],['nodata_value -999.0']]
outfile ='/home/liufei/Data/Bay/Result/'+name_select+'_mark.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,mark2
close,lun
free_lun,lun

For caseID=0,1 do begin

case caseID of 
0: year=2005
1: year=2009
endcase
Yr4= string(year,format='(i4.4)')

density1=dblarr(nlon,nlat)
density2=dblarr(nlon,nlat)
density3=dblarr(nlon,nlat)
density4=dblarr(nlon,nlat)
sum=dblarr(nlon,nlat)

filename ='/home/liufei/Data/Bay/MEIC_Emission/'+Yr4+'/'+Yr4+'__residential__NOx.nc'
fid=NCDF_OPEN(filename)
varid1=NCDF_VARID(fid,'z')
NCDF_VARGET, fid, varid1, density1
NCDF_CLOSE, fid

filename ='/home/liufei/Data/Bay/MEIC_Emission/'+Yr4+'/'+Yr4+'__industry__NOx.nc'
fid=NCDF_OPEN(filename)
varid2=NCDF_VARID(fid,'z')
NCDF_VARGET, fid, varid2, density2
NCDF_CLOSE, fid

filename ='/home/liufei/Data/Bay/MEIC_Emission/'+Yr4+'/'+Yr4+'__power__NOx.nc'
fid=NCDF_OPEN(filename)
varid3=NCDF_VARID(fid,'z')
NCDF_VARGET, fid, varid3, density3
NCDF_CLOSE, fid

filename ='/home/liufei/Data/Bay/MEIC_Emission/'+Yr4+'/'+Yr4+'__transportation__NOx.nc'
fid=NCDF_OPEN(filename)
varid4=NCDF_VARID(fid,'z')
NCDF_VARGET, fid, varid4, density4
NCDF_CLOSE, fid

;z:nodata_value = -9999.
density1[where(density1 lt 0)]=0
density2[where(density2 lt 0)]=0
density3[where(density3 lt 0)]=0
density4[where(density4 lt 0)]=0

sum=density1+density2+density3+density4
outfile ='/home/liufei/Data/Bay/Result/'+Yr4+'_MEIC_Emissions.asc'
openw,lun,outfile,/get_lun
printf,lun,header_output,sum
close,lun
free_lun,lun

print,'total of MEIC', total(sum)
print,name_select,'sum',Yr4,total(sum[where(mark eq 1)])
endfor;caseID

OBJ_DESTROY,Obj
endfor;area_ID

OBJ_DESTROY,oshp

end
