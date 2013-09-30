pro no2_column_region_average_loc

;read shapefile of target area
;shapefile ='/home/liufei/Data/Bay/Prov_Boundary/Prov_Boundary.shp'
;shapefile ='/home/liufei/Data/Bay/Prov_Boundary_Reg/Prov_Boundary.shp'

shapefile ='/home/liufei/Data/Bay/Boundary/Hebei/Hebei.shp'
;shapefile ='/home/liufei/Data/Bay/Boundary/HebeiN/HebeiN.shp'

oshp=OBJ_NEW('IDLffShape',shapefile)
oshp->getproperty,n_entities=n_ent,Attribute_info=attr_info,n_attributes=n_attr,Entity_type=ent_type
print,'n_entities',n_ent
print,'n_attributes',n_attr
print, attr_info[*].name
;name = ['Beijing','Tianjin','Hebei','Shanxi','Liaoning','Shandong','Jiangsu','Shanghai','Zhejiang','Henan','Anhui']
;ProvCode=[11,12,13,140000,21,37,320000,310000,330000,410000,340000]
;Prov_Boundary.shp should be changed to Prov_Boundary_Reg.shp
;name = ['Beijing','Tianjin','HebeiN','HebeiS','Liaoning','Shandong']
;RegCode=[1,2,3,33,6,15]

name = ['Hebei','Liaoning','Shandong']
CatCode=[1,6,15]
;name = ['HebeiN','HebeiS','Liaoning','Shandong']
;CatCode=[1,2,6,15]

;&&&&&&&&&&&&&&input&&&&&&&&&&&&&&&&&
;cycle for area
;For area_ID=0,2 do begin
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



For area_ID=0,2 do begin
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

flag=1
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

;read NO2 Column
nlon = 640
nlat = 400
lon = fltarr(nlon)
lat = fltarr(nlat)
lon = 72+indgen(nlon)*0.1
lat = 55-indgen(nlat)*0.1
density1=fltarr(nlon,nlat)
density2=fltarr(nlon,nlat)

mark=fltarr(nlon,nlat)

header=strarr(1,6)
filename = '/home/liufei/Data/Bay/2005_anual_Average_China_0.1X0.1.asc'
openr,lun,filename,/get_lun
readf,lun,header,density1
close,lun
Free_LUN,lun

filename = '/home/liufei/Data/Bay/2009_anual_Average_China_0.1X0.1.asc'
openr,lun,filename,/get_lun
readf,lun,header,density2
close,lun
Free_LUN,lun


For I=0,nlon-1 do begin
	For J=0,nlat-1 do begin
		mark[I,J] = Obj->ContainsPoints(lon[I],lat[J])
	endfor
endfor

print,name_select,'density1 & density2',mean(density1[where(mark eq 1)]),mean(density2[where(mark eq 1)])

OBJ_DESTROY,Obj
endfor;area_ID

OBJ_DESTROY,oshp

end
