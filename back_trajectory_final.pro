pro back_trajectory_final


TotalGridTE=fltarr(700,700)
xmid = fltarr(3600)
ymid = fltarr(1800)
for lonnum=0,3599 do begin
  xmid[lonnum]=-179.95+0.1*lonnum
  endfor
for latnum=0,1799 do begin
  ymid[latnum]=-89.95+0.1*latnum
  endfor

x1=where( ( xmid ge 70.0) and (xmid le 70+0.05)) 
x2=where( ( xmid ge 140-0.05) and (xmid le 140)) 
y1=where( ( ymid ge 0.0) and (ymid le 0.0+0.1)) 
y2=where( ( ymid ge 70-0.05) and (ymid le 70))

transxmid=fltarr(700)
transymid=fltarr(700)
for transx=x1[0],x2[0] do begin
  transxmid[transx-x1[0]]=xmid[transx]
  endfor
for transy=y1[0],y2[0] do begin
  transymid[transy-y1[0]]=ymid[transy]
  endfor
    
hdffile='D:\IDL\ENVI5.0\IDL82\IDLWorkspace82\zyx\back_trajectoy\deposition_parameters\parameter_for_Lu_201006.hdf'

sd_id=HDF_SD_START(hdfFile)

sds0_id=HDF_SD_SELECT(sd_id,0)
HDF_SD_GETDATA,sds0_id,lsfdata
kwdata=lsfdata*1.5e-3
HDF_SD_ENDACCESS,sds0_id


sds9_id=HDF_SD_SELECT(sd_id,9)
HDF_SD_GETDATA,sds9_id,boxheight
HDF_SD_ENDACCESS,sds9_id


sds10_id=HDF_SD_SELECT(sd_id,10)
HDF_SD_GETDATA,sds10_id,PBLheight
HDF_SD_ENDACCESS,sds10_id



HDF_SD_END,sd_id


hdfxmid = fltarr(144)
hdfymid = fltarr(91)
for hdflonnum=0,143 do begin
  hdfxmid[hdflonnum]=-178.25+2.5*hdflonnum
  endfor
for hdflatnum=0,90 do begin
  hdfymid[hdflatnum]=-89.0+2.0*hdflatnum
  endfor
    
daynumber=indgen(31)
daysub=daynumber[1:30]
dayjune=daynumber[daysub]
day=string(dayjune,format='(I02)')
hour=['00','06','12','18']
for trjday=0,29 do begin
  for trjhour=0,3 do begin
    
filename='D:\IDL\ENVI5.0\IDL82\IDLWorkspace82\zyx\back_trajectoy\raw data\'+'tdump1006'+day[trjday]+hour[trjhour]

 line='' 
openr,iluns, filename,/GET_LUN 
readf,iluns,line
readf,iluns,line
readf,iluns,line
readf,iluns,line
readf,iluns,line
readf,iluns,line
readf,iluns,line
readf,iluns,line
readf,iluns,line
readf,iluns,line
readf,iluns,line

trajpoint=fltarr(13)
linenumber=-1L
lonmax=0L
latmax=0L
lonmin=4000L
latmin=2000L
timedata=fltarr(73,700,700)
timedata[*,*,*]=999.0
gridtime=fltarr(700,700)
gridtime[*,*]=999.0
timediff=fltarr(700,700)
heightdata=fltarr(73,700,700)
gridheight=fltarr(700,700)
gridlon=lonarr(73)
gridlat=lonarr(73)
timestorage=fltarr(300)
timestorage[*]=999.0


while ~eof(iLuns) do begin
  
ReadF, iluns, trajpoint

lon = where( ( transxmid ge trajpoint[10]-0.05) and (transxmid le trajpoint[10]+0.05) ,count1)
lat = where( ( transymid ge trajpoint[9]-0.05) and (transymid le trajpoint[9]+0.05) ,count2)
  if ( count1 eq 2 ) then begin
     lon = min(lon)
   endif
  if ( count2 eq 2) then begin
    lat = min(lat)
  endif

  if (lon gt lonmax) then begin
    lonmax=lon
  endif
  if (lat gt latmax) then begin
    latmax=lat
  endif
  if (lon lt lonmin) then begin
    lonmin=lon
  endif
  if (lat lt latmin) then begin
    latmin=lat
  endif

linenumber+=1 
id=linenumber
timedata[id,lon,lat]=trajpoint[8]
heightdata[id,lon,lat]=trajpoint[11]
gridlon[id]=lon[0]
gridlat[id]=lat[0]
endwhile

close,iluns
Free_LUN,iluns

for x=lonmin[0],lonmax[0] do begin
  for y=latmin[0],latmax[0] do begin
    number=0
    totaltime=0.0
    totalheight=0.0
    
    for i=0,72 do begin
      value1=timedata[i,x,y]
      value2=heightdata[i,x,y]
      if(value1 le 0) then begin
        number+=1
        totaltime+=timedata[i,x,y]
        totalheight+=heightdata[i,x,y]
       endif
    endfor
   
       
       if (number ge 2) then begin
        timediff[x,y]=number-1.0
        endif else begin
        timediff[x,y]=number+0.0
        endelse
       if (number ne 0) then begin
        gridtime[x,y]=totaltime/number
        gridheight[x,y]=totalheight/number
        endif
       
   endfor
endfor

;interpolation when the adjacent two points are not in adjacent two grids
for d=0,71 do begin
  gn1=gridlon[d] 
  gn2=gridlon[d+1]
  gt1=gridlat[d]
  gt2=gridlat[d+1]
  if((abs((gn2-gn1)gt 1)) or (abs((gt2-gt1)gt 1))) then begin
    xinsert=abs((gn2-gn1))-1 ;the numbers of gid inserted in the longitude direction 
    yinsert=abs((gt2-gt1))-1 ;the numbers of gid inserted in the latitude direction
    gridinsert=(xinsert ge yinsert)?xinsert:yinsert ;the numbers of inserted gid  
    lowgn=(gn1 lt gn2)?gn1:gn2
    lowgt=(gt1 lt gt2)?gt1:gt2
    lowtime=(timedata[d,gn1,gt1] lt timedata[d+1,gn2,gt2])?timedata[d,gn1,gt1]:timedata[d+1,gn2,gt2]
    for z=1,gridinsert do begin
      m=floor(lowgn+xinsert*z/gridinsert);the longitude index of each inserted grid 
      n=floor(lowgt+yinsert*z/gridinsert);the latitude index of each inserted grid   
      gridheight[m,n]=(heightdata[d,gn1,gt1]+heightdata[d+1,gn2,gt2])/2;the height of trjectory in the inserted grids
      
      timediff[m,n]=1.0/(gridinsert+1)
      gridtime[m,n]=lowtime+z*timediff[m,n]
      
     endfor
   endif
 endfor
 
 ;筛选出轨迹栅格，并赋予每个栅格时间
w=-1L
for u=lonmin[0],lonmax[0] do begin
  for v=latmin[0],latmax[0] do begin
    if(gridtime[u,v] le 0) then begin
      w+=1
      timestorage[w]=gridtime[u,v]    
      endif
  endfor
 endfor
 
 ;时间由小到大排序 
 trjstorage=fltarr(w+1)
 for s=0,w do begin
    trjstorage[s]=timestorage[s]
    endfor   
 
 timeorder=fltarr(w+1)
 timeorder=trjstorage[sort(trjstorage)]    
 

;搜索每个时间所对应的栅格，从而将轨迹栅格按轨迹时间从小到大排列
 rt=fltarr(w+1)
 th=fltarr(w+1)
 hdfheight=fltarr(47)
 kw=fltarr(w+1)
 PBL=fltarr(w+1)
for p=lonmin[0],lonmax[0] do begin
  for q=latmin[0],latmax[0] do begin
     for r=0,w do begin
      if (timeorder[r] eq gridtime[p,q] ) then begin  
       rt[r]=timediff[p,q]*3600
       th[r]=gridheight[p,q]
       hdflon = where( ( hdfxmid ge transxmid[p]-1.25) and ( hdfxmid le transxmid[p]+1.25) ,count1)
       hdflat = where( ( hdfymid ge transymid[q]-1.0) and (hdfymid le transymid[q]+1.0) ,count2)
        if ( count1 eq 2 ) then begin
          hdflon = min(hdflon)
         endif
        if ( count2 eq 2) then begin
          hdflat = min(hdflat)
         endif

       hdfheight=boxheight[*,hdflat,hdflon]
       hdfhabs=abs(hdfheight-th[r])
       hdfhsub=sort(hdfhabs);hdfnabs从小到大排序数组的下标索引
       hdfh=hdfhsub[0];轨迹高度最接近boxheight的下标索引
       kw[r]=kwdata[hdfh,hdflat,hdflon]
       PBL[r]=PBLheight[hdflat,hdflon]
       endif
      endfor
  endfor
endfor 


;筛选出高度小于PBL的栅格作为排放栅格，计算其传输效率
kd=4.25e-7
kc=1.01e-5
FHPO=fltarr(w+1)
FHPI=fltarr(w+1)
Fdep=fltarr(w+1)

GridTE=fltarr(700,700)
TEstorage=fltarr(w+1)
TElon=fltarr(w+1)
TElat=fltarr(w+1)
TEnumber=-1L
for g=1,w do begin
   if (th[g-1] le PBL[g-1]) then begin
      TEnumber+=1
      for h=g,w do begin
        FHPO[g-1]=0.8
        FHPI[g-1]=0.2
        FHPO[h]=FHPO[h-1]*exp(-(kd+kc)*rt[h])
        FHPI[h]=(FHPI[h-1]-kc*FHPO[h-1]/(kw[h]-kc))*exp(-(kd+kw[h])*rt[h])+(kc*FHPO[h-1]/(kw[h]-kc))*exp(-(kd+kc)*rt[h])
        Fdep[h]=FHPO[h-1]+FHPI[h-1]-FHPO[h]-FHPI[h]
       endfor
       for j=lonmin[0],lonmax[0] do begin
         for k=latmin[0],latmax[0] do begin
          
            if (timeorder[g-1] eq gridtime[j,k] ) then begin
              GridTE[j,k]=(FHPO[w]+FHPI[w])+Fdep[w]
              TotalGridTE[j,k]+=GridTE[j,k]
              TEstorage[TEnumber]=(FHPO[w]+FHPI[w])+Fdep[w]
              TElon[TEnumber]=transxmid[j]
              TElat[TEnumber]=transymid[k]
              endif
           
          endfor
        endfor
     endif
endfor

for lonm=lonmin[0],lonmax[0] do begin
  for latm=latmin[0],latmax[0] do begin
          
            if (timeorder[w] eq gridtime[lonm,latm] ) then begin
              TEn=TEnumber+1
              GridTE[lonm,latm]=1.0
              TotalGridTE[lonm,latm]+=GridTE[lonm,latm]
              TEstorage[TEn]=1.0
              TElon[TEn]=transxmid[lonm]
              TElat[TEn]=transymid[latm]
              endif
    endfor
endfor           
               
  endfor
endfor;输入120各文件的循环   

TNT=120
nb=0

TED=fltarr(700,700)
for lonmap=0,699 do begin
  for latmap=0,699 do begin
    TED[lonmap,latmap]=TotalGridTE[lonmap,latmap]/TNT  
    
    endfor
  endfor 
    
    
TEDGIS=reverse(TED,2)

Out_TED_file = 'D:\IDL\ENVI5.0\IDL82\IDLWorkspace82\zyx\back_trajectoy\outfile\TED.asc'
 ; open 
 openw, ilun_asc, Out_TED_file ,/GET_LUN
 ; write array into it
 printf, ilun_asc, TEDGIS
 ; Free ilun
 Free_LUN, ilun_asc
 ; close file
 close, ilun_asc
 
 

 
ascfile='D:\IDL\ENVI5.0\IDL82\IDLWorkspace82\zyx\back_trajectoy\emission\2010_BC.asc'

 line='' 
emissiondata = fltarr(140,140) 
openr,iluns, ascfile,/GET_LUN 
readf,iluns,line
readf,iluns,line
readf,iluns,line
readf,iluns,line
readf,iluns,line 
readf,iluns,line

readf,iluns,emissiondata 

close,iluns
Free_LUN,iluns

emissionIDL=reverse(emissiondata,2)

for emln=0,139 do begin
  for emlt=0,139 do begin
    if (emissionIDL[emln,emlt] lt 0) then begin
      emissionIDL[emln,emlt]=0.0
      endif
   endfor
endfor

emissionxmid = fltarr(140)
emissionymid = fltarr(140)
for emissionlon=0,139 do begin
  emissionxmid[emissionlon]=113.025+0.05*emissionlon
  endfor
for emissionlat=0,139 do begin
  emissionymid[emissionlat]=36.025+0.05*emissionlat
  endfor

emission=fltarr(700,700)
for emissionlonnum=0,139 do begin
  for emissionlatnum=0,139 do begin  
       emlon = where( ( transxmid ge emissionxmid[emissionlonnum]-0.05) and ( transxmid le emissionxmid[emissionlonnum]+0.05) ,count1)
       emlat = where( (transymid ge emissionymid[emissionlatnum]-0.05) and (transymid le emissionymid[emissionlatnum]+0.05) ,count2)
        if ( count1 eq 2 ) then begin
          emlon = min(hdflon)
         endif
        if ( count2 eq 2) then begin
          emlat = min(hdflat)
         endif
       emission[emlon,emlat]+=emissionIDL[emissionlonnum,emissionlatnum]
     endfor
  endfor


EEI=fltarr(700,700)
for EEIlon=0,699 do begin
  for EEIlat=0,699 do begin
    EEI[EEIlon,EEIlat]=TED[EEIlon,EEIlat]* emission[EEIlon,EEIlat] 
    
    endfor
  endfor 
    
    
EEIGIS=reverse(EEI,2)

 Out_EEI_file = 'D:\IDL\ENVI5.0\IDL82\IDLWorkspace82\zyx\back_trajectoy\outfile\EEI.asc'
 ; open 
 openw, ilun_asc, Out_EEI_file ,/GET_LUN
 ; write array into it
 printf, ilun_asc, EEIGIS
 ; Free ilun
 Free_LUN, ilun_asc
 ; close file
 close, ilun_asc
 
 
end


