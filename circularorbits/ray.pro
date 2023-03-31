pro  ray

theta0=0.
phi0=0
WINSIZE=20
bkg='wwhite'
backg='ny'
nm=1uL
ln=20001ul
nmstart=0l
lee=1ul
lbb=0ul
;@nlff_block

if not keyword_set(bkg) then bkg='white' else bkg=bkg(0)

 Openr,xcoor,'rayx1.txt',/Get_Lun
  Point_lun,xcoor,0
  xcoord=fltarr(ln,nm)
  ReadF,xcoor,xcoord
  
   Openr,ycoor,'rayy1.txt',/Get_Lun
  Point_lun,ycoor,0
  ycoord=fltarr(ln,nm)
  ReadF,ycoor,ycoord
  
   Openr,zcoor,'rayz1.txt',/Get_Lun
  Point_lun,zcoor,0
  zcoord=fltarr(ln,nm)
  ReadF,zcoor,zcoord

   ;Openr,xcoor,'rayx3.txt',/Get_Lun
  ;Point_lun,xcoor,0
  ;xcoord3=fltarr(ln,nm)
  ;ReadF,xcoor,xcoord3

  ; Openr,ycoor,'rayy3.txt',/Get_Lun
  ;Point_lun,ycoor,0
  ;ycoord3=fltarr(ln,nm)
  ;ReadF,ycoor,ycoord3

  ; Openr,zcoor,'rayz3.txt',/Get_Lun
  ;Point_lun,zcoor,0
 ;zcoord3=fltarr(ln,nm)
  ;ReadF,zcoor,zcoord3
  
  ;Openr,lunAo102,'bg.txt',/Get_Lun
  ;Point_lun,lunAo102,0
  disk2=fltarr(801,801)
  ;ReadF,lunAo102,disk2
   

  bz=disk2
  ;bz=background

  
  free_lun,xcoor,ycoor,zcoor;,lunAo102
;print,xcoord(*,73),ycoord(*,73),'zzzzzzzzzzz=',zcoord(*,73)


  length=4
xmin=-length
xmax=length
ymin=-length
ymax=length
zmin=-length
zmax=length


if(phi0 lt 0) then begin
   phi1=phi0+360
   phi=!dtor*phi1
endif else begin
   phi=phi0*!dtor
endelse
theta=theta0*!dtor

st=sin(theta) & ct=cos(theta)
sp=sin(phi)   & cp=cos(phi)

ok=0
aga_flag=0
addnum=''

oldn=!D.name & set_plot,'ps'
 
loadct,0  ;  loadct,3 also looks nice too
    

l=16  & xxss=l & yyss=l 

device,filename='./raylinesp_1.ps',color=1,xsize=xxss,ysize=yyss,$
   xoff=(21.-xxss)/2. , yoff=(21-yyss)/2.,bits_per_pixel=8;decomposed=0


 

;plot,[-1,1],[-1,1],$
;pos=[0,1,0,1],$
;    /noerase,/device,/nodata,$
;    xsty=5,ysty=5,xr=[xmin,xmax],yr=[ymin,ymax]
 
	xlen=0.8
	ylen=0.8
	xslitlen=0.
	yslitlen=0.
	mx=1 ;figure numbers on horizon direction
	my=1 ;figure numbers on vertical direction.
	deltax=[(xlen-xslitlen)/mx+xslitlen,0,(xlen-xslitlen)/mx+xslitlen,0]
	deltay=[0,(ylen-yslitlen)/my+yslitlen,0,(ylen-yslitlen)/my+yslitlen]
	deltxy=[xlen/mx,ylen/my,xlen/mx,ylen/my]
	pos_llft=[(1.-xlen)/2.,(1.-ylen)/2.,(1.-xlen)/2.+(xlen-xslitlen)/mx,(1.-ylen)/2.+(ylen-yslitlen)/my]

        nx=0 & ny=0
		plot,[xmin,xmax],[ymin,ymax],pos=[pos_llft+deltax*nx+deltay*ny],xrange=[xmin,xmax],yrange=[ymin,ymax],$
		/noerase,/device,/ynozero,/normal,xstyle=4+1,ystyle=4+1,charsize=0.5,xtickv=[0,0.5,1,1.5],$
		xticks=10,xminor=10,xtickname=replicate(' ',10),/nodata
 
   pos=[pos_llft+deltax*nx+deltay*ny]
 im = bytscl(bz)
xs=1 & ys=1
;tvscl,im,pos[0],pos[1],xsize=xlen,ysize=ylen,/norm

loadct,3
tvlct,re,gr,bl,/get
re(249:255)=[0b,255b,0b,0b,255b,0b,255b]
gr(249:255)=[0b,0b,255b,0b,255b,255b,255b]
bl(249:255)=[0b,0b,0b,255b,0b,255b,255b]
tvlct,re,gr,bl
 

;this part to draw a out circle of the disk 
	diskl=1.0632139225171164
	delta=360*!dtor/200.	
	aaa=fltarr(200) & disklx=aaa & diskly=aaa & disklz=aaa
	thetax=aaa & j=indgen(211)
	thetax=j*delta
	disklx=diskl*cos(thetax)
	diskly=diskl*sin(thetax)
	
	;plots,disklx,diskly,color=1b,thick=8b
	 ;plots,[0,0],[0,20],color=1b,thick=4b
;print,thetax,j,delta


for i=nmstart,nm-1 do begin
;print,nn,i

x1=xcoord[i*ln+lbb:(i+1)*ln-lee]*cp-ycoord[i*ln+lbb:(i+1)*ln-lee]*sp
y1=-xcoord[i*ln+lbb:(i+1)*ln-lee]*sp+ycoord[i*ln+lbb:(i+1)*ln-lee]*cp
z1=zcoord[i*ln+lbb:(i+1)*ln-lee]

;x3=xcoord3[i*ln+lbb:(i+1)*ln-lee]*cp-ycoord3[i*ln+lbb:(i+1)*ln-lee]*sp
;y3=-xcoord3[i*ln+lbb:(i+1)*ln-lee]*sp+ycoord3[i*ln+lbb:(i+1)*ln-lee]*cp
;z3=zcoord3[i*ln+lbb:(i+1)*ln-lee]

wh1=where(abs(x1) le 20 and abs(y1) le 20,nwh1)

;x1=x1(wh1)
;y1=y1(wh1)
;z1=z1(wh1)

x2=x1
;y2=y1*ct+z1*st
;z2=-y1*st+z1*ct
y2=y1*ct-z1*st
z2=y1*st+z1*ct

;x4=x3
;y4=y3*ct+z3*st
;z4=-y3*st+z3*ct
;y4=y3*ct-z3*st
;z4=y3*st+z3*ct

 ra = sqrt(1.0632139225171164D0^2+0.998D0^2)-0.01
print,ra
;  hide line segments that are behind disk
        ;wh1=where(x2^2+z2^2 le ra,nwh1)
        ;wh2=where((-y2*st+ct*z2 lt 0) and (((x2/diskl)^2+(z2/st/diskl)^2) gt 1.^2),nwh2)
        ;wh2=where( (y2 gt 0.) and (x2^2+z2^2 lt ra*ra),nwh2 )
        ;wh3=where( (y2 lt 0.) ,nwh3 )
        ;case 1 of
        ;  (nwh1 gt 0) and (nwh2 gt 0): wh=[wh1,wh2];wh=union(wh1,wh2)
        ;  (nwh1 gt 0) and (nwh2 eq 0): wh=wh1
       ;   (nwh1 eq 0) and (nwh2 gt 0): wh=wh2
        ;  (nwh1 eq 0) and (nwh2 eq 0): doline=0
        ;endcase
        ;if (nwh1+nwh2) gt 0 then doline=1 
        ;if doline then begin
print,0
           ; select the visible coordinates of the line
        ;  x3=x2
        ;  y3=y2
        ;  z3=z2
        ;  x3(wh2)=sqrt(ra-100.);ra*x3(wh2)/sqrt(x3(wh2)^2+z3(wh2)^2)
        ;  ;y3(wh2)=ra*x3(wh2)/sqrt(x3(wh2)^2+z3(wh2)^2)
        ;  z3(wh2)=sqrt(ra-100.);ra*z3(wh2)/sqrt(x3(wh2)^2+z3(wh2)^2)

	  ;x2=-1*x1
        ;  col=249b
          oplot,x2,z2,col=1b,thick=1.5
          ;oplot,x4,z4,col=col,thick=1.2,psym=3
	;endif 
endfor
            ;******************************************************** 
		colors=1
		chsize=0.7
		chth=2
		tickth=2
		alp1=-length
		alp2=length
		beta1=-length
		beta2=length
		xminor =1 
                yminor =1
                xticks =10
                yticks =10
                xtickslen=0.01
                ytickslen=0.01
		axis,xaxis=0,xticks=xticks,xminor=xminor,xrange=[alp1,alp2],xstyle=1,charsize=chsize,charthick=chth,font=-1,$
		xthick=tickth,xtitle=textoidl('X [GM/c^2]'),color=colors,xticklen=xtickslen;,xtickname=replicate(' ',11);,$
			;xtickname=['-4','-2','0','2','4','6','8'],

		axis,xaxis=1,xticks=xticks,xminor=xminor,xrange=[alp1,alp2],xstyle=1,xtickname=replicate(' ',11),font=-1,$
		charsize=chsize,charthick=chth,xthick=tickth,color=colors,xticklen=xtickslen

		axis,yaxis=0,yticks=yticks,yminor=yminor,yrange=[beta1,beta2],ystyle=1,charsize=chsize,charthick=chth,font=-1,$
		ythick=tickth,ytitle=textoidl('Y [GM/c^2]'),color=colors,yticklen=ytickslen;ytickname=replicate(' ',11),$
			;ytickname=['-6','-4','-2','0','2','4','6'],
		;,ytickv=[0,0.5,1]

	 	axis,yaxis=1,yticks=yticks,yminor=yminor,yrange=[beta1,beta2],ystyle=1,ytickname=replicate(' ',11),font=-1,$
		charsize=chsize,charthick=chth,ythick=tickth,color=colors,yticklen=ytickslen;
            ;*******************************************************

  If(theta0 eq 0.)then begin
	;oplot,[-11,-1.23],[0,0],col=250,thick=thick_CG 
	;oplot,[1.23,11],[0,0],col=250,thick=thick_CG 
  ENdif
  ;plots,disklx,diskly,color=249b,thick=2b

 device,/close 
 
 set_plot,oldn
 
 end
