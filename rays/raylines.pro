pro  raylines

theta0=90.
phi0=0
WINSIZE=20
bkg='wwhite'
backg='ny'
nm=21*21uL
ln=201ul
nmstart=0ul
lee=1ul
lbb=0ul
;@nlff_block

if not keyword_set(bkg) then bkg='white' else bkg=bkg(0)

 Openr,xcoor,'rayx2.txt',/Get_Lun
  Point_lun,xcoor,0
  xcoord=fltarr(ln,nm)
  ReadF,xcoor,xcoord
  
   Openr,ycoor,'rayy2.txt',/Get_Lun
  Point_lun,ycoor,0
  ycoord=fltarr(ln,nm)
  ReadF,ycoor,ycoord
  
   Openr,zcoor,'rayz2.txt',/Get_Lun
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
  
  ;Openr,lunAo102,'ray_diskg.txt',/Get_Lun
  ;Point_lun,lunAo102,0
  disk2=fltarr(401,401)
  ;ReadF,lunAo102,disk2
  
  bz=congrid(-disk2*0.1,100,100)
  
  free_lun,xcoor,ycoor,zcoor;,lunAo102
;print,xcoord(*,73),ycoord(*,73),'zzzzzzzzzzz=',zcoord(*,73)


  length=5
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



loadct,33  ;  loadct,3 also looks nice too

A=[xmin,ymin,zmin]
B=[xmin,ymax,zmin]
C=[xmax,ymax,zmin]
D=[xmax,ymin,zmin]
DZ=[0,0,zmax-zmin]
E=A+DZ
F=B+DZ
G=C+DZ
H=D+DZ

A=[A[0]*cp+A[1]*sp,-A[0]*sp+A[1]*cp,A[2]]
B=[B[0]*cp+B[1]*sp,-B[0]*sp+B[1]*cp,B[2]]
C=[C[0]*cp+C[1]*sp,-C[0]*sp+C[1]*cp,C[2]]
D=[D[0]*cp+D[1]*sp,-D[0]*sp+D[1]*cp,D[2]]

E=[E[0]*cp+E[1]*sp,-E[0]*sp+E[1]*cp,E[2]]
F=[F[0]*cp+F[1]*sp,-F[0]*sp+F[1]*cp,F[2]]
G=[G[0]*cp+G[1]*sp,-G[0]*sp+G[1]*cp,G[2]]
H=[H[0]*cp+H[1]*sp,-H[0]*sp+H[1]*cp,H[2]]

A=[A[0],A[1]*ct+A[2]*st,-A[1]*st+A[2]*ct]
B=[B[0],B[1]*ct+B[2]*st,-B[1]*st+B[2]*ct]
C=[C[0],C[1]*ct+C[2]*st,-C[1]*st+C[2]*ct]
D=[D[0],D[1]*ct+D[2]*st,-D[1]*st+D[2]*ct]

E=[E[0],E[1]*ct+E[2]*st,-E[1]*st+E[2]*ct]
F=[F[0],F[1]*ct+F[2]*st,-F[1]*st+F[2]*ct]
G=[G[0],G[1]*ct+G[2]*st,-G[1]*st+G[2]*ct]
H=[H[0],H[1]*ct+H[2]*st,-H[1]*st+H[2]*ct]

x11=min([A[0],B[0],C[0],D[0],E[0],F[0],G[0],H[0]])
x22=max([A[0],B[0],C[0],D[0],E[0],F[0],G[0],H[0]])
z11=min([A[2],B[2],C[2],D[2],E[2],F[2],G[2],H[2]])
z22=max([A[2],B[2],C[2],D[2],E[2],F[2],G[2],H[2]])

rt=1.05
aga_flag=0

xs0=round(1.2*WINSIZE*(x22-x11)) & ys0=round(1.2*WINSIZE*(z22-z11))


BOX_XSIZE=round(WINSIZE*(x22-x11))*rt 
BOX_YSIZE=round(WINSIZE*(z22-z11))*rt
X0=(xs0-BOX_XSIZE)/2.
Y0=(ys0-BOX_YSIZE)/2.

l=16  & xxss=l & yyss=l*ys0/xs0

device,filename='raylinesp_2.ps',color=1,xsize=xxss,ysize=yyss,$
   xoff=(21.-xxss)/2. , yoff=(21-yyss)/2.,bits_per_pixel=8;decomposed=0

plot,[x11*rt,x22*rt],[z11*rt,z22*rt],$
pos=[(X0/xs0)*xxss*1000,(Y0/ys0)*yyss*1000,xxss*1000*(X0+BOX_XSIZE)/xs0,yyss*1000*(Y0+BOX_YSIZE)/ys0],$
    /noerase,/device,/nodata,$
    xsty=5,ysty=5,xr=[x11*rt,x22*rt],yr=[z11*rt,z22*rt]


XX=round(WINSIZE*(xmax-xmin))
YY=round(WINSIZE*(ymax-ymin))

AA=[WINSIZE*(A[0]-x11*rt)+X0,WINSIZE*(A[2]-z11*rt)+Y0]
BB=[WINSIZE*(B[0]-x11*rt)+X0,WINSIZE*(B[2]-z11*rt)+Y0]
CC=[WINSIZE*(C[0]-x11*rt)+X0,WINSIZE*(C[2]-z11*rt)+Y0]
DD=[WINSIZE*(D[0]-x11*rt)+X0,WINSIZE*(D[2]-z11*rt)+Y0]

EE=[WINSIZE*(E[0]-x11*rt)+X0,WINSIZE*(E[2]-z11*rt)+Y0]
FF=[WINSIZE*(F[0]-x11*rt)+X0,WINSIZE*(F[2]-z11*rt)+Y0]
GG=[WINSIZE*(G[0]-x11*rt)+X0,WINSIZE*(G[2]-z11*rt)+Y0]
HH=[WINSIZE*(H[0]-x11*rt)+X0,WINSIZE*(H[2]-z11*rt)+Y0]

If (bkg eq 'white') then begin
;^^^^^^^^^^^^white backgroud^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
im=bytscl(bz)
im=congrid(im,XX,YY,/interp)


xi=[0, 0,XX,XX]
yi=[YY,0,0, YY]


xo=[AA[0],BB[0],CC[0],DD[0]]
yo=[AA[1]+EE[1],BB[1]+FF[1],CC[1]+GG[1],DD[1]+HH[1]]*0.5
;yo=[AA[1],BB[1],CC[1],DD[1]]

if(theta ne 0)then begin
warped_image = WARP_TRI(xo, yo, xi, yi, im, $
  OUTPUT_SIZE=[xs0,ys0]);,  /QUINTIC);OUTPUT_SIZE=,

xc=0 & yc=0.
xs=1 & ys=1
warped_imgae=bytscl(warped_image)
If (backg eq 'y') then begin
   ;tvscl,warped_image,xc,yc,xsize=xs,ysize=ys,/norm
endif else begin
  ; tvscl,warped_image*0,xc,yc,xsize=xs,ysize=ys,/norm
endelse
endif

Endif else begin

;^^^^^^^^^^^^^^^^^^^^^^^^^^^black background^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ims=size(bz,/dim)
rhs = [0, 0, ims(0), 0, 0, ims(1)]
u2=[BB[0],CC[0],AA[0],DD[0]]
v2=[BB[1],CC[1],AA[1],DD[1]]

m = [[u2[0], v2[0], 0, 0, 1, 0], [0, 0, u2[0], v2[0], 0, 1], $
       [u2[1], v2[1], 0, 0, 1, 0], [0, 0, u2[1], v2[1], 0, 1], $
       [u2[2], v2[2], 0, 0, 1, 0], [0, 0, u2[2], v2[2], 0, 1]]

ludc, m, index, /double
  sol = lusol(m, index, rhs, /double)
 kx = fltarr(2, 2)
  ky = fltarr(2, 2)
  kx[0, 1] = sol[0]
  kx[1, 0] = sol[1]
  ky[0, 1] = sol[2]
  ky[1, 0] = sol[3]
  kx[0, 0] = sol[4]
  ky[0, 0] = sol[5]
miss=0

    imm = poly_2d(bytscl(bz), kx, ky, 2, xs0,ys0,$
                missing = miss) ;Warp it

    
xc=0 & yc=0
xs=1 & ys=1
imm=bytscl(imm)
;tvscl,imm,xc,yc,xsize=xs,ysize=ys,/norm
print,'we are here'
endelse
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


loadct,3
tvlct,re,gr,bl,/get
re(249:255)=[0b,255b,0b,0b,255b,0b,255b]
gr(249:255)=[0b,0b,255b,0b,255b,255b,255b]
bl(249:255)=[0b,0b,0b,255b,0b,255b,255b]
tvlct,re,gr,bl

if (phi ge 0) and (phi lt !pi/2.) then begin
   thick_AB=3
   thick_BC=3
   thick_DC=1
   thick_AD=1

   thick_AE=3  
   thick_BF=3
   thick_CG=3 
   thick_HD=1

   thick_EF=3
   thick_EH=3
   thick_FG=3
   thick_GH=3 
endif 
if (phi ge !pi/2.) and (phi lt !pi) then begin
   thick_AB=3
   thick_BC=1
   thick_DC=1
   thick_AD=3

   thick_AE=3  
   thick_BF=3
   thick_CG=1 
   thick_HD=3

   thick_EF=3
   thick_EH=3
   thick_FG=3
   thick_GH=3
endif

if (phi ge !pi) and (phi lt !pi*3/2) then begin
   thick_AB=1
   thick_BC=1
   thick_DC=3
   thick_AD=3

   thick_AE=3  
   thick_BF=1
   thick_CG=3 
   thick_HD=3

   thick_EF=3
   thick_EH=3
   thick_FG=3
   thick_GH=3
endif
if (phi ge !pi*3/2) and (phi lt !pi*2) then begin
   thick_AB=1
   thick_BC=3
   thick_DC=3
   thick_AD=1

   thick_AE=1  
   thick_BF=3
   thick_CG=3 
   thick_HD=3

   thick_EF=3
   thick_EH=3
   thick_FG=3
   thick_GH=3
endif

;this part to draw a out circle of the disk 
	diskl=1.35
	delta=360*!dtor/200.	
	aaa=fltarr(200) & disklx=aaa & diskly=aaa & disklz=aaa
	thetax=aaa & j=indgen(211)
	thetax=j*delta
	disklx=diskl*cos(thetax)
	diskly=diskl*sin(thetax)
	
	plots,disklx,diskly,color=1b,thick=8b
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

;  hide line segments that are behind disk
       ; wh1=where(z2^2+y2^2 le 20,nwh1)
       ; wh2=where((-y2*st+ct*z2 lt 0) and (((x2/diskl)^2+(z2/st/diskl)^2) gt 1.^2),nwh2)
       ; case 1 of
       ;   (nwh1 gt 0) and (nwh2 gt 0): wh=[wh1,wh2];wh=union(wh1,wh2)
       ;   (nwh1 gt 0) and (nwh2 eq 0): wh=wh1
       ;   (nwh1 eq 0) and (nwh2 gt 0): wh=wh2
       ;   (nwh1 eq 0) and (nwh2 eq 0): doline=0
       ; endcase
       ; if (nwh1+nwh2) gt 0 then doline=1

        ;if doline then begin

           ; select the visible coordinates of the line
         ; x2=x2(wh1)
         ; y2=y2(wh)
         ; z2=z2(wh1)

	  ;x2=-1*x1
          col=249b
          oplot,x2,z2,col=col,thick=1.2 
          ;oplot,x4,z4,col=col,thick=1.2,psym=3
	;endif 
endfor
            ;******************************************************** 
		colors=1
		chsize=0.7
		chth=2
		tickth=2
		alp1=-5
		alp2=5
		beta1=-5
		beta2=5
		xminor =2 
                yminor =2
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
	oplot,[-11,-1.23],[0,0],col=250,thick=thick_CG 
	oplot,[1.23,11],[0,0],col=250,thick=thick_CG 
  ENdif
  ;plots,disklx,diskly,color=249b,thick=2b
 device,/close 
 
 set_plot,oldn
 
 end
