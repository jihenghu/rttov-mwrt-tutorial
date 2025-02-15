begin

files=systemfunc("ls ../output/202308??/Emissivity_Clear_FY3G_MWRIA_2023080[2345]_????_V0.HDF")

print(""+files)

do iF=0,dimsizes(files)-1

emissio=addfile(files(iF),"r")
Emissivity:=emissio->Emissivity
TELSEM2:=emissio->TELSEM2_Emissivity
cfr:=emissio->MERSI_Cloud_Fraction
lsm:=emissio->Land_Sea_Mask

lat2d:=emissio->Latitude
lon2d:=emissio->Longitude



lon2d@_FillValue=-999.9
lat2d@_FillValue=-999.9
lon2d=where(lon2d.gt.-181.and.lon2d.lt.181,lon2d,lon2d@_FillValue )
lat2d=where(lat2d.gt.-90.and.lat2d.lt.90,lat2d,lat2d@_FillValue )

Emissivity@lon2d=lon2d
Emissivity@lat2d=lat2d
TELSEM2@lon2d=lon2d
TELSEM2@lat2d=lat2d

cfr@lon2d=lon2d
cfr@lat2d=lat2d
Emissivity@_FillValue=-999.9


Emissivity=where(Emissivity.lt.0.or.Emissivity.ge.1,Emissivity@_FillValue,Emissivity)
TELSEM2@_FillValue=-999.9
printVarSummary(Emissivity)
printVarSummary(lsm)

; lsm3d:=conform(Emissivity,lsm,(/0,1/))
; printVarSummary(lsm3d)

Emissivity=where(cfr.gt.20,Emissivity@_FillValue,Emissivity)
; Emissivity=where(lsm3d.ne.1,Emissivity@_FillValue,Emissivity)

TELSEM2=where(ismissing(Emissivity), TELSEM2@_FillValue ,TELSEM2)

filename=systemfunc("echo "+files(iF)+" | cut -d '/' -f4 | cut -d '.' -f1")
time=systemfunc("echo "+files(iF)+" | cut -d '/' -f4 | cut -d '_' -f6")
date=systemfunc("echo "+files(iF)+" | cut -d '/' -f3")

; print(filename)
wks:=gsn_open_wks("png",filename)
gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")
plot:=new(3,graphic)

res=True
res@gsnDraw=False
res@gsnFrame=False

res@cnFillOn=True
res@cnFillMode="CellFill"
res@cnLinesOn=False
res@cnLineLabelsOn=False
res@cnInfoLabelOn=False
res@lbLabelBarOn=True
res@lbBoxEndCapStyle="TriangleBothEnds"
res@lbOrientation="vertical"
res@pmLabelBarWidthF=0.07

res@gsnSpreadColorStart=20

res@mpFillOn=True
res@mpOutlineBoundarySets="National"
res@mpDataSetName="Earth..4"

res@mpOceanFillColor="transparent"
res@mpLandFillColor="transparent"

res@gsnAddCyclic=True
res@mpMinLatF=min(lat2d)-10
res@mpMaxLatF=max(lat2d)+10

res@cnLevelSelectionMode="ExplicitLevels"
res@cnLevels:=(/0.85,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)
res@cnMissingValFillColor="lightgray"
; res@lbLabelAlignment="BoxCenters"

res@gsnLeftString="(a) Emssivity 10V"
; res@gsnLeftStringFontHeightF=0.03
plot(0)=gsn_csm_contour_map_ce(wks, Emissivity(:,:,0), res)

res@gsnLeftString="(b) TELSEM2 Monthly Emssivity 10V"
; res@gsnLeftStringFontHeightF=0.03
plot(1)=gsn_csm_contour_map_ce(wks, TELSEM2(:,:,0), res)


res@cnLevels:=(/5,10,15,20,30,40,50,60,70,80,90/)
res@gsnLeftString="(c) Cloud Fraction"
plot(2)=gsn_csm_contour_map_ce(wks, cfr(:,:,0), res)


resP=True
resP@gsnFrame    = False  
resP@gsnPanelMainString    = "FY3G MWRI Orbit: "+date+" "+time
gsn_panel(wks, plot, (/3,1/), resP)

; Set up some histogram resources.
;
histres = True
histres@gsnFrame=False

histres@vpXF      = 0.14
histres@vpYF      = 0.84
histres@vpWidthF  = 0.18
histres@vpHeightF = 0.15

; Turn off X and Y tick marks and X/Y borders.
;
;   histres@tmXBOn       = False
  histres@tmYROn       = False
  histres@tmYLOn       = False
;   histres@tmXBBorderOn = False
  histres@tmXTBorderOn = False
  histres@tmYLBorderOn = False
  histres@tmYRBorderOn = False
;
; Turn off X axis label.
;
;   histres@tiYAxisOn = False
;   histres@gsnLeftString=" Frequency"

  ; opt=True
  ; opt@bin_bound_min =0.7
  ; opt@bin_bound_max =1
  ; pdf=pdfx(Emissivity(:,:,0),10,opt)
histres@tmYLLabelFontHeightF = 0.01
histres@tmXBLabelFontHeightF = 0.01

histres@tmXBLabelStride = 2


histres@tiYAxisOn = False
histres@tiXAxisOn =  False


histres@gsnHistogramBinIntervals := (/0.5,0.7,0.85,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.00/)
; Draw histogram.
emiss:=ndtooned(Emissivity(:,:,0))
emissx:=emiss(ind(emiss.gt.0.7.and.emiss.lt.1))
histogram1 = gsn_histogram(wks,emissx,histres)



histres@vpXF      = 0.14
histres@vpYF      = 0.53
emiss:=ndtooned(TELSEM2(:,:,0))
emissx:=emiss(ind(emiss.gt.0.7.and.emiss.lt.1))
histogram2 = gsn_histogram(wks,emissx,histres)

frame(wks)

end do
end
