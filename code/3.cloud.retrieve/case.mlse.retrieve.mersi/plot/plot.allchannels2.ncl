begin

files=systemfunc("ls ../output/202308??/Emissivity_Clear_FY3G_MWRIA_2023080[1]_????_V0.HDF")

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
wks:=gsn_open_wks("png",filename+"channel-5")
gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")
; plot:=new(3,graphic)

res=True
res@gsnDraw=False
res@gsnFrame=False

res@cnFillOn=True
res@cnFillMode="CellFill"
res@cnLinesOn=False
res@cnLineLabelsOn=False
res@cnInfoLabelOn=False
res@lbLabelBarOn=False
res@lbBoxEndCapStyle="TriangleBothEnds"
; res@lbOrientation="vertical"
; res@pmLabelBarWidthF=0.07

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

res@cnMissingValFillColor="lightgray"
; res@lbLabelAlignment="BoxCenters"

plot=new(4,graphic)
res@cnLevelSelectionMode="ExplicitLevels"

res@gsnLeftStringFontHeightF=0.025
res@lbLabelFontHeightF=0.02
res@gsnLeftString="Retrieved_Emissivity_ 23V"
res@cnLevels:=(/0.85,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)
plot(0)=gsn_csm_contour_map(wks,  Emissivity(:,:,4), res)

res@cnLevels:=(/0.4,0.5,0.55,0.6,0.8,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)
res@gsnLeftString="Retrieved_Emissivity_ 23H"
plot(1)=gsn_csm_contour_map(wks,  Emissivity(:,:,5), res)


res@lbLabelBarOn=True
res@cnLevels:=(/0.85,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)
res@gsnLeftString="Retrieved_Emissivity_ 89V"
plot(2)=gsn_csm_contour_map(wks,  Emissivity(:,:,8), res)

res@cnLevels:=(/0.4,0.5,0.55,0.6,0.8,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)

res@gsnLeftString="Retrieved_Emissivity_ 89H"
plot(3)=gsn_csm_contour_map(wks,  Emissivity(:,:,9), res)



resp=True
; resp@gsnPanelLabelBar=True
; resp@lbBoxEndCapStyle="TriangleBothEnds"
resp@gsnPanelMainString= date+" "+time
gsn_panel(wks,plot,(/2,2/),resp)
end do
end
