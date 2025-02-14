begin

day="20230627"
files:=systemfunc("ls ../output/RTTOV_TELSEM_CLOUDY_SIMU_Swath_FY3G_MWRI-_ORBA_L1_"+day+"_*_7000M_V0.HDF")

print(files)

; exit
fios=addfiles(files,"r")
ListSetType(fios, "cat")

Longitude_2d:=fios[:]->Longitude
Latitude_2d:=fios[:]->Latitude
Longitude_2d@_FillValue=-999.9
Latitude_2d@_FillValue=-999.9

Longitude_2d@units="degree_east"
Latitude_2d@units="degree_north"

lat2d = Latitude_2d
lon2d = Longitude_2d

tbs_raw:=fios[:]->Emissivity_retrieval; 	Emissivity_retrieval	Geo2D
tbs_simu:=fios[:]->Emissivity_simu; 	Emissivity_retrieval	Geo2D
tbs := tbs_raw
; printVarSummary(tbs)
tbs@_FillValue=-999.9
; tbs@units="K"
; tbs@description="Description :	10-89GHz Earth observation brightness temperature datasets (10V,10H,18V,18H,23V,23H,36V,36H,89V,89H)"
; tbs=where(tbs.lt.10.or.tbs.gt.350,tbs@_FillValue,tbs)

; printMinMax(tbs,0)


tbs@lon2d=lon2d
tbs@lat2d=lat2d
tbs_simu@lon2d=lon2d
tbs_simu@lat2d=lat2d
; printVarSummary(tbs)

; do fr=0,9

wks=gsn_open_wks("png", "telsem_Emissivity_"+day) 
gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")
; gsn_define_colormap(wks, "NCV_Bright")
res=True
res@gsnDraw=False
res@gsnFrame=False
res@cnFillOn=True
res@cnFillMode="CellFill"
res@cnInfoLabelOn=False
res@cnLinesOn=False
res@cnLineLabelsOn=False

res@lbLabelBarOn=True
res@lbOrientation="Vertical"
res@pmLabelBarHeightF=0.25
res@pmLabelBarWidthF=0.06
res@lbBoxEndCapStyle="TriangleBothEnds"


res@mpFillOn=True

; res@gsnAddCyclic=True
res@gsnAddCyclic=False
res@mpMinLatF=-60
res@mpMaxLatF=60
res@mpMaxLonF=180
res@mpMinLonF=-180

plot=new(6,graphic)
res@cnLevelSelectionMode="ExplicitLevels"
res@cnLevels=(/0.4,0.5,0.55,0.6,0.8,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99/)
res@gsnLeftString="Retrieved_Emissivity_ 10V"
; res@gsnLeftString="Retrieved_Emissivity_Chn: "+fr
; plot=gsn_csm_contour_map(wks, tbs(:,:,fr), res)
plot(0)=gsn_csm_contour_map(wks, tbs(:,:,0), res)
; res@cnLevels:=(/170,175,180,190,200,210,220,230,240,250,260,265,270,275,280,285,290,295,300,305,310,320/)
res@gsnLeftString="Retrieved_Emissivity_ 19V"
plot(1)=gsn_csm_contour_map(wks, tbs(:,:,2), res)
; res@cnLevels:=(/170,175,180,190,200,210,220,230,240,250,260,265,270,275,280,285,290,295,300,305,310,320/)
res@gsnLeftString="Retrieved_Emissivity_ 37V"
plot(2)=gsn_csm_contour_map(wks, tbs(:,:,6), res)

res@gsnLeftString="Modeled & Telsem Emissivity_10V"
plot(3)=gsn_csm_contour_map(wks, tbs_simu(:,:,0), res)
; res@cnLevels:=(/170,175,180,190,200,210,220,230,240,250,260,265,270,275,280,285,290,295,300,305,310,320/)
res@gsnLeftString="Modeled & Telsem Emissivity_19V"
plot(4)=gsn_csm_contour_map(wks, tbs_simu(:,:,2), res)
; res@cnLevels:=(/170,175,180,190,200,210,220,230,240,250,260,265,270,275,280,285,290,295,300,305,310,320/)
res@gsnLeftString="Modeled & Telsem Emissivity_37V"
plot(5)=gsn_csm_contour_map(wks, tbs_simu(:,:,6), res)


resp=True
; resp@gsnPanelLabelBar=True
; resp@lbBoxEndCapStyle="TriangleBothEnds"

gsn_panel(wks,plot(3:),(/3,1/),resp)


; end do
end 