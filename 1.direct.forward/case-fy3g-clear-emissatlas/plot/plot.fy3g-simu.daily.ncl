begin

day="20230627"
files:=systemfunc("ls ../output/RTTOV_CLEAR_SIMU_Swath_FY3G_MWRI-_ORBA_L1_"+day+"_*_7000M_V0.HDF")

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

tbs_raw:=fios[:]->TB_FY3G_simu
tbs := tbs_raw
; printVarSummary(tbs)
tbs@_FillValue=-999.9
tbs@units="K"
tbs@description="Description :	10-89GHz Earth observation brightness temperature datasets (10V,10H,18V,18H,23V,23H,36V,36H,89V,89H)"
tbs=where(tbs.lt.10.or.tbs.gt.350,tbs@_FillValue,tbs)

; printMinMax(tbs,0)


tbs@lon2d=lon2d
tbs@lat2d=lat2d
printVarSummary(tbs)


wks=gsn_open_wks("png", "FY3G_TBs_RTTOV-TELSEM2_"+day) 
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
; res@mpMaxLonF=145
; res@mpMinLonF=70

plot=new(3,graphic)
res@cnLevelSelectionMode="ExplicitLevels"
res@cnLevels=(/150,155,160,162,165,168,170,175,180,240,250,260,265,270,275,280,285,290,295,300,305,310,320/)
res@gsnLeftString="RTTOV Cloudy TB 10V (TELSEM2-LSE + Fastem-SSE)"
plot(0)=gsn_csm_contour_map(wks, tbs(:,:,0), res)
res@cnLevels:=(/170,175,180,190,200,210,220,230,240,250,260,265,270,275,280,285,290,295,300,305,310,320/)
res@gsnLeftString="RTTOV Cloudy TB 19V (TELSEM2-LSE + Fastem-SSE)"
plot(1)=gsn_csm_contour_map(wks, tbs(:,:,2), res)
res@cnLevels:=(/170,175,180,190,200,210,220,230,240,250,260,265,270,275,280,285,290,295,300,305,310,320/)
res@gsnLeftString="RTTOV Cloudy TB 37V (TELSEM2-LSE + Fastem-SSE)"
plot(2)=gsn_csm_contour_map(wks, tbs(:,:,6), res)


resp=True
; resp@gsnPanelLabelBar=True
; resp@lbBoxEndCapStyle="TriangleBothEnds"

gsn_panel(wks,plot,(/3,1/),resp)
end 