begin

mersi=addfile("Cloud_Fraction_5km_MWRIORB_20230801_0055.HDF","r")
cfr=mersi->Cloud_Fraction_5km
lon=mersi->Longitude
lat=mersi->Latitude


lon@units="degree_east"
lat@units="degree_north"

cfr!0="lat"
cfr!1="lon"
cfr&lon=lon
cfr&lat=lat
cfr@units="%"
cfr@_FillValue=-999.9
cfr=where(cfr.lt.0,cfr@_FillValue,cfr)


mwri=addfile("/data04/0/MWRI/FY3G_IOT/ASCEND/20230801/FY3G_MWRI-_ORBA_L1_20230801_0055_7000M_V0.HDF","r")
TBs=mwri->/S1/Data/EARTH_OBSERVE_BT_10_to_89GHz*0.01+327.68
Lon2d=mwri->/S1/Geolocation/Longitude
Lat2d=mwri->/S1/Geolocation/Latitude

Lon2d@units="degree_east"
Lat2d@units="degree_north"
TBs@lon2d=Lon2d
TBs@lat2d=Lat2d
TBs@units="K"
printVarSummary(TBs)


setvalues NhlGetWorkspaceObjectId()
"wsMaximumSize": 2000000000
end setvalues

wks=gsn_open_wks("newps","MWRIORB_20230801_0055-zoom")
gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")
plot=new(2,graphic)

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

res@mpFillOn=True
res@mpOutlineBoundarySets="National"
res@mpDataSetName="Earth..4"

res@mpOceanFillColor="lightblue"


res@gsnAddCyclic=False
res@mpMinLatF=20
res@mpMaxLatF=60
res@mpMinLonF=70
res@mpMaxLonF=120
res@cnLevelSelectionMode="ExplicitLevels"
res@cnLevels=(/5,10,15,20,30,40,50,60,70,80,90/)
res@cnMissingValFillColor="lightgray"
; res@lbLabelAlignment="BoxCenters"

res@gsnLeftString="MERSI CFR 0.05 Grid"
; res@gsnLeftStringFontHeightF=0.03
plot(0)=gsn_csm_contour_map(wks, cfr, res)

respl=True
respl@gsLineColor="Green"
respl@gsLineThicknessF=2

dims=dimsizes(Lon2d)
plR=gsn_add_polyline(wks, plot(0), Lon2d(:,0), Lat2d(:,0), respl)
plL=gsn_add_polyline(wks, plot(0), Lon2d(:,dims(1)-1), Lat2d(:,dims(1)-1), respl)

plS=gsn_add_polyline(wks, plot(0), Lon2d(0,:), Lat2d(0,:), respl)
plE=gsn_add_polyline(wks, plot(0), Lon2d(dims(0)-1,:), Lat2d(dims(0)-1,:), respl)



res@cnLevels:=(/150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310/)
res@gsnLeftString="MWRI TB10V "
plot(1)=gsn_csm_contour_map(wks, TBs(:,:,0), res)


resP=True
gsn_panel(wks, plot, (/2,1/), resP)

end
