rm *.ps
ps=Basin_tseries.ps
gmt set PS_MEDIA=Custom_18.3cx8.5c
gmt set PS_PAGE_ORIENTATION=portrait
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_FRAME_PEN=0.4p,40/40/40
gmt set MAP_FRAME_TYPE=plain
gmt set MAP_LABEL_OFFSET=1.7p
gmt set MAP_TICK_LENGTH_PRIMARY=0.03c
gmt set MAP_TICK_PEN=thinnest,40,40,40
gmt set MAP_GRID_PEN_PRIMARY=thinnest,lightgrey
gmt set MAP_ANNOT_OFFSET_PRIMARY=1.5p
gmt set PS_LINE_CAP=butt
gmt set PS_LINE_JOIN=round

gmt set FONT_ANNOT_PRIMARY             = 7p,Helvetica,black
gmt set FONT_ANNOT_SECONDARY           = 7p,Helvetica,black
gmt set FONT_LABEL                     = 7p,Helvetica,black
gmt set FONT_LOGO                      = 7p,Helvetica,black
gmt set FONT_TITLE                     = 7p,Helvetica,black

Jbasin=X5.6c/3.5c

xright=5.8c
yright=0.0c

xleft=-11.6c
yleft=-3.9c

color0=20/20/20 # Budget
color1=#0072B2 # Observed
color2=#ff7f0e # Steric
color3=#d62728 # Mass
color4=#56B4E9 # Altimetry
color5=#9467bd # Glaciers
color6=#8c564b # GrIS
color7=#e377c2 # AIS
color8=#bcbd22 # TWS 
color5=#2ca02c #TWS natural
color10=#17becf # Dams
color11=#7f7f7f # GWD


rbasin=1900/2018/-360/80

rtrend=0/14/-0.5/5.5
Jtrend=X2.0c/1.5c

xstep=25g12.5
ystep=100g50

basin=0
gmt psbasemap -K -R$rbasin -J$Jbasin -X1.0c -Y4.7c -Bx$xstep -By$ystep+l'Height (mm)' -BWesn+t"Subpolar North Atlantic (5.8%)" > $ps
gmt psxy  -R -J -L -t60 -N -G$color5 -O -K "gia_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color3 -O -K "grd_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color2 -O -K "steric_${basin}_tseries_ci.txt" >> $ps
gmt psxy  -R -J -L -t60 -N -G$color0 -O -K "budget_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color4 -O -K "alt_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color1 -O -K "obs_${basin}_tseries_ci.txt"    >> $ps
gmt set PS_LINE_CAP=round
gmt psxy  -R -J -W1.5p,$color5,1_2:1 -O -K "gia_${basin}_tseries_m.txt"    >> $ps
gmt set PS_LINE_CAP=butt
gmt psxy  -R -J -W1.5p,$color3 -O -K "grd_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color2 -O -K "steric_${basin}_tseries_m.txt" >> $ps
gmt psxy  -R -J -W1.5p,$color0 -O -K "budget_${basin}_tseries_m.txt"  >> $ps
gmt psxy  -R -J -W1.5p,$color1 -O -K "obs_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color4 -O -K "alt_${basin}_tseries_m.txt"    >> $ps
echo "1902 -350 a" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLB -O -K >> $ps

echo 1
basin=2
gmt psbasemap -O -K -R$rbasin -J$Jbasin -X$xright -Y$yright -Bx$xstep -By$ystep -Bwesn+t"Subtropical North Atlantic (5.0%)" >> $ps
gmt psxy  -R -J -L -t60 -N -G$color5 -O -K "gia_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color3 -O -K "grd_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color2 -O -K "steric_${basin}_tseries_ci.txt" >> $ps
gmt psxy  -R -J -L -t60 -N -G$color0 -O -K "budget_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color4 -O -K "alt_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color1 -O -K "obs_${basin}_tseries_ci.txt"    >> $ps
gmt set PS_LINE_CAP=round
gmt psxy  -R -J -W1.5p,$color5,1_2:1 -O -K "gia_${basin}_tseries_m.txt"    >> $ps
gmt set PS_LINE_CAP=butt
gmt psxy  -R -J -W1.5p,$color3 -O -K "grd_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color2 -O -K "steric_${basin}_tseries_m.txt" >> $ps
gmt psxy  -R -J -W1.5p,$color0 -O -K "budget_${basin}_tseries_m.txt"  >> $ps
gmt psxy  -R -J -W1.5p,$color1 -O -K "obs_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color4 -O -K "alt_${basin}_tseries_m.txt"    >> $ps
echo "1902 -350 b" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLB -O -K >> $ps

basin=4
gmt psbasemap -O -K -R$rbasin -J$Jbasin -X$xright -Y$yright  -Bx$xstep -By$ystep -Bwesn+t"South Atlantic (21.8%) " >> $ps
gmt psxy  -R -J -L -t60 -N -G$color5 -O -K "gia_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color3 -O -K "grd_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color2 -O -K "steric_${basin}_tseries_ci.txt" >> $ps
gmt psxy  -R -J -L -t60 -N -G$color0 -O -K "budget_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color4 -O -K "alt_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color1 -O -K "obs_${basin}_tseries_ci.txt"    >> $ps
gmt set PS_LINE_CAP=round
gmt psxy  -R -J -W1.5p,$color5,1_2:1 -O -K "gia_${basin}_tseries_m.txt"    >> $ps
gmt set PS_LINE_CAP=butt
gmt psxy  -R -J -W1.5p,$color3 -O -K "grd_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color2 -O -K "steric_${basin}_tseries_m.txt" >> $ps
gmt psxy  -R -J -W1.5p,$color0 -O -K "budget_${basin}_tseries_m.txt"  >> $ps
gmt psxy  -R -J -W1.5p,$color1 -O -K "obs_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color4 -O -K "alt_${basin}_tseries_m.txt"    >> $ps
echo "1902 -350 c" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLB -O -K >> $ps

basin=1
gmt psbasemap -O -K -R$rbasin -J$Jbasin -X$xleft -Y$yleft -Bx$xstep -By$ystep+l'Height (mm)' -BWeSn+t" Indian Ocean - South Pacific (37.6%)" >> $ps
gmt psxy  -R -J -L -t60 -N -G$color5 -O -K "gia_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color3 -O -K "grd_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color2 -O -K "steric_${basin}_tseries_ci.txt" >> $ps
gmt psxy  -R -J -L -t60 -N -G$color0 -O -K "budget_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color4 -O -K "alt_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color1 -O -K "obs_${basin}_tseries_ci.txt"    >> $ps
gmt set PS_LINE_CAP=round
gmt psxy  -R -J -W1.5p,$color5,1_2:1 -O -K "gia_${basin}_tseries_m.txt"    >> $ps
gmt set PS_LINE_CAP=butt
gmt psxy  -R -J -W1.5p,$color3 -O -K "grd_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color2 -O -K "steric_${basin}_tseries_m.txt" >> $ps
gmt psxy  -R -J -W1.5p,$color0 -O -K "budget_${basin}_tseries_m.txt"  >> $ps
gmt psxy  -R -J -W1.5p,$color1 -O -K "obs_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color4 -O -K "alt_${basin}_tseries_m.txt"    >> $ps
echo "1902 -350 d" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLB -O -K >> $ps

#
basin=3
gmt psbasemap -O -K -R$rbasin -J$Jbasin -X$xright -Y$yright -Bx$xstep -By$ystep -BweSn+t" East Pacific (15.9%)" >> $ps
gmt psxy  -R -J -L -t60 -N -G$color5 -O -K "gia_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color3 -O -K "grd_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color2 -O -K "steric_${basin}_tseries_ci.txt" >> $ps
gmt psxy  -R -J -L -t60 -N -G$color0 -O -K "budget_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color4 -O -K "alt_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color1 -O -K "obs_${basin}_tseries_ci.txt"    >> $ps
gmt set PS_LINE_CAP=round
gmt psxy  -R -J -W1.5p,$color5,1_2:1 -O -K "gia_${basin}_tseries_m.txt"    >> $ps
gmt set PS_LINE_CAP=butt
gmt psxy  -R -J -W1.5p,$color3 -O -K "grd_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color2 -O -K "steric_${basin}_tseries_m.txt" >> $ps
gmt psxy  -R -J -W1.5p,$color0 -O -K "budget_${basin}_tseries_m.txt"  >> $ps
gmt psxy  -R -J -W1.5p,$color1 -O -K "obs_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color4 -O -K "alt_${basin}_tseries_m.txt"    >> $ps
echo "1902 -350 e" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLB -O -K >> $ps

basin=5
gmt psbasemap -O -K -R$rbasin -J$Jbasin -X$xright -Y$yright -Bx$xstep -By$ystep -BweSn+t"Northwest Pacific (13.9%)"      >> $ps
gmt psxy  -R -J -L -t60 -N -G$color5 -O -K "gia_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color3 -O -K "grd_${basin}_tseries_ci.txt"    >> $ps
gmt psxy  -R -J -L -t60 -N -G$color2 -O -K "steric_${basin}_tseries_ci.txt" >> $ps
gmt psxy  -R -J -L -t60 -N -G$color0 -O -K "budget_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color4 -O -K "alt_${basin}_tseries_ci.txt"  >> $ps
gmt psxy  -R -J -L -t60 -N -G$color1 -O -K "obs_${basin}_tseries_ci.txt"    >> $ps
gmt set PS_LINE_CAP=round
gmt psxy  -R -J -W1.5p,$color5,1_2:1 -O -K "gia_${basin}_tseries_m.txt"    >> $ps
gmt set PS_LINE_CAP=butt
gmt psxy  -R -J -W1.5p,$color3 -O -K "grd_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color2 -O -K "steric_${basin}_tseries_m.txt" >> $ps
gmt psxy  -R -J -W1.5p,$color0 -O -K "budget_${basin}_tseries_m.txt"  >> $ps
gmt psxy  -R -J -W1.5p,$color1 -O -K "obs_${basin}_tseries_m.txt"    >> $ps
gmt psxy  -R -J -W1.5p,$color4 -O -K "alt_${basin}_tseries_m.txt"    >> $ps
echo "1902 -350 f" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLB -O -K >> $ps

gmt psbasemap -O -K -R0/0.81/0/1 -JX16.2c/0.3c -X-11.1c -Y-0.75c -T-100/100/0.1 >> $ps
echo -e "0.000 0.5 \n 0.030 0.5" | psxy -R -J -O -K -W6p,$color1 -t70 >> $ps
echo -e "0.000 0.5 \n 0.030 0.5" | psxy -R -J -O -K -W1.5p,$color1 >> $ps
echo "0.033 0.5 Observed sea level" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.160 0.5 \n 0.190 0.5" | psxy -R -J -O -K -W6p,$color4 -t70 >> $ps
echo -e "0.160 0.5 \n 0.190 0.5" | psxy -R -J -O -K -W1.5p,$color4 >> $ps
echo "0.193 0.5 Altimetry" | gmt pstext -R -J -F+f7,Helvetica+jLM -O -K >> $ps
echo -e "0.260 0.5 \n 0.290 0.5" | psxy -R -J -O -K -W6p,$color0 -t70 >> $ps
echo -e "0.260 0.5 \n 0.290 0.5" | psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "0.293 0.5 Sum of contributors" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.420 0.5 \n 0.450 0.5" | psxy -R -J -O -K -W6p,$color2 -t70 >> $ps
echo -e "0.420 0.5 \n 0.450 0.5" | psxy -R -J -O -K -W1.5p,$color2 >> $ps
echo "0.453 0.5 Steric" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.500 0.5 \n 0.530 0.5" | psxy -R -J -O -K -W6p,$color3 -t70 >> $ps
echo -e "0.500 0.5 \n 0.530 0.5" | psxy -R -J -O -K -W1.5p,$color3 >> $ps
echo "0.533 0.5 Ocean mass" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.620 0.5 \n 0.650 0.5" | psxy -R -J -O -K -W6p,$color5 -t70 >> $ps
gmt set PS_LINE_CAP=round
echo -e "0.620 0.5 \n 0.650 0.5" | psxy -R -J -O -K -W1.5p,$color5,1_2:0.5 >> $ps
echo "0.653 0.5 Glacial Isostatic Adjustment" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps



gmt psconvert $ps -Tf
