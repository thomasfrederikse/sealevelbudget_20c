
rm *.ps
ps=Global_tseries_trends.ps
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_MEDIA=Custom_18.3cx8.8c
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_FRAME_PEN=0.4p,40/40/40
gmt set MAP_FRAME_TYPE=plain
gmt set MAP_LABEL_OFFSET=0.6p
gmt set MAP_TICK_LENGTH_PRIMARY=0.03c
gmt set MAP_TICK_PEN=thinnest,40,40,40
gmt set MAP_GRID_PEN_PRIMARY=thinnest,lightgrey
gmt set MAP_ANNOT_OFFSET_PRIMARY=1.0p
gmt set PS_LINE_CAP=butt
gmt set PS_LINE_JOIN=round

gmt set FONT_ANNOT_PRIMARY             = 7p,Helvetica,black
gmt set FONT_ANNOT_SECONDARY           = 7p,Helvetica,black
gmt set FONT_LABEL                     = 7p,Helvetica,black
gmt set FONT_LOGO                      = 7p,Helvetica,black
gmt set FONT_TITLE                     = 7p,Helvetica,black

# color0=20/20/20 # Budget
# color1=#4c78a8 # Observed
# color2=#f58518 # Steric
# color3=#54a24b # Mass
# color4=#9ecae9 # Altimetry
# color5=#72b7b2 # Glaciers
# color6=#b79a20 # GrIS
# color7=#f2cf5b # AIS
# color8=#b279a2 # TWS 
# color9=#e45756 #TWS natural
# color10=#9c755f # Dams
# color11=#79706e # GWD

color0=20/20/20 # Budget
color1=#0072B2 # Observed
color2=#2ca02c # Steric
color3=#D55E00 # Mass
color4=#9467bd # Altimetry
color5=#8c564b # Glaciers




Jopt=X6.9c/4c
Xjump=7.1c
Yjump=-4.3c

barprop=b0.13cb0
errprop=y+a+w2.5p+p0.4p

gmt psbasemap -K -R1900/2018/-200/50 -J$Jopt -X0.85c -Y4.6c -BWesn -Bx20g20 -By40g40+l'Height (mm)' > $ps
gmt psxy -R -J -O -K -L -t70 -G$color4 alt_glb_tseries_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color2 steric_glb_tseries_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color3 grd_glb_tseries_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color0 budget_tseries_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color1 obs_glb_tseries_ci.txt >> $ps

gmt psxy -R -J -O -K -W1.5p,$color2 steric_glb_tseries_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color3 grd_glb_tseries_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color0 budget_tseries_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color1 obs_glb_tseries_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color4 alt_glb_tseries_m.txt >> $ps
echo "1902 -195 a" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLB -O -K >> $ps

gmt psbasemap -O -K -R1900/2018/-200/50 -J$Jopt -X$Xjump -Bwesn -Bx20g20 -By40g40 >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color8 tws_glb_tseries_ci.txt  >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color7 AIS_glb_tseries_ci.txt  >> $ps
gmt psxy -R -J -O -K -L -t75 -G$color9 nat_glb_tseries_ci.txt  >> $ps
gmt psxy -R -J -O -K -L -t75 -G$color11 gwd_glb_tseries_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t75 -G$color10 dam_glb_tseries_ci.txt  >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color6 GrIS_glb_tseries_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color5 glac_glb_tseries_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color3 grd_glb_tseries_ci.txt  >> $ps

gmt psxy -R -J -O -K -W1.5p,$color7 AIS_glb_tseries_m.txt >> $ps
gmt set PS_LINE_CAP=round
gmt psxy -R -J -O -K -W1.5p,$color9,1_2:1 nat_glb_tseries_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color11,1_2:1 gwd_glb_tseries_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color10,1_2:1  dam_glb_tseries_m.txt >> $ps
gmt set PS_LINE_CAP=butt
gmt psxy -R -J -O -K -W1.5p,$color8  tws_glb_tseries_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color6 GrIS_glb_tseries_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color5 glac_glb_tseries_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color3 grd_glb_tseries_m.txt >> $ps
echo "1902 -195 b" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLB -O -K >> $ps

gmt psbasemap -O -K -R1900/2018/-1/4.0 -J$Jopt -X-$Xjump -Y$Yjump -BWeSn -Bx20g20 -By0.8g0.8+l'Trend (mm yr@+-1@+)' >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color2 steric_glb_sliding_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color3 grd_glb_sliding_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color0 budget_sliding_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color1 obs_glb_sliding_ci.txt >> $ps

gmt psxy -R -J -O -K -W1.5p,$color2 steric_glb_sliding_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color3 grd_glb_sliding_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color0 budget_sliding_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color1 obs_glb_sliding_m.txt >> $ps
echo "1902 -0.86 c" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLB -O -K >> $ps

gmt psbasemap -O -K -R1900/2018/-1/4.0 -J$Jopt -X$Xjump -BweSn -Bx20g20 -By0.8g0.8 >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color8 tws_glb_sliding_ci.txt  >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color7 AIS_glb_sliding_ci.txt  >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color6 GrIS_glb_sliding_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color5 glac_glb_sliding_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color3 grd_glb_sliding_ci.txt  >> $ps
gmt psxy -R -J -O -K -W1.5p,$color7 AIS_glb_sliding_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color8  tws_glb_sliding_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color6 GrIS_glb_sliding_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color5 glac_glb_sliding_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color3 grd_glb_sliding_m.txt >> $ps
echo "1902 -0.86 d" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLB -O -K >> $ps



gmt psbasemap -O -K -R0/1/-0.5/11.5 -JX3.2c/6.1c -Y1.1c -X$Xjump -T-100/-100/0.1 >> $ps

echo -e "0.0 11 \n 0.12 11" | psxy -R -J -O -K -W6p,$color1 -t70 >> $ps
echo -e "0.0 11 \n 0.12 11" | psxy -R -J -O -K -W1.5p,$color1 >> $ps
echo "0.13 11 Observed sea level" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 10 \n 0.12 10" | psxy -R -J -O -K -W6p,$color4 -t70 >> $ps
echo -e "0.0 10 \n 0.12 10" | psxy -R -J -O -K -W1.5p,$color4 >> $ps
echo "0.13 10 Altimetry" | gmt pstext -R -J -F+f7,Helvetica+jLM -O -K >> $ps

echo -e "0.0 9 \n 0.12 9" | psxy -R -J -O -K -W6p,$color0 -t70 >> $ps
echo -e "0.0 9 \n 0.12 9" | psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "0.13 9 Sum of contributors" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 8 \n 0.12 8" | psxy -R -J -O -K -W6p,$color2 -t70 >> $ps
echo -e "0.0 8 \n 0.12 8" | psxy -R -J -O -K -W1.5p,$color2 >> $ps
echo "0.13 8 Thermosteric" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 7 \n 0.12 7" | psxy -R -J -O -K -W6p,$color3 -t70 >> $ps
echo -e "0.0 7 \n 0.12 7" | psxy -R -J -O -K -W1.5p,$color3 >> $ps
echo "0.13 7 Barystatic" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps


echo -e "0.0 6 \n 0.12 6" | psxy -R -J -O -K -W6p,$color5 -t70 >> $ps
echo -e "0.0 6 \n 0.12 6" | psxy -R -J -O -K -W1.5p,$color5 >> $ps
echo "0.13 6 Glaciers" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 5 \n 0.12 5" | psxy -R -J -O -K -W6p,$color6 -t70 >> $ps
echo -e "0.0 5 \n 0.12 5" | psxy -R -J -O -K -W1.5p,$color6 >> $ps
echo "0.13 5 Greenland Ice Sheet" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 4 \n 0.12 4" | psxy -R -J -O -K -W6p,$color7 -t70 >> $ps
echo -e "0.0 4 \n 0.12 4" | psxy -R -J -O -K -W1.5p,$color7 >> $ps
echo "0.13 4 Antarctic Ice Sheet" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 3 \n 0.12 3" | psxy -R -J -O -K -W6p,$color8 -t70 >> $ps
echo -e "0.0 3 \n 0.12 3" | psxy -R -J -O -K -W1.5p,$color8 >> $ps
echo "0.13 3 Terrestrial water storage" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 2 \n 0.12 2" | psxy -R -J -O -K -W6p,$color9 -t70 >> $ps
gmt set PS_LINE_CAP=round
echo -e "0.0 2 \n 0.12 2" | psxy -R -J -O -K -W1.5p,$color9,1_2:0.5 >> $ps
gmt set PS_LINE_CAP=butt
echo "0.13 2 Natural TWS" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 1 \n 0.12 1" | psxy -R -J -O -K -W6p,$color10 -t70 >> $ps
gmt set PS_LINE_CAP=round
echo -e "0.0 1 \n 0.12 1" | psxy -R -J -O -K -W1.5p,$color10,1_2:0.5 >> $ps
gmt set PS_LINE_CAP=butt
echo "0.13 1 Dam impoundment" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 0 \n 0.12 0" | psxy -R -J -O -K -W6p,$color11 -t70 >> $ps
gmt set PS_LINE_CAP=round
echo -e "0.0 0 \n 0.12 0" | psxy -R -J -O -K -W1.5p,$color11,1_2:0.5 >> $ps
gmt set PS_LINE_CAP=butt
echo "0.13 0 Groundwater depletion" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps

gmt psconvert $ps -Tf
