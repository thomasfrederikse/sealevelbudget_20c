
rm *.ps
ps=Mass_steric_ratio.ps
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_MEDIA=Custom_18.3cx4.6c
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_FRAME_PEN=0.4p,40/40/40
gmt set MAP_FRAME_TYPE=plain
gmt set MAP_ANNOT_OFFSET_PRIMARY=1.0p
gmt set MAP_LABEL_OFFSET=0.7p
gmt set MAP_TITLE_OFFSET=0.05c
gmt set MAP_TICK_LENGTH_PRIMARY=0.03c
gmt set MAP_TICK_PEN=thinnest,40,40,40
gmt set MAP_GRID_PEN_PRIMARY=thinnest,lightgrey
gmt set PS_LINE_CAP=butt
gmt set PS_LINE_JOIN=round

gmt set FONT_ANNOT_PRIMARY             = 7p,Helvetica,black
gmt set FONT_ANNOT_SECONDARY           = 7p,Helvetica,black
gmt set FONT_LABEL                     = 7p,Helvetica,black
gmt set FONT_LOGO                      = 7p,Helvetica,black
gmt set FONT_TITLE                     = 7p,Helvetica,black

color0=20/20/20 # Budget
color1=#0072B2 # Observed
color2=#ff7f0e # Steric
color3=#d62728 # Mass
color4=#56B4E9 # Altimetry
color5=#9467bd # Glaciers
color6=#8c564b # GrIS
color7=#e377c2 # AIS
color8=#bcbd22 # TWS 
color9=#2ca02c #TWS natural
color10=#17becf # Dams
color11=#7f7f7f # GWD

Jopt=X6.7c/4c
Xjump=7.3c

gmt psbasemap -K -R1900/2018/-1.5/1.5 -J$Jopt -X0.9c -Y0.3c -BWeSn+t'All terms included' -Bx20g20 -By0.4g0.4+l'Fraction (-)' > $ps
gmt psxy -R -J -O -K -L -t70 -G$color8 frac_tws_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color7 frac_AIS_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color6 frac_GrIS_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color5 frac_glac_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color3 frac_total_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color2 frac_steric_ci.txt >> $ps

gmt psxy -R -J -O -K -W1.5p,$color8 frac_tws_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color7 frac_AIS_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color6 frac_GrIS_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color5 frac_glac_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color3 frac_total_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color2 frac_steric_m.txt >> $ps
echo "1902 -1.35 a" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLM -O -K >> $ps

gmt psbasemap -O -K -R1900/2018/-0.5/1.5 -J$Jopt -X$Xjump -BWeSn+t'No TWS' -Bx20g20 -By0.3g0.3 >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color7 frac_AIS_notws_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color6 frac_GrIS_notws_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color5 frac_glac_notws_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color3 frac_total_notws_ci.txt >> $ps
gmt psxy -R -J -O -K -L -t70 -G$color2 frac_steric_notws_ci.txt >> $ps

gmt psxy -R -J -O -K -W1.5p,$color7 frac_AIS_notws_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color6 frac_GrIS_notws_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color5 frac_glac_notws_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color3 frac_total_notws_m.txt >> $ps
gmt psxy -R -J -O -K -W1.5p,$color2 frac_steric_notws_m.txt >> $ps
echo "1902 -0.4 b" | gmt pstext -R -J -F+f8,Helvetica-Bold+jLM -O -K >> $ps

gmt psbasemap -O -K -R0/1/-0.5/5.5 -JX3.2c/2.5c -Y0.75c -X6.9c -T-100/-100/0.1 >> $ps
echo -e "0.0 5 \n 0.12 5" | psxy -R -J -O -K -W6p,$color2 -t70 >> $ps
echo -e "0.0 5 \n 0.12 5" | psxy -R -J -O -K -W1.5p,$color2 >> $ps
echo "0.13 5 Thermosteric" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 4 \n 0.12 4" | psxy -R -J -O -K -W6p,$color3 -t70 >> $ps
echo -e "0.0 4 \n 0.12 4" | psxy -R -J -O -K -W1.5p,$color3 >> $ps
echo "0.13 4 Barystatic" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 3 \n 0.12 3" | psxy -R -J -O -K -W6p,$color5 -t70 >> $ps
echo -e "0.0 3 \n 0.12 3" | psxy -R -J -O -K -W1.5p,$color5 >> $ps
echo "0.13 3 Glaciers" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 2 \n 0.12 2" | psxy -R -J -O -K -W6p,$color6 -t70 >> $ps
echo -e "0.0 2 \n 0.12 2" | psxy -R -J -O -K -W1.5p,$color6 >> $ps
echo "0.13 2 Greenland Ice Sheet" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 1 \n 0.12 1" | psxy -R -J -O -K -W6p,$color7 -t70 >> $ps
echo -e "0.0 1 \n 0.12 1" | psxy -R -J -O -K -W1.5p,$color7 >> $ps
echo "0.13 1 Antarctic Ice Sheet" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
echo -e "0.0 0 \n 0.12 0" | psxy -R -J -O -K -W6p,$color8 -t70 >> $ps
echo -e "0.0 0 \n 0.12 0" | psxy -R -J -O -K -W1.5p,$color8 >> $ps
echo "0.13 0 Terrestrial water storage" | gmt pstext -R -J -F+f7+jLM -O -K >> $ps
gmt psconvert $ps -Tf
