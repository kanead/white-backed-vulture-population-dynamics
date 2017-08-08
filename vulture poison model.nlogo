;;-------------------------------------------------------------;;
;;-------------------------------------------------------------;;
; http://jasss.soc.surrey.ac.uk/20/1/3.html
breed [adults adult]
breed [subadults subadult]
breed [juveniles juvenile]

breed [carcasses carcass]

adults-own [energy x0 y0 xcar ycar]
subadults-own [energy x0 y0 xcar ycar target-patch]
juveniles-own [energy x0 y0 xcar ycar target-patch]
carcasses-own [mass decay]

globals [day poison]

;;-------------------------------------------------------------;;
;;----------------------- SETUP COMMANDS---------------------- ;;
;;-------------------------------------------------------------;;
to setup

  clear-all
ask patches [set pcolor green - 1]
 ask patch 99.5 99.5
  [
   set pcolor yellow
    ask patches in-radius 8 [set pcolor yellow]
    ask patches in-radius 50 with [pcolor != yellow] [set pcolor green - 3]
    ask  patches in-radius 80 with [pcolor != yellow and pcolor != green - 3] [set pcolor green]
  ]

ask n-of roosts patches [set pcolor brown]

ask n-of n-adults patches with [pcolor = yellow] [sprout-adults 1 [
   set color red
   set size 2
  set energy 10000
  set x0 xcor
 set y0 ycor]
]


ask n-of n-subadults patches  [sprout-subadults 1 [
    set color brown
    set size 2
    set energy 10000
    set x0 xcor
    set y0 ycor]
]


ask n-of n-juveniles patches  [sprout-juveniles 1 [
    set color blue
    set size 2
    set energy 10000
    set x0 xcor
    set y0 ycor]
]



  reset-ticks
end

;;-------------------------------------------------------------;;
;;----------------------- GO COMMANDS--------------------------;;
;;-------------------------------------------------------------;;
to go
  if ticks = day-length  [set day day + 1 create-next-day]


;  if ticks = 0 [
;create-carcasses n-carcasses [
;setxy random-xcor random-ycor
;set color white
;set size 2
;set shape "circle"
;]]
  if ticks = 0 [
  ask patches with [pcolor != green - 1] [
    while [sum [mass] of carcasses  < random-normal 3000 100] [
sprout-carcasses 1 [
set mass random-gamma alpha beta ;; 1.2 0.004
move-to one-of  patches with [pcolor != green - 1]
;setxy random-pxcor random-pycor
set color white
set size 2
set shape "circle"
]]]
  ]


    if ticks = 0 [
  ask patches with [pcolor = green - 1] [
    while [sum [mass] of carcasses < random-normal 9000 100] [
sprout-carcasses 1 [
set mass random-normal mu std ; 500 100
move-to one-of  patches with [pcolor = green - 1]
;setxy random-pxcor random-pycor
set color white
set size 2
set shape "circle"
]]]
  ]

    ask carcasses [
  check
  if distancexy 99.5 99.5 < 80 and mass > 1000 [set color red]
]


if ticks = 1[
  ask carcasses [ifelse distancexy 99.5 99.5 > 80 [if random outside-rate = 1 [set color cyan]]
    [if random inside-rate = 1 [set color cyan]]
      ]]

if ticks = 1[
ask carcasses
[ifelse color != red [set decay decay + 2] [set decay decay + 1]
]]





  ask adults
 [forage-vul
  feed-vul
  social-vul
  rtb-vul
  territory-vul
  ]


 ask subadults
 [forage-vul
   feed-vul
   social-vul
   rtb-sub
  set x0 xcor
  set y0 ycor
  ]


 ask juveniles
 [forage-vul
   feed-vul
   social-vul
   rtb-juvenile
   set x0 xcor
  set y0 ycor
  ]


;if not any? adults or not any? juveniles [ print "day" print day print "n-adults" print count adults print "n-juveniles" print count juveniles stop]

 tick
end


;;-------------------------------------------------------------;;
;;--------------------- ADULT COMMANDS-------------------------;;
;;-------------------------------------------------------------;;
to forage-vul
    set energy  energy  - 1
  fd v
   if random 600 = 1 ;; frequency of turn
  [ ifelse random 2 = 0 ;; 50:50 chance of left or right
    [ rt 30 ] ;; could add some variation to this with random-normal 45 5
    [ lt 30 ]] ;; so that it samples from a dist with mean 45 SD 5
 end
  to territory-vul
    while [[pcolor] of patch-here = green ]
     [
       face min-one-of patches with [pcolor = green - 3  ] [ distance myself ]
       forage-vul
     ]

  end

  to feed-vul
    if energy > 0 [
let target-food min-one-of (carcasses with [shape = "circle"] in-radius vision) [distance myself]
  if target-food != nobody  [
    move-to target-food
    set xcar xcor
    set ycar ycor
    if [color] of target-food = cyan  [die]
  ]]
  end


    to social-vul
    if energy > 0 [
let target-food min-one-of (carcasses with [shape = "target"] in-radius 7) [distance myself]
  if target-food != nobody  [
    move-to target-food
    set xcar xcor
    set ycar ycor
    if [color] of target-food = cyan  [die]
  ]]
  end

to rtb-vul
  if energy <= 0 [
   face patch x0 y0
   fd v * 2
  ]
end

;;-------------------------------------------------------------;;
;;------------------- SUBADULT COMMANDS------------------------;;
;;-------------------------------------------------------------;;

to rtb-sub
  if energy = 0 [ ifelse random 2 = 1 [set target-patch min-one-of (patches  with [pcolor = yellow]) [distance myself]]
    [
  set target-patch min-one-of (patches  with [pcolor = brown]) [distance myself]
    ]
  ]
  if energy < 0   [ face target-patch
  fd v * 2
  ]
end


;;-------------------------------------------------------------;;
;;------------------- JUVENILE COMMANDS------------------------;;
;;-------------------------------------------------------------;;

to rtb-juvenile
  if energy = 0 [
  set target-patch min-one-of (patches  with [pcolor = brown]) [distance myself]
  ]
  if energy < 0   [ face target-patch
      fd v * 2
    ]
end

;;-------------------------------------------------------------;;
;;------------------- CARCASS COMMANDS-------------------------;;
;;-------------------------------------------------------------;;
to check
ifelse any? turtles-here with [shape = "default"] [set shape "target"][set shape "circle"]
end

;;-------------------------------------------------------------;;
;;------------------- GENERAL COMMANDS-------------------------;;
;;-------------------------------------------------------------;;
to create-next-day
  clear-links
  reset-ticks
     ask adults [set energy 10000
        face patch xcar ycar]
  ask subadults [set energy 10000
        face patch xcar ycar]
  ask juveniles [set energy 10000
        face patch xcar ycar]
  ask carcasses with [decay = 2] [die]
  go
end
@#$#@#$#@
GRAPHICS-WINDOW
211
15
721
546
-1
-1
2.5
1
10
1
1
1
0
1
1
1
0
199
0
199
0
0
1
ticks
30.0

BUTTON
15
16
82
49
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
88
16
151
49
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
6
60
178
93
N-adults
N-adults
0
100
0
1
1
NIL
HORIZONTAL

MONITOR
96
238
153
283
NIL
day
17
1
11

INPUTBOX
8
236
86
296
day-length
20000
1
0
Number

INPUTBOX
5
366
78
426
v
0.00667
1
0
Number

SLIDER
7
130
179
163
N-juveniles
N-juveniles
0
100
0
1
1
NIL
HORIZONTAL

SLIDER
6
203
178
236
vision
vision
0
10
6
1
1
NIL
HORIZONTAL

MONITOR
87
457
202
502
N poisoned carcasses
count carcasses with [color = cyan]
17
1
11

MONITOR
100
359
193
404
NIL
count adults
17
1
11

MONITOR
100
409
196
454
NIL
count juveniles
17
1
11

TEXTBOX
1090
10
1779
790
NOTES\n\nHABITAT\n- Total area is 200 x 200km = 40,000km^2\n- 2 habitat types, inside + outside of protected area, adults don't go beyond dark green radius (r=50km) (Spiegel et al 2013)\n- Larger light green circle is the extent of the park, represents Kruger with area = 20,000km^2\n- Remainder of area is non-protected = 20,000km^2\n- small yellow circle = adult roost\n\nANIMALS\n- N-adults = number of adult birds\n- N-juveniles = number of juvenile birds\n- N-carcasses = number of carcasses\n- vision = distance (km) at which a bird can detect a carcass (Kane & Kendall 2017)\n- day-length = length of one foraging day in seconds (Spiegel et al 2013)\n- Distance travelled = 120 km so mean speed = 120/5 = 24km/hr where 5 is the time between successive roosts (Spiegel et al 2013)\n- v = speed in km/s; 0.0067km/s = 24km/hr\n\nCARCASSES\n0.15kg of carcass per km^2 (Murn & Anderson 2008) * area of park (20,000km^2) = 3,000kg carrion in area. Distributed according to a Gamma dist. which allows for the occasional large carcass. \n\nOutside of the park there is the same area. Campbell reckons density is higher, let's say double so 0.3kg of carcass per km^2 = 6,000 kg of carrion. Here it is distributed according to a normal distribution as most carrion will be domestic animals like cows. \n\nPOISON\n- inside-rate and outside-rate are the rates at which a carcass can be poisoned inside and outside the green circle respectively. 5 means a 1 in 5 chance of a carcass being poisoned.\n\nREFERENCES\n(a) Kane, A., & Kendall, C. J. (2017). Understanding how mammalian scavengers use information from avian scavengers: cue from above. Journal of Animal Ecology, 86(4), 837-846.\n(b) Spiegel, O., Getz, W. M., & Nathan, R. (2013). Factors influencing foraging search efficiency: why do scarce lappet-faced vultures outperform ubiquitous white-backed vultures?. The American Naturalist, 181(5), E102-E115.\n(c) Murn, C., & Anderson, M. D. (2008). Activity patterns of African White-backed Vultures Gyps africanus in relation to different land-use practices and food availability. Ostrich-Journal of African Ornithology, 79(2), 191-198.\n\n
14
0.0
1

CHOOSER
6
308
98
353
inside-rate
inside-rate
2 5 10 15 20 25 30 35 40 50 75 100 200 500 1000
13

CHOOSER
102
309
194
354
outside-rate
outside-rate
2 5 10 15 20 25 30 35 40 50 75 100 200 500 1000
13

PLOT
733
18
1061
375
Time series of bird population
time
Count birds
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"adults" 1.0 0 -2674135 true "" "plot count adults"
"juveniles" 1.0 0 -13345367 true "" "plot count juveniles"
"subadults" 1.0 0 -6459832 true "" "plot count subadults"

MONITOR
90
507
201
552
total carrion mass
sum [mass] of carcasses
1
1
11

SLIDER
8
96
180
129
N-subadults
N-subadults
0
100
0
1
1
NIL
HORIZONTAL

MONITOR
93
560
203
605
max carcass size
max [mass] of carcasses
1
1
11

INPUTBOX
250
581
300
641
alpha
1.2
1
0
Number

INPUTBOX
310
581
378
641
beta
0.0040
1
0
Number

INPUTBOX
552
584
602
644
mu
500
1
0
Number

INPUTBOX
612
583
662
643
std
100
1
0
Number

TEXTBOX
222
556
435
590
Parameters of Gamma distribution
14
0.0
1

TEXTBOX
516
557
774
591
Parameters of normal distribution
14
0.0
1

SLIDER
5
166
177
199
roosts
roosts
0
100
10
1
1
NIL
HORIZONTAL

TEXTBOX
222
660
713
702
http://homepage.divms.uiowa.edu/~mbognar/applets/gamma.html
11
0.0
1

@#$#@#$#@

## Vulture Data
20–50 WBVs per 100 km2; Mundy et al. 1992
Assume 50 per 100km^2

Mean Distance Travelled 120 km
Therefore average speed is 120/5 = 24km/hr (0.00666667 km/second)
Daily flight duration ~ 5 hours

Home range for Cape Vulture ~ 38,000 km^2, radius of 110 km to get this area in a circle, Kendall & Virani 2012


## Carcass Data

4.38kg of carrion per km^2 per day
0.12-0.15 kg km-2 day-1 from Campbell

##Area Data
Kruger is 20,000 km2
square with side of length 141 km

Circle with the same area has a radius = 80km

@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
0
Rectangle -7500403 true true 151 225 180 285
Rectangle -7500403 true true 47 225 75 285
Rectangle -7500403 true true 15 75 210 225
Circle -7500403 true true 135 75 150
Circle -16777216 true false 165 76 116

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.2.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>report mean distances</metric>
    <enumeratedValueSet variable="vision">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="excretion-rate">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-hyenas">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="day-length">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-carcasses">
      <value value="22"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-vultures">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="diffusion-rate">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nutrient-heat">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v-hyena">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v-vulture">
      <value value="0.0125"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
