# image-restoration
######Alec Bunnell, Anne Grosse, Jialun Luo, Sam Spaeth

######Carleton College 2015-16 Math/Stats Comps

#### Use:
File imageRestorationSUPERSPEEDY.R should contain everything needed to process an image from start to finish; experimental color code is in imageRestorationSUPERSPEEDYcolor.R. Call `R CMD SHLIB speedy.c` to compile C dependencies (pre-built speedy.dll works on 64-bit Windows).

#### TODO:
+ Nothing woo!

#### Think about:
+ Pre-analyze (8-bit grids)
+ More JPEG-oriented trust term/consistency measure (8x8?)
+ Parameter optimization
+ (More advanced) color!
+ Continuous jitter?

#### Code optimizations (low priority):
+ Parallelize?
+ Don't check who's on the edge on every single energy function call
