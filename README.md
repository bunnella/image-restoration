# image-restoration
######Alec Bunnell, Anne Grosse, Jialun Luo, Sam Spaeth

######Carleton College 2015-16 Math/Stats Comps

#### Use:
File imageRestorationSUPERSPEEDY.R should contain everything needed to process an image from start to finish. Call `R CMD SHLIB speedy.c` to compile C dependencies (pre-built speedy.dll works on 64-bit Windows).

#### TODO:
+ Tailor data term more to JPEG
+ Compare output to original quantitatively
+ Parameter optimization (mega-run script?)
+ Get rid of random black/white spaces
+ Figure out why it crashes sometimes on repeated calls to the same chain...
+ Color!

#### Think about:
+ Jitter gamma (microedge penalty)
+ Continuous sampling (have Mathematica try to integrate H)
+ Pre-analyze (8-bit grids)
+ Texture segmentation
+ More parameters with Greek letter names wooo

#### Code optimizations (low priority):
+ Parallelize?
+ Don't check who's on the edge on every single energy function call
+ Come up with some clever algorithm for sampling without having to brute force every (discrete) level
