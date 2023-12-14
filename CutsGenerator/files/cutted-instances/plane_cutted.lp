\ File written by HiGHS .lp file handler
min
 obj: +1 y1 +1 y2 +1 y3 +1 y4 
st
 a1: +1 x1 -1 x2 -1 x3 = +0
 a2: +1 x2 -1 x3 +1 x4 = +0
 a3: +1 x1 -1 x3 = +0
 b1: -82 y1 +1 x1 <= +1
 b2: +82 y2 -1 x2 >= -1
 b3: -82 y3 -1 x3 <= +1
 b4: -82 y4 +1 x4 >= +1
bounds
 y1 <= 1
 y2 <= 1
 y3 <= 1
 y4 <= 1
bin
 y1
 y2
 y3
 y4
gen
semi
end
