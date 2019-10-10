\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
 
Subject To
 R0: - 4 x0 - 4 x1 - y0 - y1 + alpha >= -8
 R1: - x1 <= 0
 R2: - 2 x0 - 2 x1 - y0 - y1 + alpha >= -2
 R3: - 2 x0 - x1 <= -3
 R4: x0 - 2 y0 >= 0
 R5: x0 - x1 + 3 y0 <= 3
 R6: x0 + y0 >= 1
 R7: x1 - y1 >= 0
 R8: x0 + x1 - 3 y0 >= 0
 R9: y0 + y1 >= 1
Bounds
 x0 <= 4
 x1 <= 4
Binaries
 y0 y1
End
