import math
 
def pointsOnSphere(N):
    N = float(N) # in case we got an int which we surely got
    pts = []
 
    inc = math.pi * (3 - math.sqrt(5))
    off = 2 / N
    for k in range(0, N):
        y = k * off - 1 + (off / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        pts.append([math.cos(phi)*r, y, math.sin(phi)*r])
 
    return pts

points=pointsOnSphere(256)

print points
