import numpy as np

def rk4( dx, dy, dz, x0, y0, z0, t , n):
    h   = np.linspace(t[0],t[1],n)
    h0  = abs(h[1] - h[0])
    x   = np.zeros(n)
    x[0]= x0
    y   = np.zeros(n)
    y[0]= y0
    z   = np.zeros(n)
    z[0]= z0
    kx, kx1, kx2, kx3, kx4, = [0,0,0,0,0]
    ky, ky1, ky2, ky3, ky4, = [0,0,0,0,0]
    kz, kz1, kz2, kz3, kz4, = [0,0,0,0,0]
    
    for i in range(n-1):
        kx1 = dx( h[i], x[i], y[i], z[i] )
        ky1 = dy( h[i], x[i], y[i], z[i] )
        kz1 = dz( h[i], x[i], y[i], z[i] )
        
        kx2 = dx( h[i] + h0/2 , x[i] + kx1*h0/2, y[i] +ky1*h0/2, z[i] + kz1*h0/2 )
        ky2 = dy( h[i] + h0/2 , x[i] + kx1*h0/2, y[i] +ky1*h0/2, z[i] + kz1*h0/2 )
        kz2 = dz( h[i] + h0/2 , x[i] + kx1*h0/2, y[i] +ky1*h0/2, z[i] + kz1*h0/2 )
        
        kx3 = dx( h[i] + h0/2 , x[i] + kx2*h0/2, y[i] +ky2*h0/2, z[i] + kz2*h0/2 )
        ky3 = dy( h[i] + h0/2 , x[i] + kx2*h0/2, y[i] +ky2*h0/2, z[i] + kz2*h0/2 )
        kz3 = dz( h[i] + h0/2 , x[i] + kx2*h0/2, y[i] +ky2*h0/2, z[i] + kz2*h0/2 )
        
        kx4 = dx( h[i] + h0   , x[i] + kx3*h0  , y[i] +ky3*h0  , z[i] + kz3*h0   )
        ky4 = dy( h[i] + h0   , x[i] + kx3*h0  , y[i] +ky3*h0  , z[i] + kz3*h0   )
        kz4 = dz( h[i] + h0   , x[i] + kx3*h0  , y[i] +ky3*h0  , z[i] + kz3*h0   )
        
        kx   =( kx1 + 2*kx2 + 2*kx3 +kx4 )/6
        ky   =( ky1 + 2*ky2 + 2*ky3 +ky4 )/6
        kz   =( kz1 + 2*kz2 + 2*kz3 +kz4 )/6
        
        x[i+1] = x[i] + kx * h0
        y[i+1] = y[i] + ky * h0
        z[i+1] = z[i] + kz * h0

    return h, x ,y ,z
        
        
        