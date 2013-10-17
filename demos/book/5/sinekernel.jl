#sinekernel.jl

n=200
x=(0:.025:1)*3*pi
include("hermitekernel.jl")
K=hermitekernel(n-1,x/sqrt(2*n))*pi/sqrt(2*n)

## Plot
using Winston
module Plot3D
    using Winston
    # Warning: the 3D interface is very new at the time
    # of writing this and is likely to not be available
    # and/or have changed
    if Winston.output_surface == :gtk
        include(Pkg.dir("Winston","src","canvas3d_gtk.jl"))
    else
        include(Pkg.dir("Winston","src","canvas3d.jl"))
    end
end
xx = Float64[xx for xx = x, yy = x]
yy = Float64[yy for xx = x, yy = x]
Plot3D.plot3d(Plot3D.surf(xx,K,yy))

Plot3D.plot3d(Plot3D.surf(
    (u,v)->u,
    (u,v)->sin(u-v)./(u-v),
    (u,v)->v,
    x, x))
