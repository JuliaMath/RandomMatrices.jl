#bellcurve2.jl
#Algorithm 2.1 of Random Eigenvalues by Alan Edelman

#Experiment:    Generate random samples from the normal distribution
#Plot:          Histogram of random samples
#Theory:        The normal distribution curve


## Experiment
t = 1000000     # Trials
dx = 0.25       # binsize
v = randn(t,2)  # samples
x = -4:dx:4   # range

## Plot
using Winston
module Plot3D
    # Warning: the 3D interface is very new at the time
    # of writing this and is likely to not be available
    # and/or have changed
    if Winston.output_surface == :gtk
        include(Pkg.dir("Winston","src","canvas3d_gtk.jl"))
    else
        include(Pkg.dir("Winston","src","canvas3d.jl"))
    end
end
edg1, edg2, count = hist2(v, x-dx/2)
cntrs = (first(x)+dx/2, last(x)-dx/2, length(x)-1)
x,y = meshgrid(cntrs, cntrs)
Plot3D.plot3d(Plot3D.surf(x,y,count./(t*dx^2) *10*length(cntrs)))
Plot3D.plot3d(Plot3D.surf(x,y,exp(-(x.^2+y.^2)/2)/(2*pi) *10*length(cntrs)))
