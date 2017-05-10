using SIT
using PyPlot

# 2D
xin= Array(linspace(1,20,10)); 
xout=Array(linspace(1,20,20));
zin= Array(linspace(1,20,15)); 
zout=Array(linspace(1,20,30));
yin=randn(length(zin), length(xin));
yout=zeros(length(zout),length(xout));
SIT.Interpolation.interp_spray!(xin,zin, yin, xout,zout, yout, :interp,:B2)
subplot(131);
imshow(yin, aspect="auto"); 
subplot(132);
imshow(yout, aspect="auto")
subplot(133);
imshow(itp[zout, xout], aspect="auto")


# 1D
using SIT
using PyPlot
reload("SIT")
xin= Array(linspace(1,10,10)); 
xout=Array(linspace(1,10,40));
yin=randn(length(xin));
# to get kernel
#yin = zeros(length(xin));yin[6]=1.0
plot(xin, yin, "r*-")
yout=zeros(length(xout));
SIT.Interpolation.interp_spray!(xin, yin, xout, yout, :interp)
plot(xout, yout, "b*-")
SIT.Interpolation.interp_spray!(xin, yin, xout, yout, :interp, :B2)
plot(xout, yout, "g*-")

