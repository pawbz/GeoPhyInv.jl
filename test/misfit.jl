
using PyPlot

range=10*pi;
x = linspace(-range, range, 1000);
sig1 = sin(x);
sig2 = sin(x-0.33*pi)

subplot(411)
plot(sig1,color="red");
plot(sig2,color="blue");
title("two sinusoids with a phase shift")

subplot(412)
plot(sig1.*sig2)
title("product of sinusoids; sum along horizontal axis for LSC ")


subplot(413)
sigdiff  =  sig1 - sig2;
plot(sigdiff,color="red"); 
plot((-1.0 * (sigdiff).^2.0), color="blue");

title("their difference and squared difference; needed for NLSC")

subplot(414)
plot(exp(-1.0 * (sigdiff).^2.0), color="green");
plot(exp(-4.0 * (sigdiff).^2.0), color="orange");
plot(exp(-16.0 * (sigdiff).^2.0), color="black");
plot(exp(-1000.0 * (sigdiff).^2.0), color="black");
title("exponent of squared difference for different values of sigma for NLSC")

show()


