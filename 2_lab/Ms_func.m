function Ms_func = Ms_func(h,m,z_p,c1,r2,w1,r1,x_s1,x_s2,U_work,M_se)
    kr = h*(sinh(2*h)+sin(2*h))/(cosh(2*h)-cos(2*h));
    kx = 3/(2*h)*(sinh(2*h)-sin(2*h))/(cosh(2*h)-cos(2*h));
    Ms_func = z_p*m*(U_work)^2*r2*kr/(w1*((r1+c1*r2*kr)^2+(x_s1+c1*x_s2*kx)^2)) - M_se;
end