function res=Q(a,b,c,d,e,f,g)
res=c^(-d)*f^(-g)*(b/c)^(-d)*(e/f)^(a+1-g)*beta(a+1,d-(a+1)+g)*hypergeom([d,a+1],d+g,1-e*c/(f*b));
% syms v;
% fv=v^a*(b+c*v)^(-d)*(e+f*v)^(-g);
% res=double(int(fv,0,inf));
end