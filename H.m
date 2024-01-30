function res=H(a,b,c,d,e,f,g,h)
if((c/d)==(f/g))
    res=0;
    for k2=0:a
        res=res+d^(-e)*g^(-h)*nchoosek(a,k2)*(-c/d)^(a-k2)*exp(b*c/d)*G(c/d,k2-e-h,b);
    end
    
else
    res1=0;
    res2=0;
    for n1=1:e
        syms v;
        fv=1/(v+f/g)^h;
        fv_d=diff(fv,v,e-n1);
        fai1=1/factorial(e-n1)*subs(fv_d,v,-c/d);
       for t3=0:a
          res1=res1+fai1*nchoosek(a,t3)*exp(b*c/d)*(-c/d)^(a-t3)*G(c/d,t3-n1,b); 
       end
    end
    for n2=1:h
        syms v;
        fv=1/(v+c/d)^e;
        fv_d=diff(fv,v,h-n2);
        fai2=1/factorial(h-n2)*subs(fv_d,v,-f/g);
       for t2=0:a
           res2=res2+fai2*nchoosek(a,t2)*exp(b*f/g)*(-f/g)^(a-t2)*G(f/g,t2-n2,b); 
       end
    end
    res=d^(-e)*g^(-h)*(res1+res2);
end
end