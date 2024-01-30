function res=G(l,f1,b)
if(f1>=0)
   res=b^(-f1-1)*gamma(f1+1)*gammainc(b*l,f1+1,'upper'); 
elseif(f1==-1)
    res=-ei(-b*l);
else
    res1=(-1)^(-f1)*ei(-b*l)/factorial(-f1-1)*b^(-f1-1);
    res2=0;
    for k1=0:-f1-2
       res2=res2+exp(-b*l)/(l^(-f1-1))*(-1)^k1*b^k1*l^k1*factorial(-f1-2-k1)/(factorial(-f1-1));
    end
    res=res1+res2;
end
end