clear all;
clc;


L=3000000;

 m_AU=2;
 omiga_AU=1;
 lambda_AU=m_AU/omiga_AU;
 beta_AU=lambda_AU^m_AU/gamma(m_AU);
 
 m_AE=2;
 omiga_AE=1;
 lambda_AE=m_AE/omiga_AE;
 beta_AE=lambda_AE^m_AE/gamma(m_AE);
 
 m_JU=2;
 omiga_JU=0.1;
 lambda_JU=m_JU/omiga_JU;
 beta_JU=lambda_JU^m_JU/gamma(m_JU);
 
 m_JE=2;
 omiga_JE=0.2;
 lambda_JE=m_JE/omiga_JE;
 beta_JE=lambda_JE^m_JE/gamma(m_JE);
 
 
 tao=0.3;
 Rs=1;
 h=waitbar(0,'please wait');
kk=1;
for SNR=0:4:32
    str=['运行中...',num2str(kk/10*100),'%'];
    waitbar(kk/10,h,str)
Gamma_B=10^(SNR/10);

u=0.5;
 M=2;
N=2;
   Gamma_J=Gamma_B*u;

     h_bu = gamrnd(m_AU,omiga_AU/m_AU,N,L);
    h_be = gamrnd(m_AE,omiga_AE/m_AE,M,L);   %
    h_ju = gamrnd(m_JU,omiga_JU/m_JU,1,L);
    h_je = gamrnd(m_JE,omiga_JE/m_JE,M,L);
    delta1=(2^Rs-1)/Gamma_B;
    
    afa1=2^Rs*omiga_AU/Gamma_J;
    delta2=(2^Rs-1)*u;
   
 %%%simulation
 %%NJ
 
 Pout_NJ_s(kk)=0;
 for i=1:L
     [max_A,index_A]=max(h_bu(:,i));
     h_AE=h_be(:,i);
     [max_E,index_E]=max(h_AE);
     if(max_A<delta1+2^Rs*max_E)
         Pout_NJ_s(kk)=Pout_NJ_s(kk)+1/L;
     end
 end
    
 
 %%%AJ
 
 
   Pout_AJ_s(kk)=0;
   for i=1:L
       [max_A,index_A]=max(h_bu(:,i));
       h_AE=h_be(:,i);
     [max_E,index_E]=max(h_AE);
    C_bu_aj=log2(1+Gamma_B*max_A/(Gamma_J*h_ju(i)+1));
   C_be_aj=log2(1+Gamma_B*h_be(:,i)./(Gamma_J*h_je(:,i)+1));
        [max_E2,index_E2]=max(C_be_aj);
      if(h_ju(i)<tao)
          if(C_bu_aj-max_E2<Rs)
             Pout_AJ_s(kk)=Pout_AJ_s(kk)+1/L;
          end
      else
          if(max_A<delta1+2^Rs*max_E)
              Pout_AJ_s(kk)=Pout_AJ_s(kk)+1/L;
          end
      end
   end
   
   %%%IC
   Pout_IC_s(kk)=0;
   Pout_IC_s_high(kk)=0;
   for i=1:L
       [max_A,index_A]=max(h_bu(:,i));
           C_bu_ic=log2(1+(Gamma_B-omiga_JU/omiga_AU*Gamma_J)*max_A);
   C_be_ic=log2(1+(Gamma_B*omiga_AU-Gamma_J*omiga_JU)*h_be(:,i)./(Gamma_J*(h_be(:,i)*h_ju(i)+h_je(:,i)*max_A)+omiga_AU));
        [max_E3,index_E3]=max(C_be_ic);
      if(C_bu_ic-max_E3<Rs)
          Pout_IC_s(kk)=Pout_IC_s(kk)+1/L;
      end
      [min_E4,index_E4]=min(h_je(:,i)./h_be(:,i));
      if(max_A*(h_ju(i)+min_E4*max_A)<afa1)
          Pout_IC_s_high(kk)=Pout_IC_s_high(kk)+1/L;
      end
   end
 
 %%%%numerical
 
 %% NJ
 Pout_NJ_t(kk)=0;

 for n1=0:N 
     for n2=0:N
         for n3=0:N
            if (n1+n2+n3==N)
               B1=lambda_AU*(N-n3);
               C1=n2;
               D1=factorial(N)/(factorial(n1)*factorial(n2)*factorial(n3))*(-1)^n1*(-lambda_AU)^n2;
               for m1=0:M
                  for m2=0:M
                     for m3=0:M
                        if(m1+m2+m3==M)
                           B2=lambda_AE*(M-m3);
                           C2=m2;
                           D2=factorial(M)/(factorial(m1)*factorial(m2)*factorial(m3))*(-1)^m1*(-lambda_AE)^m2;
                           for t=0:C1
                               temp1=D1*nchoosek(C1,t)*delta1^(C1-t)*2^(t*Rs)*exp(-B1*delta1);
                               temp2=0;
                               temp3=0;
                               if(B2~=0)
                                    temp2=D2*B2*factorial(t+C2)*(B2+B1*2^Rs)^(-t-C2-1);
                               end
                               if(C2~=0)
                                   temp3=D2*C2*factorial(t+C2-1)*(B2+B1*2^Rs)^(-t-C2);
                               end
                               Pout_NJ_t(kk)=Pout_NJ_t(kk)+temp1*(temp3-temp2);
                           end
                           
                        end
                     end
                  end
               end
            end
         end
     end
 end
%  
 
   %% AJ
 P2_1=0;
 for k=0:m_JU-1
    P2_1=P2_1+lambda_JU^k/factorial(k)*tao^k*exp(-lambda_JU*tao);
 end
P2=Pout_NJ_t(kk)*P2_1;
%  P1_1=0;
%  P1_2=0;
%  P1_3=0;
%  P1_4=0;
%  P1_5=0;
%  P1_6=0;
%  P1_7=0;
%   for n1=0:N 
%      for n2=0:N
%          for n3=0:N
%             if (n1+n2+n3==N)
%                B1=lambda_AU*(N-n3);
%                C1=n2;
%                D1=factorial(N)/(factorial(n1)*factorial(n2)*factorial(n3))*(-1)^n1*(-lambda_AU)^n2;
%                for t=0:C1
%                    E1=D1*nchoosek(C1,t)*Gamma_J^t*beta_JU;
%                    
%                    if(B1==0)
%                        P1_7=P1_7+E1*lambda_JU^(-t-m_JU)*gamma(t+m_JU)*gammainc(lambda_JU*tao,t+m_JU,'lower');
%                    else
%                        for m0=0:M
%                           for m00=0:M
%                              for m10=0:M
%                                 for m11=0:M
%                                     if(m0+m00+m10+m11==M)
%                                            B3=lambda_AE*(M-m0);
%                                            C3=m10+m11;
%                                            E00=beta_JE*factorial(m_JE-1);
%                                            E10=lambda_AE*nchoosek(1,0)*beta_JE*factorial(m_JE-1);
%                                            E11=lambda_AE*nchoosek(1,1)*Gamma_J*beta_JE*factorial(1+m_JE-1);
%                                            D3=factorial(M)/(factorial(m0)*factorial(m00)*factorial(m10)*factorial(m11))*(-E00)^m00*(-E10)^m10*(-E11)^m11;
%                                            G3=m00*(m_JE)+m10*m_JE+m11*(1+m_JE);
%                                            if(G3~=0)
%                                               for t1=0:C1
%                                                   H1_1=H(t1+C3,B1*2^Rs+B3,B1*Gamma_J*delta1+lambda_JU,B1*Gamma_J*2^Rs,t+m_JU,lambda_JE,lambda_AE*Gamma_J,G3);
%                                                   P1_1=P1_1+E1*factorial(t+m_JU-1)*nchoosek(C1,t1)*delta1^(C1-t1)*2^(t1*Rs)*exp(-B1*delta1)*D3*B3*H1_1;
%                                                   if(C3~=0)
%                                                       H1_2=H(t1+C3-1,B1*2^Rs+B3,B1*Gamma_J*delta1+lambda_JU,B1*Gamma_J*2^Rs,t+m_JU,lambda_JE,lambda_AE*Gamma_J,G3);
%                                                       P1_2=P1_2+E1*factorial(t+m_JU-1)*nchoosek(C1,t1)*delta1^(C1-t1)*2^(t1*Rs)*exp(-B1*delta1)*D3*C3*H1_2;
%                                                   end
%                                                   
%                                                   H1_3=H(t1+C3,B1*2^Rs+B3,B1*Gamma_J*delta1+lambda_JU,B1*Gamma_J*2^Rs,t+m_JU,lambda_JE,lambda_AE*Gamma_J,G3+1);
%                                                   P1_3=P1_3+E1*factorial(t+m_JU-1)*nchoosek(C1,t1)*delta1^(C1-t1)*2^(t1*Rs)*exp(-B1*delta1)*D3*G3*lambda_AE*Gamma_J*H1_3;
%                                                   for k=0:t+m_JU-1
%                                                      E2=E1*exp(-lambda_JU*tao)*tao^k/factorial(k)*factorial(t+m_JU-1); 
%                                                      H1_4=H(t1+C3,(B1+B1*Gamma_J*tao)*2^Rs+B3,B1*Gamma_J*delta1+lambda_JU,B1*Gamma_J*2^Rs,-k+t+m_JU,lambda_JE,lambda_AE*Gamma_J,G3);
%                                                       P1_4=P1_4+E2*nchoosek(C1,t1)*delta1^(C1-t1)*2^(t1*Rs)*exp(-(B1+B1*Gamma_J*tao)*delta1)*D3*B3*H1_4;
%                                                       if(C3~=0)
%                                                          H1_5=H(t1+C3-1,(B1+B1*Gamma_J*tao)*2^Rs+B3,B1*Gamma_J*delta1+lambda_JU,B1*Gamma_J*2^Rs,-k+t+m_JU,lambda_JE,lambda_AE*Gamma_J,G3);
%                                                       P1_5=P1_5+E2*nchoosek(C1,t1)*delta1^(C1-t1)*2^(t1*Rs)*exp(-(B1+B1*Gamma_J*tao)*delta1)*D3*C3*H1_5; 
%                                                       end
%                                                       H1_6=H(t1+C3,(B1+B1*Gamma_J*tao)*2^Rs+B3,B1*Gamma_J*delta1+lambda_JU,B1*Gamma_J*2^Rs,-k+t+m_JU,lambda_JE,lambda_AE*Gamma_J,G3+1);
%                                                       P1_6=P1_6+E2*nchoosek(C1,t1)*delta1^(C1-t1)*2^(t1*Rs)*exp(-(B1+B1*Gamma_J*tao)*delta1)*D3*G3*lambda_AE*Gamma_J*H1_6;
%                                                    end
%                                                   
%                                               end
%                                            end 
%                                     end
%                                 end
%                              end
%                           end
%                        end 
%                    end
%                end
% 
%             end
%          end
%      end
%   end
%  P1=-P1_1+P1_2-P1_3+P1_4-P1_5+P1_6+P1_7;
%  Pout_AJ_t(kk)=P1+P2;
syms v x;
Fv=0;
for k=0:m_AE-1
   for t=0:k
      E_kt=lambda_AE^k/factorial(k)*nchoosek(k,t)*Gamma_J^t*beta_JE*factorial(t+m_JE-1);
      Fv=Fv+E_kt*v^k*exp(-lambda_AE*v)*(lambda_JE+lambda_AE*Gamma_J*v)^(-t-m_JE);
   end
end
Fv=(1-Fv)^M;
fv=diff(Fv,v);
Fw=0;
fx=beta_JU*exp(-lambda_JU*x)*x^(m_JU-1);
 for n1=0:N 
     for n2=0:N
         for n3=0:N
            if (n1+n2+n3==N)
               B1=lambda_AU*(N-n3);
               C1=n2;
               D1=factorial(N)/(factorial(n1)*factorial(n2)*factorial(n3))*(-1)^n1*(-lambda_AU)^n2;
               

                Fw=Fw+D1*(delta1*(x*Gamma_J+1)+2^Rs*v*(x*Gamma_J+1))^C1*exp(-B1*(delta1*(x*Gamma_J+1)+2^Rs*v*(x*Gamma_J+1)));
            end
         end
     end
 end
  fun=eval(['@(x,v)',vectorize(Fw*fv*fx)]);
 P1=double(integral2(fun,0,tao,0,inf));
 Pout_AJ_t(kk)=P1+P2;


%%%asy

 %% NJ
 Pout_NJ_s_asy(kk)=0;
 for i=1:L
     [max_A,index_A]=max(h_bu(:,i));
     h_AE=h_be(:,i);
     [max_E,index_E]=max(h_AE);
     if(max_A<2^Rs*max_E)
         Pout_NJ_s_asy(kk)=Pout_NJ_s_asy(kk)+1/L;
     end
 end

 Pout_NJ_t_asy(kk)=0;

 for n1=0:N 
     for n2=0:N
         for n3=0:N
            if (n1+n2+n3==N)
               B1=lambda_AU*(N-n3);
               C1=n2;
               D1=factorial(N)/(factorial(n1)*factorial(n2)*factorial(n3))*(-1)^n1*(-lambda_AU)^n2;
               for m1=0:M
                  for m2=0:M
                     for m3=0:M
                        if(m1+m2+m3==M)
                           B2=lambda_AE*(M-m3);
                           C2=m2;
                           D2=factorial(M)/(factorial(m1)*factorial(m2)*factorial(m3))*(-1)^m1*(-lambda_AE)^m2;
                           temp0=D1*2^(C1*Rs);
                           temp1=0;
                           temp2=0;
                           if(B2~=0)
                               temp1=D2*B2*factorial(C1+C2)*(B1*2^Rs+B2)^(-C1-C2-1);
                           end
                         if(C2~=0)
                            temp2=D2*C2*factorial(C1+C2-1)*(B1*2^Rs+B2)^(-C1-C2);
                         end
                         Pout_NJ_t_asy(kk)=Pout_NJ_t_asy(kk)+temp0*(temp2-temp1);
                        end
                     end
                  end
               end
            end
         end
     end
 end
     
 
 
 %% AJ
 
  %%%AJ
 
 
   Pout_AJ_s_asy(kk)=0;
   for i=1:L
       [max_A,index_A]=max(h_bu(:,i));
       h_AE=h_be(:,i);
     [max_E,index_E]=max(h_AE);
    C_bu_aj=log2(1+max_A/(u*h_ju(i)));
   C_be_aj=log2(1+h_be(:,i)./(u*h_je(:,i)));
        [max_E2,index_E2]=max(C_be_aj);
      if(h_ju(i)<tao)
          if(C_bu_aj-max_E2<Rs)
             Pout_AJ_s_asy(kk)=Pout_AJ_s_asy(kk)+1/L;
          end
      else
          if(max_A<2^Rs*max_E)
              Pout_AJ_s_asy(kk)=Pout_AJ_s_asy(kk)+1/L;
          end
      end
   end

   
   
   P2_asy=Pout_NJ_t_asy(kk)*P2_1;
   
%  P1_1_asy=0;
%  P1_2_asy=0;
%  P1_3_asy=0;
%  P1_4_asy=0;
%  P1_5_asy=0;
%   for n1=0:N 
%      for n2=0:N
%          for n3=0:N
%             if (n1+n2+n3==N)
%                B1=lambda_AU*(N-n3);
%                C1=n2;
%                D1=factorial(N)/(factorial(n1)*factorial(n2)*factorial(n3))*(-1)^n1*(-lambda_AU)^n2;
%                E3=D1*beta_JU*factorial(C1+m_JU-1);
% 
%                if(B1==0)
%                   P1_5_asy=P1_5_asy+D1*beta_JU*lambda_JU^(-m_JU)*gamma(m_JU)*gammainc(lambda_JU*tao,m_JU,'lower');
%                else
%                    for m1=0:M
%                       for m2=0:M
%                          for m3=0:M
%                             if(m1+m2+m3==M)
%                                C4=m2;
%                                G4=m_JE*m1+(1+m_JE)*m2;
%                                D4=factorial(M)/(factorial(m1)*factorial(m2)*factorial(m3))*(-beta_JE*factorial(m_JE-1))^m1*(-beta_JE*lambda_AE*factorial(2-1+m_JE-1))^m2;
%                                if(G4~=0)
%                                   for t=0:C1
%                                     if(C4~=0)
%                                         Q1=Q(t+C4-1,B1*delta2+lambda_JU,B1*2^Rs,C1+m_JU,lambda_JE,lambda_AE,G4);
%                                         P1_1_asy=P1_1_asy+E3*nchoosek(C1,t)*delta2^(C1-t)*2^(t*Rs)*D4*C4*Q1;
%                                     end
%                                     Q2=Q(t+C4,B1*delta2+lambda_JU,B1*2^Rs,C1+m_JU,lambda_JE,lambda_AE,G4+1);
%                                     P1_2_asy=P1_2_asy+E3*nchoosek(C1,t)*delta2^(C1-t)*2^(t*Rs)*D4*G4*lambda_AE*Q2;
%                                     for k=0:C1+m_JU-1
%                                          E4=D1*beta_JU* factorial(C1+m_JU-1)*exp(-lambda_JU*tao)*tao^k/factorial(k);
%                                          if(C4~=0)
%                                             H1= H(t+C4-1,B1*tao*2^Rs,B1*delta2+lambda_JU,B1*2^Rs,-k+C1+m_JU,lambda_JE,lambda_AE,G4);
%                                             P1_3_asy=P1_3_asy+E4*nchoosek(C1,t)*delta2^(C1-t)*2^(t*Rs)*exp(-B1*tao*delta2)*D4*C4*H1;
%                                          end
%                                          H2= H(t+C4,B1*tao*2^Rs,B1*delta2+lambda_JU,B1*2^Rs,-k+C1+m_JU,lambda_JE,lambda_AE,G4+1);
%                                          P1_4_asy=P1_4_asy+E4*nchoosek(C1,t)*delta2^(C1-t)*2^(t*Rs)*exp(-B1*tao*delta2)*D4*G4*lambda_AE*H2;
%                                     end
%                                   end
%                                end 
%                                
%                             end
%                          end
%                       end
%                    end
% 
%                end
%                
%             end
%          end
%      end
%  end
% P1_asy=P1_1_asy-P1_2_asy-P1_3_asy+P1_4_asy+P1_5_asy;
% Pout_AJ_t_asy(kk)=P1_asy+P2_asy;

 %%%
P1_asy=0;
               syms v x;
               Fv=0;
               for k=0:m_AE-1
                Fv=Fv+lambda_AE^k/factorial(k)*beta_JE*factorial(k+m_JE-1)*v^k*(lambda_AE*v+lambda_JE)^(-k-m_JE);
               end
               Fv2=(1-Fv)^M;
               fv=diff(Fv2,v);
               fx=beta_JU*exp(-lambda_JU*x)*x^(m_JU-1);
               Fw=0;
 for n1=0:N 
     for n2=0:N
         for n3=0:N
            if (n1+n2+n3==N)
               B1=lambda_AU*(N-n3);
               C1=n2;
               D1=factorial(N)/(factorial(n1)*factorial(n2)*factorial(n3))*(-1)^n1*(-lambda_AU)^n2;

               
               Fw=Fw+D1*(delta2*x+2^Rs*v*x)^C1*exp(-B1*(delta2*x+2^Rs*v*x));
            end
         end
     end
 end
 fun=eval(['@(x,v)',vectorize(Fw*fv*fx)]);
 P1_asy=double(integral2(fun,0,tao,0,inf));
 Pout_AJ_t_asy(kk)=P1_asy+P2_asy;

 
%% IC
LL=30;
J1=0;
J2=0;
for l=1:LL
   sita_l=cos((2*l-1)*pi/(2*LL));
   w_l=(sita_l+1)*pi/4;
    v_l=tan(w_l);
    fx=0;
    for n1=0:N 
     for n2=0:N
         for n3=0:N
            if (n1+n2+n3==N)
               B1=lambda_AU*(N-n3);
               C1=n2;
               D1=factorial(N)/(factorial(n1)*factorial(n2)*factorial(n3))*(-1)^n1*(-lambda_AU)^n2;
               if(B1~=0)
                   fx=fx-D1*B1*tan(w_l)^C1*exp(-B1*tan(w_l));
               end
               if(C1~=0)
                   fx=fx+D1*C1*tan(w_l)^(C1-1)*exp(-B1*tan(w_l));
               end
            end
         end
     end
    end
    for k=0:m_JU-1
        J1=J1+lambda_JU^k/factorial(k)*pi^2/4/LL*sqrt(1-sita_l^2)*(sec(w_l))^2*(afa1/tan(w_l))^k*exp(-lambda_JU*afa1/tan(w_l))*fx;
    end
    
    for m1=0:M
       for m2=0:M
          if(m1+m2==M)
             C5=m2;
             G5=m1*(m_AE)+m2*(1+m_AE);
             D5=factorial(M)/(factorial(m1)*factorial(m2))*(factorial(m_AE-1)*beta_AE)^m1*(lambda_JE*beta_AE*factorial(2-1+m_AE-1))^m2;
             E5=(-1/tan(w_l))^C5*(-lambda_JE/tan(w_l))^(-G5);
             afax=lambda_AE*tan(w_l)/lambda_JE+afa1/tan(w_l);
             for t1=0:C5
                for t2=0:m_JU-1
                   T1=E5*nchoosek(C5,t1)*(lambda_AE/lambda_JE*tan(w_l))^(C5-t1)*beta_JU*nchoosek(m_JU-1,t2)*afax^(m_JU-1-t2)*exp(-lambda_JU*afax);
                   h1=t1-G5+t2;
                   
                   J2_2=0;
                   for i=0:30
                       if(i==-1-h1)
                           J2_2=J2_2+(-lambda_JU)^(-1-h1)/factorial(-1-h1)*(log(abs(-lambda_AE*tan(w_l)/lambda_JE))-log(abs(-afax)));
                       else
                       J2_2=J2_2+(-lambda_JU)^i/factorial(i)/(i+h1+1)*((-lambda_AE*tan(w_l)/lambda_JE)^(i+h1+1)-(-afax)^(i+h1+1));
                       end
                   end
                   faix=T1*(J2_2);
                   J2=J2+D5*pi^2/4/LL*sqrt(1-sita_l^2)*(sec(w_l))^2*faix*fx;
                end
             end
          end
       end
    end
    
end
J1=1-J1;
Pout_IC_t(kk)=J1-J2;
    clear h_ud h_ue v;
     kk=kk+1;
end
delete(h);

SNR=0:4:32

plot(SNR,Pout_NJ_t,'k-o');
 hold on;
  plot(SNR,Pout_AJ_t,'r-d')
 hold on;
plot(SNR,Pout_IC_s,'b-p')
hold on;
%  
 plot(SNR,Pout_NJ_t_asy,'r--')
 hold on;
 plot(SNR,Pout_NJ_s,'k*');
hold on;
 plot(SNR,Pout_NJ_s_asy,'k*')
 hold on;
plot(SNR,Pout_AJ_s,'k*')
hold on;
plot(SNR,Pout_AJ_s_asy,'k*')
hold on;
% 

plot(SNR,Pout_AJ_t_asy,'r--')
hold on;

plot(SNR,Pout_IC_s_high,'k*')
hold on;
plot(SNR,Pout_IC_t,'r--')
hold on;
% 
% 

legend('NJ scheme (t.)','AJ scheme (t.)','IC scheme (s., Eq. (35))','asymptotic','s.');
xlabel('SNR $\mathop \gamma \nolimits_B $(dB)','interpreter','latex');    
ylabel('Secrecy outage probability');
      set(gca,'yscale','log');
  axis([0 32 1e-5 1]);