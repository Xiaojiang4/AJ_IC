clear all;
clc;


L=200000;

 m_AU=2;
 omiga_AU=1;
 lambda_AU=m_AU/omiga_AU;
 beta_AU=lambda_AU^m_AU/gamma(m_AU);
 
 m_AE=2;
 omiga_AE=1;
 lambda_AE=m_AE/omiga_AE;
 beta_AE=lambda_AE^m_AE/gamma(m_AE);
 
 m_JU=2;
 omiga_JU=0.5;
 lambda_JU=m_JU/omiga_JU;
 beta_JU=lambda_JU^m_JU/gamma(m_JU);
 
 m_JE=2;
 omiga_JE=0.5;
 lambda_JE=m_JE/omiga_JE;
 beta_JE=lambda_JE^m_JE/gamma(m_JE);
 
 
 tao=0.5;
 Rs=1;
 h=waitbar(0,'please wait');
kk=1;
for SNR=-5:2:15
    str=['运行中...',num2str(kk/10*100),'%'];
    waitbar(kk/10,h,str)
Gamma_B=10^(SNR/10);

u=0.5;
 M=2;
N=2;
   Gamma_J=Gamma_B*u;

     h_bu = gamrnd(m_AU,omiga_AU/m_AU,N,L);
    h_be = gamrnd(m_AE,omiga_AE/m_AE,M,L);   
    h_ju = gamrnd(m_JU,omiga_JU/m_JU,1,L);
    h_je = gamrnd(m_JE,omiga_JE/m_JE,M,L);
    delta1=(2^Rs-1)/Gamma_B;
    
    afa1=2^Rs*omiga_AU/Gamma_J;
    delta2=(2^Rs-1)*u;
   
 %%%simulation
 %%%NJ
%  
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

%%% 
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


    clear h_ud h_ue v;
     kk=kk+1;
end
delete(h);


SNR=-5:2:15

 plot(SNR,Pout_NJ_t,'k-o');
 hold on;
 plot(SNR,Pout_AJ_t,'r-d')
 hold on;
  plot(SNR,Pout_IC_s,'b-p')
 hold on;
% % %  



 plot(SNR,Pout_AJ_s,'k+')
 hold on;


 plot(SNR,Pout_NJ_s,'k+');
hold on;



legend('NJ scheme (t.)','AJ scheme (t.)','IC scheme (s., Eq. (35))','s.');
xlabel('SNR $\mathop \gamma \nolimits_B $(dB)','interpreter','latex');    
ylabel('Secrecy outage probability');
      set(gca,'yscale','log');
     axis([-5 15 1e-2 1]);