%3D random media with a single layer of homogeneous media. This code is an
%iterative approach for the single layer problem. After this step is done,
%we can decide either to go with the gradually change profile or
%multi-layer eigen method appraoch. 
eps0=8.854e-12;
miu0=4*pi*1e-7;
c=1/sqrt(eps0*miu0);

fGhz=0.5;
k0=fGhz*1e9*2*pi*sqrt(eps0*miu0);
lambda=2*pi/k0;

Niter=180;

kb=1.38;

d=12.3;
dz=0.01;
z_a=0:-dz:-d;
% z_a=z_a;
lza=length(z_a);
% T0=265;
z_s=[0;-0.45;-0.95;-1.54;-2;-2.5;-2.9;-4;-5;-6;-6.9;-8;-9;-10;-11.5;-12.3];
T_s=[-19;-16;-12.8;-11;-9.4;-8.67;-7.7;-6.22;-5.06;-4.13;-3.2;-2.44;-1.69;-0.81;-0.35;0];
T0=interp1(z_s,T_s,z_a)+273;
T2=273;




% epsrm=(1.64+0.001i)*ones(lza,1);
% delta1_a=0.002*ones(lza,1);%0.008;
% corl_z_a=0.02*ones(lza,1);
% corl_rho_a=0.5*ones(lza,1);


rhom_a=0.922-0.5761*exp(0.0941*z_a);
delta_rho_a=0.088*ones(lza,1);
corl_z_a=0.085*ones(lza,1);
corl_rho_a=0.21*ones(lza,1);
epsrm=zeros(lza,1);
delta1_a=zeros(lza,1);
for irho=1:1:lza
    delta_rho=delta_rho_a(irho);
    rhom=rhom_a(irho);
    T=T0(irho);
    [epsrm_MT,delta_Matzler1]=Delta_calculator_MT_func(delta_rho,rhom,fGhz,T);
    epsrm(irho)=epsrm_MT;
    delta1_a(irho)=delta_Matzler1;
end



k1mp_a=real(k0*sqrt(epsrm));
kappa_a=2*imag(k0*sqrt(epsrm));
epsrm0=epsrm(1);
epsrmd=epsrm(end);
epsr_aquifer=7.6+0.25i;%7.6+0.25i;%9.5+0.42i;%%4.8+0.05i;%7.6+0.25i;%9.5+1i;%4.1;%7.6+0.25i%8.4+1i;
epsr0=1;



C=1;

Nquad=64;
Ndeg=Nquad/2;
[x,wi]=GLNodeWt(Nquad);
miu0=-x(1:Ndeg);
miu=[-x(1:Ndeg);x(1:Ndeg)];
miu_up=[-x(1:Ndeg);-x(1:Ndeg)];
miu_up_inv=1./miu_up;
Miu_up_inv=diag(miu_up_inv);

a=wi(1:Ndeg);
wii=[wi(1:Ndeg);wi(1:Ndeg)];
Wm=diag(wii);%weighting value matrix a
theta_a=acos(miu0);
theta_ob=theta_a;
T2_vec=T2*ones(Nquad,1);
% T0_vec=T0*ones(1,lza);
T0_vec=T0;
% T0_m=T0*ones(Nquad,lza);
T0_m=repmat(T0,Nquad,1);
%reflectivity for the 2 boundary 
rv10_Ndeg=zeros(Ndeg,1);
rh10_Ndeg=zeros(Ndeg,1);
rv12_Ndeg=zeros(Ndeg,1);
rh12_Ndeg=zeros(Ndeg,1);

for ideg=1:1:Ndeg
    tai=theta_a(ideg);
    [rv10,rh10,~,~]=Fresnel(epsrm0,epsr0,tai);%reflectivity at snow air interface
    [rv12,rh12,~,~]=Fresnel(epsrmd,epsr_aquifer,tai);
    rv10_Ndeg(ideg)=rv10;
    rh10_Ndeg(ideg)=rh10;
    rv12_Ndeg(ideg)=rv12;
    rh12_Ndeg(ideg)=rh12;
end
w_crt=100;
%obtain the scattering coefficient and the phase matrix for each depth
F_dp=zeros(Nquad,Nquad,lza);
B_dp=zeros(Nquad,Nquad,lza);
kappa_sv_dbp=zeros(Ndeg,lza);
kappa_sh_dbp=zeros(Ndeg,lza);
kappa_ev_dbp=zeros(Ndeg,lza);
fprintf('Start formulating phase matrix....')
kappa_a_a=kappa_a.';
kappa_a_m=repmat(kappa_a_a,Nquad/2,1);
kappa_a_m2=repmat(kappa_a_a,Nquad,1);
tic
for izp=1:1:lza
    klmp=k1mp_a(izp);
    delta1=delta1_a(izp);
    corl_z=corl_z_a(izp);
    corl_rho=corl_rho_a(izp);
    [F,B,kappa_sv_a, kappa_sh_a]=phase_matrix_3D_rand_med(klmp,miu,wii,delta1,corl_z,corl_rho,Nquad,w_crt);
    F_dp(:,:,izp)=F;
    B_dp(:,:,izp)=B;
    kappa_sv_dbp(:,izp)=kappa_sv_a;
    kappa_sh_dbp(:,izp)=kappa_sh_a;
    kappa_ev_dbp=kappa_sv_dbp+kappa_a_m;
    kappa_eh_dbp=kappa_sh_dbp+kappa_a_m;
    if mod(izp,100)==0
        fprintf('iz=%d in total %d\n',izp,lza);
    end
end
toc
%Here we write down the 0th order solution

%First step is to obatin the factor that accounts for the double bounce
int_kappa_ev=trapz(kappa_ev_dbp,2)*dz;%sum(kappa_ev_dbp(:,2:lza),2)*dz;%does not use the first point for inetral
int_kappa_eh=trapz(kappa_eh_dbp,2)*dz;
int_kappa_e_vec=[int_kappa_ev;int_kappa_eh];


%here implement the 0th order solution
R12_vec=[rv12_Ndeg;rh12_Ndeg];
R10_vec=[rv10_Ndeg;rh10_Ndeg];
R12=diag(R12_vec);
R10=diag(R10_vec);
beta_f2=exp(-2*miu_up_inv.*int_kappa_e_vec);
beta_f=exp(-miu_up_inv.*int_kappa_e_vec);
Beta_f=diag(beta_f);
zeta=1-beta_f2.*R12_vec.*R10_vec;
Zeta_inv=diag(1./zeta);

Iu0=zeros(Nquad,lza);
Id0=zeros(Nquad,lza);
uniI=eye(Nquad,Nquad);
fprintf('Start preparing 0th order soltuion...\n')
term2_u=zeros(Nquad,lza);
%0th order solution at z=-d
Tdd_a=zeros(Nquad,lza);
for iz3=1:1:lza
    kappa_a=kappa_a_a(iz3);
    T0=T0_vec(iz3);
    indz_u32=iz3:1:lza;
    ke_term3_u2=[kappa_ev_dbp(:,indz_u32);kappa_eh_dbp(:,indz_u32)];
    ke_term3_u2_int=trapz(ke_term3_u2,2)*dz;
%     Ke_term3_u2_int=diag(ke_term3_u2_int);
%     Beta_zd32=exp(-Miu_up_inv*Ke_term3_u2_int);
    beta_zd32=exp(-miu_up_inv.*ke_term3_u2_int);
    Beta_zd32=diag(beta_zd32);
    Tdd_a(:,iz3)=Beta_zd32*kappa_a*T0*ones(Nquad,1);
end
Tdd=Miu_up_inv*trapz(Tdd_a,2)*dz;
%0th order solution for upward emission at z=0
Tuu_a=zeros(Nquad,lza);
for iz4=1:1:lza
    kappa_a=kappa_a_a(iz4);
    T0=T0_vec(iz4);
    indz_u42=1:1:iz4;
    ke_term4_u2=[kappa_ev_dbp(:,indz_u42);kappa_eh_dbp(:,indz_u42)];
    ke_term4_u2_int=trapz(ke_term4_u2,2)*dz;
%     Ke_term4_u2_int=diag(ke_term4_u2_int);
%     Beta_zd42=exp(-Miu_up_inv*Ke_term4_u2_int);
    beta_zd42=exp(-miu_up_inv.*ke_term4_u2_int);
    Beta_zd42=diag(beta_zd42);
    Tuu_a(:,iz4)=Beta_zd42*kappa_a*T0*ones(Nquad,1);
end
Tuu=Miu_up_inv*trapz(Tuu_a,2)*dz;

%calculate the vector of beta_zd and beta_0z for different values of z
%before starting all the calculations 
beta_zd_a=zeros(Nquad,lza);
beta_0z_a=zeros(Nquad,lza);
for izb=1:1:lza
    ke_term1_u=[kappa_ev_dbp(:,izb:lza);kappa_eh_dbp(:,izb:lza)];
    ke_term1_u_int=trapz(ke_term1_u,2)*dz;
    beta_zd_tp=exp(-miu_up_inv.*ke_term1_u_int);
    beta_zd_a(:,izb)=beta_zd_tp;
    
    ke_term1_d=[kappa_ev_dbp(:,1:izb);kappa_eh_dbp(:,1:izb)];
    ke_term1_d_int=trapz(ke_term1_d,2)*dz;
    beta_0z_tp=exp(-miu_up_inv.*ke_term1_d_int);
    beta_0z_a(:,izb)=beta_0z_tp;
end

for iz=1:1:lza
    if mod(iz,100)==0
        fprintf('iz=%d in total %d\n',iz,lza);
    end
    %implement the upward going intensity
    %term1 direction emission from lower layer.

    beta_zd=beta_zd_a(:,iz);
    Beta_zd=diag(beta_zd);
    term1_u=Zeta_inv*(uniI-R12)*Beta_zd*T2_vec;
    %term2 direct emission
    term2_u_a=zeros(Nquad,lza-iz+1);
    %use cumtrapz to replace the double integral 
    
     indz_u2=iz:lza;
     ke_term2_u=[kappa_ev_dbp(:,indz_u2);kappa_eh_dbp(:,indz_u2)];
     ke_term2_u_int_a=cumtrapz(ke_term2_u,2)*dz;
     beta_zzp_m=exp(-Miu_up_inv*ke_term2_u_int_a);
     kappa_a_mz=kappa_a_m2(:,iz:lza);
     T0_mz=T0_m(:,iz:lza);
     term2_u_cz=Miu_up_inv*trapz(beta_zzp_m.*kappa_a_mz.*T0_mz,2)*dz;
     
%      tic
%     for iiz=iz:1:lza
%         indz_u2=iz:iiz;
%         ke_term2_u=[kappa_ev_dbp(:,indz_u2);kappa_eh_dbp(:,indz_u2)];
%         ke_term2_u_int=trapz(ke_term2_u,2)*dz;
% %         Ke_term2_u_int=diag(ke_term2_u_int);
% %         Beta_zzp=exp(-Miu_up_inv*Ke_term2_u_int);
%         beta_zzp=exp(-miu_up_inv.*ke_term2_u_int);
%         Beta_zzp=diag(beta_zzp);
%         kappa_a_iiz=kappa_a_a(iiz);
%         T0_iiz=T0_vec(iiz);
%         term2_u_a(:,iiz-iz+1)=Beta_zzp*kappa_a_iiz*T0_iiz*ones(Nquad,1);
%     end
%     term2_u=Miu_up_inv*trapz(term2_u_a,2)*dz;
%     toc
    %term3 emission due to lower boundary reflection

  
    Beta_zd3=Beta_zd;

    %total upward emission
    term3_u=Zeta_inv*Beta_zd3*R12*Tdd;
    %term4 emission due to upper and lower boundary double reflection
    term4_u=Zeta_inv*Beta_zd3*R10*R12*Beta_f*Tuu;
    
    
    %implement the downward going specific intenstiy
    %term1 direct emission then reflected by the upper boundary
    beta_0zd1=beta_0z_a(:,iz);
    Beta_0zd1=diag(beta_0zd1);
%     Ke_term1_d_int=diag(ke_term1_d_int);
%     Beta_0z=exp(-Miu_up_inv*Ke_term1_d_int);
    term1_d=Zeta_inv*R10*Beta_0zd1*(uniI-R12)*Beta_f*T2_vec;
    %term2 emission from layer
    
     indz_d2=1:iz;
     ke_term2_d=[kappa_ev_dbp(:,indz_d2);kappa_eh_dbp(:,indz_d2)];
     ke_term2_d=fliplr(ke_term2_d);
     ke_term2_d_int_a=cumtrapz(ke_term2_d,2)*dz;
     ke_term2_d_int_a=fliplr(ke_term2_d_int_a);
     beta_zpz_m=exp(-Miu_up_inv*ke_term2_d_int_a);
     kappa_a_mz=kappa_a_m2(:,1:iz);
     T0_mz=T0_m(:,1:iz);
     term2_d_cz=Miu_up_inv*trapz(beta_zpz_m.*kappa_a_mz.*T0_mz,2)*dz;
    
%     term2_d_a=zeros(Nquad,iz);
%     for izd2=1:1:iz
%         indzd2=izd2:iz;
%         ke_term2_d=[kappa_ev_dbp(:,indzd2);kappa_eh_dbp(:,indzd2)];
%         ke_term2_d_int=trapz(ke_term2_d,2)*dz;
% %         Ke_term2_d_int=diag(ke_term2_d_int);
% %         Beta_zpz=exp(-Miu_up_inv*Ke_term2_d_int);
%         beta_zpz=exp(-miu_up_inv.*ke_term2_d_int);
%         Beta_zpz=diag(beta_zpz);
%         kappa_a_izd2=kappa_a_a(izd2);
%         T0_izd2=T0_vec(izd2);
%         term2_d_a(:,izd2)=Beta_zpz*T0_izd2*kappa_a_izd2*ones(Nquad,1);
%     end
%     term2_d=Miu_up_inv*trapz(term2_d_a,2)*dz;
    %term3 reflection from the top boundary

    beta_0z=beta_0z_a(:,iz);
    Beta_0zd=diag(beta_0z);
    term3_d=Zeta_inv*R10*Beta_0zd*Tuu;
    %term4 double reflection from downward emission
    term4_d=Zeta_inv*R10*R12*Beta_0zd*Beta_f*Tdd;
    
%     Iu0(:,iz)=term1_u+term2_u+term3_u+term4_u;
%     Id0(:,iz)=term1_d+term2_d+term3_d+term4_d;
    Iu0(:,iz)=term1_u+term2_u_cz+term3_u+term4_u;
    Id0(:,iz)=term1_d+term2_d_cz+term3_d+term4_d;
    
end

%start performing the iterations
toc
Iun_a=zeros(Nquad,lza,Niter);
Idn_a=zeros(Nquad,lza,Niter);

fprintf('Start iteration.....\n')
for iter=1:1:Niter
    if iter==1
       Ium=Iu0;
       Idm=Id0;
    else
        Ium=Iun;
        Idm=Idn;
    end
%intermediate results from multiplication of specific intensity and phase
%matrix
    Iumn=zeros(Nquad,lza);
    Idmn=zeros(Nquad,lza);

%
% Iun=zeros(Nquad,lza);
% Idn=zeros(Nquad,lza);

%calculate the intensity after multiplying with phase matrix
    for izz=1:1:lza
        F=F_dp(:,:,izz);
        B=B_dp(:,:,izz);
        Iumz=Ium(:,izz);
        Idmz=Idm(:,izz);
        Iumn(:,izz)=F*Wm*Iumz+B*Wm*Idmz;
        Idmn(:,izz)=B*Wm*Iumz+F*Wm*Idmz;
    end

%calculate the upward going specific intensity
    S=zeros(Nquad,lza);
    W=zeros(Nquad,lza);
    fprintf('Working on direct scattering...\n')
    for izs=1:1:lza
        if mod(izs,100)==0
        fprintf('iz=%d in total %d\n',izs,lza);
        end
        %use trapzoidal rule to 
%         Szp=zeros(Nquad,lza-izs+1);
%         Wzp=zeros(Nquad,izs);
           indz_s=izs:lza;
           ke_term_s=[kappa_ev_dbp(:,indz_s);kappa_eh_dbp(:,indz_s)];
           ke_term2_s_int_a=cumtrapz(ke_term_s,2)*dz;
           beta_zzp_m=exp(-Miu_up_inv*ke_term2_s_int_a);
           Iumn_m=Iumn(:,indz_s);
           S(:,izs)=Miu_up_inv*trapz(beta_zzp_m.*Iumn_m,2)*dz;
%         for izs2=izs:lza
%             inds2=izs:1:izs2;
%             kappa_e_s=[kappa_ev_dbp(:,inds2);kappa_eh_dbp(:,inds2)];
%             kappa_e_s_int=trapz(kappa_e_s,2)*dz;
%             beta_zzp=exp(-miu_up_inv.*kappa_e_s_int);
%             Beta_zzp=diag(beta_zzp);
%             Kappa_e_s_int=diag(kappa_e_s_int);
%             Beta_zzp=exp(-Miu_up_inv*Kappa_es_int);
%             Szp(:,izs2-izs+1)=Beta_zzp*Iumn(:,izs2);
%         end
%         S(:,izs)=Miu_up_inv*trapz(Szp,2)*dz;
         indw2=1:izs;
         ke_term_w=[kappa_ev_dbp(:,indw2);kappa_eh_dbp(:,indw2)];
         ke_term_w=fliplr(ke_term_w);
         ke_term2_w_int_a=cumtrapz(ke_term_w,2)*dz;
         ke_term2_w_int_a=fliplr(ke_term2_w_int_a);
         beta_zpz_m=exp(-Miu_up_inv*ke_term2_w_int_a);
         Idmn_m=Idmn(:,indw2);
         W(:,izs)=Miu_up_inv*trapz(beta_zpz_m.*Idmn_m,2)*dz;
        
    end
    
    % boundary reflections due to the scattering term
    fprintf('Working on boundary refelction terms....\n')
    Wzd=W(:,lza);
    Sz0=S(:,1);

    Wu=zeros(Nquad,lza);
    Sdu=zeros(Nquad,lza);
    Sd=zeros(Nquad,lza);
    Wud=zeros(Nquad,lza);
    for izr=1:1:lza
%         indrs=izr:lza;
%         kappa_e_sdb=[kappa_ev_dbp(:,indrs);kappa_eh_dbp(:,indrs)];
%         kappa_e_sdb_int=trapz(kappa_e_sdb,2)*dz;
%         beta_zd=exp(-miu_up_inv.*kappa_e_sdb_int);
        beta_zd=beta_zd_a(:,izr);
        Beta_zd=diag(beta_zd);
        Wu(:,izr)=Zeta_inv*R12*Beta_zd*Wzd;
        Sdu(:,izr)=Zeta_inv*R10*R12*Beta_zd*Beta_f*Sz0;
        

        beta_0z=beta_0z_a(:,izr);
        Beta_0z=diag(beta_0z);
        Sd(:,izr)=Zeta_inv*R10*Beta_0z*Sz0;
        Wud(:,izr)=Zeta_inv*R10*R12*Beta_0z*Beta_f*Wzd;      
    end
    Iun=S+Wu+Sdu;
    Idn=W+Sd+Wud;


    Iun_a(:,:,iter)=Iun;
    Idn_a(:,:,iter)=Idn;



    fprintf('Done with %d order iteration total order %d\n',iter,Niter)
end
Ius=zeros(Nquad,lza);
Ids=zeros(Nquad,lza);

for iorder=1:1:Niter
    Ius=Ius+Iun_a(:,:,iorder);
    Ids=Ids+Idn_a(:,:,iorder);
end
Iu=Iu0+Ius;
Id=Id0+Ids;

Iu_all_v=Iu(1:Ndeg,1);
Iu_all_h=Iu(Ndeg+1:Nquad,1);
tai_crt=asin(1/real(sqrt(epsrm0)));
tai_id=theta_a<tai_crt;
klmp0=k1mp_a(1);
degs0=asin(klmp0*sin(theta_a(tai_id))/k0)/pi*180;
TB_v=(1-rv10_Ndeg(tai_id)).*Iu_all_v(tai_id);
TB_h=(1-rh10_Ndeg(tai_id)).*Iu_all_h(tai_id);
figure;plot(degs0,TB_v,degs0,TB_h)
figure;plot(degs0,Iu_all_v(tai_id),degs0,Iu_all_h(tai_id))