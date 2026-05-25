% Verify all numerical answers for ET4260 Unit 4 homework
clear; clc; format short g
e = 1.602e-19; kB = 1.381e-23; eps0 = 8.854e-12; T = 300;

fprintf('===== PART A: Hall =====\n');
SI2 = 90; SV2 = 0.057; Rsq2 = 430;
LW2 = SI2/(SV2*Rsq2);
G2  = LW2^2/sqrt(LW2^4 + 0.5*LW2^2 + 4);
fprintf('Gen2: (L/W)eff=%.3f  G=%.4f\n', LW2, G2);
SI3 = 34; SV3 = 0.034;
Rsq3 = SI3/(SV3*LW2);
fprintf('Gen3: Rsq=%.1f ohm/sq\n', Rsq3);
ns2 = G2/(SI2*e);  mu2 = SV2*LW2/G2;
ns3 = G2/(SI3*e);  mu3 = SV3*LW2/G2;
fprintf('Gen2: n*t=%.3e cm^-2  mu=%.0f cm2/Vs\n', ns2*1e-4, mu2*1e4);
fprintf('Gen3: n*t=%.3e cm^-2  mu=%.0f cm2/Vs\n', ns3*1e-4, mu3*1e4);
fprintf('  check mu2 via Rsq: %.0f cm2/Vs\n', 1/(Rsq2*e*ns2)*1e4);
fprintf('  check mu3 via Rsq: %.0f cm2/Vs\n', 1/(Rsq3*e*ns3)*1e4);
SIsi = 300; SVsi = 0.060; LWsi = sqrt(2);
Gsi = LWsi^2/sqrt(LWsi^4 + 0.5*LWsi^2 + 4);
Rsqsi = SIsi/(SVsi*LWsi);
nssi = Gsi/(SIsi*e); musi = SVsi*LWsi/Gsi;
fprintf('Si: G=%.4f Rsq=%.0f ohm/sq  n*t=%.3e cm^-2  mu=%.0f cm2/Vs\n', Gsi, Rsqsi, nssi*1e-4, musi*1e4);

fprintf('\n===== PART C: Photodetectors =====\n');
a600=4140; a400=95200;
fprintf('600nm: d90=%.2f um d99=%.2f um\n', log(10)/a600*1e4, log(100)/a600*1e4);
fprintf('400nm: d90=%.0f nm d99=%.0f nm\n', log(10)/a400*1e7, log(100)/a400*1e7);
fprintf('Rel change QE 90->99%%: %.1f%%\n', (0.99-0.90)/0.90*100);
fprintf('eta(400nm,99%%)=%.3f  eta(600nm,99%%)=%.3f\n', (1-0.485)*0.99, (1-0.354)*0.99);
ND=1e17; mun=700; mup=240; tcm=1e-4; Wcm=10e-4; Lcm=30e-4; V=1.8;
sigma = ND*e*mun; rho = 1/sigma; R = rho*Lcm/(Wcm*tcm); Idark = V/R;
fprintf('sigma=%.3f S/cm R=%.0f ohm Idark=%.3f mA\n', sigma, R, Idark*1e3);
phi=1e16; Lm=30e-6; Wm=10e-6; Phi=phi*Lm*Wm; eta=0.515; tau=100e-9;
taue = (Lcm^2)/(mun*V);
Gph = (tau/taue)*(1+mup/mun); ip = e*eta*Phi*Gph;
fprintf('Phi=%.3e ph/s  tau_e=%.3e s  G=%.2f  ip=%.3e A (%.2f pA)\n', Phi, taue, Gph, ip, ip*1e12);
fprintf('  primary iph=%.3e A (%.3f pA)\n', e*eta*Phi, e*eta*Phi*1e12);

fprintf('\n===== PART B: Wheatstone =====\n');
aneg=-3e-4; apos=5e-4;
fprintf('S_bias=(apos-aneg)/2 = %.2e per K (%.2f mV/V/K)\n', (apos-aneg)/2, (apos-aneg)/2*1e3);
VB=5; Rr=1000; ENBW=1e3; dT=100;
Vo = VB*(apos-aneg)*dT/(2+(apos+aneg)*dT);
Rn=Rr*(1+aneg*dT); Rp=Rr*(1+apos*dT); Rout=2*Rn*Rp/(Rn+Rp);
un=sqrt(4*kB*(300+dT)*Rout*ENBW);
fprintf('T=400K: Vo=%.4f V Rout=%.1f un=%.2e V SNR=%.1f dB\n', Vo, Rout, un, 20*log10(Vo/un));

fprintf('\n===== PART D: MOS-Cap =====\n');
tox=30e-9;
o1=mosparams(3.9,3e16,tox,T);
o2=mosparams(7.5,0.5e16,tox,T);
fprintf('SiO2/3e16: Cox=%.3e CFB=%.3e Cmin=%.3e F/m2 (phiF=%.4f Wmax=%.1f nm)\n',o1(1),o1(2),o1(3),o1(4),o1(5)*1e9);
fprintf('  nF/cm2: Cox=%.1f CFB=%.1f Cmin=%.1f\n',o1(1)*1e5,o1(2)*1e5,o1(3)*1e5);
fprintf('SiN/0.5e16: Cox=%.3e CFB=%.3e Cmin=%.3e F/m2 (phiF=%.4f Wmax=%.1f nm)\n',o2(1),o2(2),o2(3),o2(4),o2(5)*1e9);
fprintf('  nF/cm2: Cox=%.1f CFB=%.1f Cmin=%.1f\n',o2(1)*1e5,o2(2)*1e5,o2(3)*1e5);
fprintf('dC swing case1=%.3e case2=%.3e  ratio=%.2f\n', o1(1)-o1(3), o2(1)-o2(3), (o2(1)-o2(3))/(o1(1)-o1(3)));
S59=2.303*kB*T/e;
fprintf('Nernst max=%.1f mV/pH; A(a=0.5)=%.1f B(a=0.9)=%.1f mV/pH ratio=%.1f\n', S59*1e3, 0.5*S59*1e3, 0.9*S59*1e3, 0.9/0.5);

function out = mosparams(epsox_rel, N_cm3, tox, T)
  e=1.602e-19; kB=1.381e-23; eps0=8.854e-12; ni=1.5e16; % m^-3
  epsSi=11.7*eps0; epsox=epsox_rel*eps0; N=N_cm3*1e6;
  Cox=epsox/tox;
  phiF=(kB*T/e)*log(N/ni);
  Wmax=sqrt(4*epsSi*phiF/(e*N));
  Cscmin=epsSi/Wmax; Cmin=Cox*Cscmin/(Cox+Cscmin);
  LD=sqrt(epsSi*(kB*T/e)/(e*N)); CscFB=epsSi/LD; CFB=Cox*CscFB/(Cox+CscFB);
  out=[Cox, CFB, Cmin, phiF, Wmax];
end
