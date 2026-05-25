% Generate all figures for ET4260 Unit 4 homework report
clear; clc; close all
set(groot,'defaultAxesFontSize',11,'defaultLineLineWidth',1.6, ...
          'defaultAxesTickLabelInterpreter','latex', ...
          'defaultLegendInterpreter','latex','defaultTextInterpreter','latex');
outdir = fullfile('..','report','figures');

%% ---------- Fig 1: Task 3b  SNR vs Temperature ----------
e=1.602e-19; kB=1.381e-23;
VB=5; R=1000; ENBW=1e3; aneg=-3e-4; apos=5e-4; T0=300;
T = linspace(300,400,1001); dT = T - T0;
Vo  = VB*(apos-aneg).*dT ./ (2 + (apos+aneg).*dT);          % exact Vo(T)
Rn  = R*(1+aneg.*dT);  Rp = R*(1+apos.*dT);
Rout= 2*Rn.*Rp./(Rn+Rp);                                    % bridge output resistance
un  = sqrt(4*kB.*T.*Rout*ENBW);
SNR = 20*log10(Vo./un);

figure('Units','centimeters','Position',[2 2 13 8.5]);
plot(T,SNR,'Color',[0 0.35 0.75]); grid on; box on
xlabel('Temperature $T$ [K]'); ylabel('SNR [dB]');
title('Wheatstone-bridge SNR vs.\ temperature (exact)');
xlim([300 400]); ylim([0 130]);
exportgraphics(gcf,fullfile(outdir,'q3b_snr.pdf'),'ContentType','vector');

% NOTE: the Part D C-V figures (Task 4a labels and Task 4b sketch) are NOT plotted
% here. They are vector overlays drawn directly on the original assignment figures
% (report/figures/cv_fig3_blank.png, cv_fig4_blank.png) via TikZ in
% report/figures/q4a_fig.tex and q4b_fig.tex.

%% ---------- Fig: Part E Q3  grid spacing on curve (c) ----------
xa = linspace(-4,-2,40);  ya = (xa+4)/2 - 1;                 % line (-4,-1)->(-2,0)
xb = linspace(-2,2,200);  yb = -sqrt(max(4-xb.^2,0));        % lower semicircle
xc = linspace(2,4,40);    yc = (xc-2);                       % up to (4,2)
xd = linspace(4,6,40);    yd = 2 - (xd-4);                   % down to (6,0)
figure('Units','centimeters','Position',[2 2 14 8]);
plot(xa,ya,'k',xb,yb,'k',xc,yc,'k',xd,yd,'k','LineWidth',1.6); hold on; grid on; box on
% proposed nodes: coarse on straight parts, fine on the semicircle, nodes at kinks
ns_line = [-4 -2];                       % straight: endpoints only
ns_semi = linspace(-2,2,11);             % curved: dense
ns_up   = [2 4]; ns_dn = [4 6];          % straight: endpoints only
plot(ns_line,(ns_line+4)/2-1,'ro','MarkerFaceColor','r','MarkerSize',6);
plot(ns_semi,-sqrt(max(4-ns_semi.^2,0)),'ro','MarkerFaceColor','r','MarkerSize',6);
plot(ns_up,(ns_up-2),'ro','MarkerFaceColor','r','MarkerSize',6);
plot(ns_dn,2-(ns_dn-4),'ro','MarkerFaceColor','r','MarkerSize',6);
xlabel('$x$'); ylabel('$y$'); title('Proposed node placement (red) for linear hat functions');
axis equal; xlim([-4.5 6.5]); ylim([-2.5 2.5]);
text(-3.9,-1.4,'coarse','FontSize',9); text(-1.6,-2.25,'fine (high curvature)','FontSize',9);
text(4.1,1.3,'coarse','FontSize',9);
exportgraphics(gcf,fullfile(outdir,'fem_grid.pdf'),'ContentType','vector');

disp('All figures written.');
