clear all
close all

%% Quasistatic comparison - September '23
% Henry Putley, PhD, Mathematics, Imperial College London -
% h.putley18@imperial.ac.uk

%% Setup variables
d=1; a1=[d 0]; a2=[0 d]; delta=pi/2; h=200; hk=300; % 200, 300 (respectively)
rc=@(f) sqrt(f*d^2/pi);
kBs=[zeros(1,h),linspace(0.001,pi/d,h);linspace(pi/d,0.001,h),zeros(1,h)]; % Path N -> Gamma -> X (principal axes)
k_perpspan=linspace(0.01,pi/d-0.001,hk);
k_perps=k_perpspan;

%% Dipole Formulae
eps_x = @(nr,f) (nr^2*(1+f)^2-f^2+1)/(nr^2+1-f^2*(nr^2-1)+2*f); 
eps_y =@(f) (1+f)/(1-f);
alpha_x = @(nr,f) sqrt(1/eps_x(nr,f));
alpha_y = @(f) sqrt(1/eps_y(f));

fs=[0.1,0.3,0.5];
nrs=[2,1/2,1];

%% Multipole - uncomment to use
% M=6; tmpt=2*M+2; mpo=M+1; tmpo=2*M+1;
% %d=1; f=fs(2); r_c=sqrt(f*d^2/pi); 
% %r_c=0.3;
% %f=pi*r_c^2/(d^2);
% f=fs(2); r_c=sqrt(f*d^2/pi);
% n_r=nrs(2);
% res=pi/(2*n_r*r_c)
% heigs=1;
% 
% % pre-compute Bessel functions
% BesselJs=zeros(hk,2*tmpt); BesselHs=zeros(hk,2*tmpt); dBesselJs=zeros(hk,2*tmpt); dBesselHs=zeros(hk,2*tmpt);
% BesselJsin=zeros(hk,2*tmpt); dBesselJsin=zeros(hk,2*tmpt);
% for j=1:length(k_perps)
%     k_perp=k_perps(j);
%     BesselJs(j,:)=besselj(-tmpt:tmpt-1,k_perp*r_c);
%     BesselHs(j,:)=besselh(-tmpt:tmpt-1,k_perp*r_c);
%     dBesselJs(j,:)=(1/2)*(besselj((-tmpt:tmpt-1)-1,k_perp*r_c)-besselj((-tmpt:tmpt-1)+1,k_perp*r_c));
%     dBesselHs(j,:)=(1/2)*(besselh((-tmpt:tmpt-1)-1,k_perp*r_c)-besselh((-tmpt:tmpt-1)+1,k_perp*r_c));
%     BesselJsin(j,:)=besselj(-tmpt:tmpt-1,n_r*k_perp*r_c);
%     dBesselJsin(j,:)=(1/2)*(besselj((-tmpt:tmpt-1)-1,n_r*k_perp*r_c)-besselj((-tmpt:tmpt-1)+1,n_r*k_perp*r_c));
% end
% for kk=1:1
% % BesselYs=imag(BesselHs); dBesselYs=imag(dBesselHs);
% xi = sqrt(sum(a1.^2)); m = 0:tmpo;
% acc = [1 5.*ones(size(m(2:end)))];
% % set up reciprocal vectors:
% A = a1(1)*a2(2) - a1(2)*a2(1);
% b1 = 2*pi/A.*[a2(2) -a2(1)];
% b2 = 2*pi/A.*[-a1(2) a1(1)];
% HMAX=20; [h1, h2] = meshgrid([-HMAX:HMAX]);
% Kh1 = reshape(h1.*b1(1) + h2.*b2(1),length(h1).^2,1);
% Kh2 = reshape(h1.*b1(2) + h2.*b2(2),length(h2).^2,1);
% 
% Khvec = [Kh1 Kh2];
% accx = repmat(acc,length(Kh1),1);
% mx = repmat(m,length(Kh1),1);
% macx = mx + accx;
% end % hide code
% 
% for k=1:length(kBs)
%     kB=[kBs(1,k),kBs(2,k)];
%     Qhvec(:,1) = Khvec(:,1) + kB(1); Qhvec(:,2) = Khvec(:,2) + kB(2);
%     [Th, Qh] = cart2pol(Qhvec(:,1),Qhvec(:,2));
%     Qhx = repmat(Qh,1,length(m)); Thx = repmat(Th,1,length(m));
%     BesselAccs=besselj(macx,Qhx.*xi);
%     clear eigsP
%     
% %     if k<=h
% %         k_perps=linspace(alpha_x(n_r,f)*abs(kB(2))-0.05,alpha_x(n_r,f)*abs(kB(2))+0.05,20);
% %     else
% %         k_perps=linspace(alpha_y(f)*abs(kB(1))-0.05,alpha_y(f)*abs(kB(1))+0.05,20);
% %     end
%     for j=1:length(k_perps)
%         k_perp=k_perps(j);
% %         if k>h && k_perp>2.09706
% %             break
% %         end
% %         elseif k>=1 && k_perp>2.35004
% %             break
% %         end
%             
%         [S_mp,Sy_mp]=LatticeSumMat_faster(M,k_perp,xi,m,mx,A,Qhx,Thx,acc,accx,BesselAccs);
%         if norm(kB)==k_perp
%             continue
%         end
%         P=TwoDMetaMatS_nr(M,delta,BesselJs(j,:),BesselHs(j,:),dBesselJs(j,:),dBesselHs(j,:),S_mp,n_r,BesselJsin(j,:),dBesselJsin(j,:));
%         [~,eigsP(j,:)]=eigs(P,1,'smallestabs');
%         eigsDim(k,j,:)=abs(eigsP(j,:));
%     end
%     [mineigs(k,1:heigs),I(k,1:heigs)]=mink(eigsP,heigs,'ComparisonMethod','abs');
%     ks(k,1:heigs)=k_perps(I(k,1:heigs));
% end
% 
% 
% figure
% hold on
% plot([1:h,h+1:2*h],ks(:,1),'-k')
% xlim([0 2*h])
% xticks([1 h 2*h])
% xticklabels(["N","\Gamma","X"])
% ylabel("k")
% 
% figure
% pcolor([1:h,h+1:2*h],k_perpspan,eigsDim.')
% load('davos.mat');
% colormap(davos)
% colorbar
% shading interp
% hold on
% plot(1:h,alpha_x(n_r,f)*abs(kBs(2,1:h)),'w--','LineWidth',2)
% plot(h+1:2*h,alpha_y(f)*abs(kBs(1,h+1:2*h)),'w--','LineWidth',2)

% save('Numerics_b_nr2.mat','ks','k_perpspan','eigsDim','h','f','n_r')


%% Comparison 
set(groot, 'DefaultFigurePosition', [400, 400, 600, 450])
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

h=200;
kBs=[zeros(1,h),linspace(0.001,pi/d,h);linspace(pi/d,0.001,h),zeros(1,h)]; % Path N -> Gamma -> X (principal axes)
txtsf=["Numerics_a_fp1.mat","Numerics_a_fp3.mat","Numerics_a_fp5.mat"];
txtsnr=["Numerics_b_nr2.mat","Numerics_b_nrp5.mat","Numerics_b_nrp05.mat"];
colors=['r','b','k'];
f_files=["fig6a_f01.fig","fig6a_f03.fig","fig6a_f05.fig"];
nr_files=["fig6a_f03_nr2.fig","fig6a_f03_nr05.fig"];

Fem_X=linspace(1,400,140);

% fhg=figure;
% hold on
for idx=1:3
    open(f_files(idx))
    D=get(gca,'Children'); %get the handle of the line object
    XData=get(D,'XData'); %get the x data
    YData=get(D,'YData'); %get the y data
    %Data=[XData' YData'];
    Y=cell2mat(YData(end));
    Y2=cell2mat(YData(length(YData)/2));
    Fem_Y=[Y,Y2];
    
    figure(2)
    hold on
    f=fs(idx); 
    data=load(txtsf(idx));
    plot([1:h,h+1:2*h],[alpha_x(1,f)*abs(kBs(2,1:h)),alpha_y(f)*abs(kBs(1,h+1:2*h))],'-','LineWidth',1,'color',colors(idx))
    scatter([2:6:2*h],data.ks(2:6:end,1),sprintf('+ %s',colors(idx)))
    scatter(Fem_X(2:3:end),Fem_Y(2:3:end),25,sprintf('sq %s',colors(idx)),"filled")
end
ylim([0 pi])
xlim([0 2*h])
xticks([1 h 2*h])
xticklabels(["$\mathrm{N}$","$\mathrm{\Gamma}$","$\mathrm{X}$"])
ylabel('$k$')
% text(405,2.85,'$f=0.1$')
% text(405,2.3,'$f=0.3$')
% text(405,1.8,'$f=0.5$')
% x1=310;
% x2=380;
% y1=0.1;
% y2=0.1+0.0105*(70/400)*300+0.1;
% x = [x1, x2, x2, x1, x1];
% y = [y1, y1, y2, y2, y1];
% plot(x, y, 'k-', 'LineWidth', 1);
% text(341,(0.1+0.0105*(70/400)*300+0.1)/2+0.03,'$\mathrm{\Gamma}$')
% text(341,(0.1+0.0105*(70/400)*300+0.1)+0.08,'$\mathrm{N}$')
% text(384,(0.1+0.0105*(70/400)*300+0.1)/2+0.04,'$\mathrm{X}$')
% annotation('rectangle',[0.725 .13 .18*450/1200 .18],'LineWidth',1,'FaceColor',[0.9,0.9,0.9]);
% annotation('rectangle',[0.725 .13 .18*450/600 .18/2],'LineWidth',1,'FaceColor',[0.9,0.9,0.9]);

annotation('rectangle',[0.725 .13 .18*450/600 .18],'LineWidth',1,'FaceColor',[0.9,0.9,0.9]);
dim = [.725+.18*450/1200-0.019 .157 .8*450/600 .09];
str = '$\mathrm{\Gamma}$';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','Interpreter','latex');
dim = [.725+.18*450/1200-0.019 .27 .18*450/600 .09];
str = '$\mathrm{N}$';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','Interpreter','latex');
dim = [.725+.18*450/1200+0.1*450/600-0.009 .157 .18*450/600 .09];
str = '$\mathrm{X}$';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','Interpreter','latex');
annotation('line',[0.725+.18*450/1200,0.725+.18*450/1200],[0.31,0.3])
annotation('line',[0.725+.18*450/600,0.725+.18*450/600-0.01*450/600],[0.13+0.09,0.13+0.09])
box on
hg = zeros(3, 1);
hg(1) = plot(NaN,NaN,'-k');
hg(2) = plot(NaN,NaN,'+k');
hg(3) = plot(NaN,NaN,'sqk','MarkerFaceColor','k');
hg(4) = plot(NaN,NaN,'.r','MarkerSize',20);
hg(5) = plot(NaN,NaN,'.b','MarkerSize',20);
hg(6) = plot(NaN,NaN,'.k','MarkerSize',20);
legend(hg, 'Dipole','Multipole','FEM','$f=0.1$','$f=0.3$','$f=0.5$','Location','SW');
legend boxoff
set(findall(gcf,'-property','FontSize'),'FontSize',12)
% axis square
% exportgraphics(gcf,'C:\Documents\PhD\Thesis\Thesis_Figures\compare_f.eps','BackgroundColor','none','ContentType','vector')
% exportgraphics(gcf,'C:\Documents\PhD\Thesis\Thesis_Figures\compare_f.png','BackgroundColor','white')



figure
hold on
for idx=1:2
    open(nr_files(idx))
    D=get(gca,'Children'); %get the handle of the line object
    XData=get(D,'XData'); %get the x data
    YData=get(D,'YData'); %get the y data
    %Data=[XData' YData'];
    Y=cell2mat(YData(end));
    Y2=cell2mat(YData(length(YData)/2));
    Fem_Y=[Y,Y2];
    
    figure(5)
    hold on
    nr=nrs(idx);
    data=load(txtsnr(idx));
    data.ks(data.ks>=3.14059)=NaN;
    plot([1:h,h+1:2*h],[alpha_x(nr,0.3)*abs(kBs(2,1:h)),alpha_y(0.3)*abs(kBs(1,h+1:2*h))],'-k','LineWidth',1,'color',colors(idx))
    scatter([2:6:2*h],data.ks(2:6:end,1),sprintf('+ %s',colors(idx)))
    scatter(Fem_X(2:3:end),Fem_Y(2:3:end),25,sprintf('sq %s',colors(idx)),"filled")
end
ylim([0 pi])
xlim([0 2*h])
xticks([1 h 2*h])
xticklabels(["$\mathrm{N}$","$\mathrm{\Gamma}$","$\mathrm{X}$"])
ylabel('$k$')
% text(-45,2.72,'$n_r=2$')
% text(1,3.25,'$n_r=\frac{1}{2}$')
% text(50,3.25,'$n_r=\frac{1}{20}$')
annotation('rectangle',[0.725 .13 .18*450/600 .18],'LineWidth',1,'FaceColor',[0.9,0.9,0.9]);
dim = [.725+.18*450/1200-0.019 .157 .8*450/600 .09];
str = '$\mathrm{\Gamma}$';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','Interpreter','latex');
dim = [.725+.18*450/1200-0.019 .27 .18*450/600 .09];
str = '$\mathrm{N}$';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','Interpreter','latex');
dim = [.725+.18*450/1200+0.1*450/600-0.009 .157 .18*450/600 .09];
str = '$\mathrm{X}$';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','Interpreter','latex');
annotation('line',[0.725+.18*450/1200,0.725+.18*450/1200],[0.31,0.3])
annotation('line',[0.725+.18*450/600,0.725+.18*450/600-0.01*450/600],[0.13+0.09,0.13+0.09])
box on
hg = zeros(3, 1);
hg(1) = plot(NaN,NaN,'-k');
hg(2) = plot(NaN,NaN,'+k');
hg(3) = plot(NaN,NaN,'sqk','MarkerFaceColor','k');
hg(4) = plot(NaN,NaN,'.r','MarkerSize',20);
hg(5) = plot(NaN,NaN,'.b','MarkerSize',20);
hg(6) = plot(NaN,NaN,'.k','MarkerSize',20);
legend(hg, 'Dipole','Multipole','FEM','$n_r=2$','$n_r=\frac{1}{2}$','Location','SW');
legend boxoff
set(findall(gcf,'-property','FontSize'),'FontSize',12)
% exportgraphics(gcf,'C:\Documents\PhD\Thesis\Thesis_Figures\compare_nr.eps','BackgroundColor','none','ContentType','vector')
% exportgraphics(gcf,'C:\Documents\PhD\Thesis\Thesis_Figures\compare_nr.png','BackgroundColor','white')


