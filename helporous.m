%计算多孔材料吸声系数
%T温度 phim 孔隙率 Alpham 曲折度 vclm黏性特征长度 tclm热特征长度 sigma0 20°下的静流阻
phim=0.9;Alpham=1.05;vclm=19.7e-6;tclm=36e-6;sigma0=489336;Rs=286.987;
gamma=1.4;H=50e-3;p0 = 101325;
i = 0;alpha =zeros(100,1);
T = [293;473;673;873];
%这个是我把comsol的有限元结果也导出来了，
%只要理论解的话把跟data有关的都注释掉就好了
data = xlsread("HelproComsolData.xlsx",'Sheet1');

%-----最后的计算--------
for j = 1:4
    rho0 = p0/(Rs*T(j));C0=sqrt(gamma*Rs*T(j));
    for f = 1:10:1000 
        [rhom,Km,sigmam] = matrix(T(j),phim,Alpham,vclm,tclm,sigma0,f,gamma);
        [rhop,Kp,phip] = hel(T(j),gamma,f);
        [rhodp,Kdp] = eff(rhom,Km,sigmam,phim,rhop,Kp,phip,f);
        Zdp = sqrt(rhodp*Kdp);
        kdp = 2*pi*f*sqrt(rhodp/Kdp);
        Zs = -1j*Zdp*cot(kdp*H);
        Z0 = rho0*C0;
        i=i+1;
        alpha(i,:) = 4*real(Zs)/Z0/((1+real(Zs)/Z0)^2+(imag(Zs)/Z0)^2);
    end
end
%------画图--------

color1 =[0, 3, 249]/255;
color2 =[1, 209, 137]/255;
color3 =[225, 175, 42]/255;
color4 =[249, 0, 4]/255;
p1=plot(1:10:1000,alpha(1:100,1),'Color',color1,'LineWidth',2.5);
hold on 
set(gcf,'Color','white');
set(gca,'linewidth',1.3);%全局效果放置于最前
set(gca,'FontSize',18);
set(gca,'XTick',0:100:1000);
set(gca,'XTicklabel',{0,[],200,[],400,[],600,[],800,[],1000},'FontSize',15);
set(gca,'YTick',0.0:0.1:1.0);
set(gca,'YTicklabel',{'0.0','',0.2,'',0.4,'',0.6,'',0.8,'','1.0'});
xlabel('Frequency (Hz)','fontsize',16,'fontname','Times','FontWeight','bold')
ylabel('Absorption coefficient','fontsize',16,'fontname','Times','FontWeight','bold')
set(gca,'XTickLabelRotation',0)
s1=scatter(data((1:2:100),1),data((1:2:100),2),60,color1,'MarkerFaceColor','white', ...
    'MarkerEdgeColor',color1,'LineWidth',2);
p2=plot(1:10:1000,alpha(101:200,1),'Color',color2,'LineWidth',2.5);
s2=scatter(data((1:2:100),1),data((1:2:100),4),60,color2,'MarkerFaceColor','white', ...
    'MarkerEdgeColor',color2,'LineWidth',2);
p3=plot(1:10:1000,alpha(201:300,1),'Color',color3,'LineWidth',2.5);
s3=scatter(data((1:2:100),1),data((1:2:100),6),60,color3,'MarkerFaceColor','white', ...
    'MarkerEdgeColor',color3,'LineWidth',2);
p4=plot(1:10:1000,alpha(301:400,1),'Color',color4,'LineWidth',2.5);
% plot(data(:,1),data(:,2),'Color',color1,'LineWidth',2.5,'Marker','o','MarkerFaceColor','white','MarkerEdgeColor',color1,'MarkerSize',10)
s4=scatter(data((1:2:100),1),data((1:2:100),8),60,color4,'MarkerFaceColor','white', ...
    'MarkerEdgeColor',color4,'LineWidth',2);

%设置图例我只会用最简单的方法设置图例 自己一个个挪动位置。。。。
legend1 = legend([p1,s1],'20°C(Theory)','20°C(FEM)','Location','SouthEast');
set(legend1, 'EdgeColor', 'none')
set(gca,'FontSize',12);
%第二个图例
a2 = axes('position',get(gca,'position'),'visible','off');
set(a2,'linewidth',1.3);
set(a2,'FontSize',12);
legend2 = legend(a2,[p2,s2],'200°C(Theory)','200°C(FEM)','Location','SouthEast');
set(legend2, 'EdgeColor', 'none')
%第三个图例
a4 = axes('position',get(gca,'position'),'visible','off');
set(a4,'linewidth',1.3);
set(a4,'FontSize',12);
legend3 = legend(a4,[p3,s3],'400°C(Theory)','400°C(FEM)','Location','SouthEast');
set(legend3, 'EdgeColor', 'none')
%第四个图例
a4 = axes('position',get(gca,'position'),'visible','off');
set(a4,'linewidth',1.3);
set(a4,'FontSize',12);
legend3 = legend(a4,[p4,s4],'600°C(Theory)','600°C(FEM)','Location','SouthEast');
set(legend3, 'EdgeColor', 'none')
%设置文本


% legendText=legend('20','20','200','200','400','400','600','600','Location','SouthEast');
% set(legendText,'Orientation','horizon')
% ah = axes('position',get(gca,'position'),'visible','off'); %新建一个坐标轴对象，位置和当前的坐标轴一致，并且设置为不可见，这样就不会覆原来绘制的图
function [rhom,Km,sigmam] = matrix(T,phim,Alpham,vclm,tclm,sigma0,f,gamma)
%计算基体的有效密度和有效体积模量
%T温度 phim 孔隙率 Alpham 曲折度 vclm黏性特征长度 tclm热特征长度 sigma0 20°下的静流阻
%f 频率 %eta动态粘度 gamma比热比 k 空气热导率 Cp 空气比热容 
p0 = 101325;
rho0 = p0/(286.987*T);
T0 = 293;
eata0 = 17.9e-6;
sigmam = sigma0*(T/T0)^0.6;
eta = eata0*(T/T0)^1.5*((120+T0)/(120+T));
An = [-0.0022758,1.1548e-4,-7.902528e-8,4.117025e-11,-7.4386433e-15];
Bn = [1047.63657,-0.372589265,9.45304214e-4,-6.0240944e-7,1.2858961e-10];
Tn = [1,T,T^2,T^3,T^4];
k = sum(An.*Tn);
Cp = sum(Bn.*Tn);
w=2*pi*f;
rhom = Alpham*rho0/phim*(1+phim*sigmam/(1i*w*rho0*Alpham) ...
        *sqrt(1+1i*4*w*eta*rho0*Alpham^2/(sigmam^2*vclm^2*phim^2)));

Km = gamma*p0/phim*(gamma-(gamma-1)*(1-1i*8*k*sqrt(1+1i*w*rho0*Cp*tclm^2/(16*k)) ...
                                                       /(w*rho0*Cp*tclm^2))^-1)^-1;
end
function [rhop,Kp,phip] = hel(T,gamma,f)
%算螺旋结构的有效密度和有效体积模量
%phip 螺旋孔隙率 Alphap 螺旋曲折度 beta 和 pr 是普朗特数
%定义变量
H=50e-3;D=35e-3;l=15e-3;p=15e-3;d=10e-3;
cphi = p/sqrt((pi*l)^2+p^2);
p0 = 101325;
rho0 = p0/(286.987*T);
T0 = 293;
eata0 = 17.9e-6;
eta = eata0*(T/T0)^1.5*((120+T0)/(120+T));
w=2*pi*f;
An = [-0.0022758,1.1548e-4,-7.902528e-8,4.117025e-11,-7.4386433e-15];
Bn = [1047.63657,-0.372589265,9.45304214e-4,-6.0240944e-7,1.2858961e-10];
Tn = [1,T,T^2,T^3,T^4];
k = sum(An.*Tn);
Cp = sum(Bn.*Tn);
Alphap = 1/cphi^2;
beta = d*sqrt(w*rho0/eta)/2;
pr = eta*Cp/k;
phip = d^2/(D^2*cphi);
rhop = Alphap*rho0/phip/(1-2*besselj(1,beta*sqrt(-1j))/(beta*sqrt(-1j)*besselj(0,beta*sqrt(-1j))));
Kp = gamma*p0/phip/(1+(gamma-1)*2*besselj(1,beta*sqrt(-1j*pr))/(beta*sqrt(-1j*pr)*besselj(0,beta*sqrt(-1j*pr))));
end

function [rhodp,Kdp] = eff(rhom,Km,sigmam,phim,rhop,Kp,phip,f)
%计算总有效密度和有效体积模量
fdw = f2(Km,sigmam,phim,phip,f);
    function fdw = f2(Km,sigmam,phim,phip,f)
        w=2*pi*f;
        p0 = 101325;
        H=50e-3;D=35e-3;l=15e-3;p=15e-3;d=10e-3;
        cphi = p/sqrt((pi*l)^2+p^2);
        D0 = (1-phip)*pi*D^2/4;
        wd = (1-phip)*p0/(phim*sigmam*D0);
        Omegam = pi*D^2*H*(1-phip)/4;
        Gammamp = (pi*D^2-pi*d^2/cphi)/4+pi*d*H/cphi;
        lambdad = 2*Omegam/Gammamp;
        Md = 8*D0/(lambdad^2*(1-phip));
        fdw = 1-1j*w*p0/(wd*phim*Km)/(1j*w*p0/(wd*phim*Km)+sqrt(1+1j*w*p0/(wd*phim*Km)*Md/2));    
    end
rhodp = (rhop^-1+(1-phip)*rhom^-1)^-1;
Kdp = (Kp^-1+fdw*(1-phip)*Km^-1)^-1;

end









