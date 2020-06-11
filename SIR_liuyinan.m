global date cases %定义真实数据为全局变量供函数调用
date=[0:70]';
cases=[7;12;15;15;16;18;23;24;24;25;27;28;28;28;28;28;29;30;31;46;80;155;345;601;762;892;1146;1595;2022;2931;3526;4212;4812;5328;5766;6284;6767;7134;7382;7513;7755;7869;7979;8086;8162;8236;8320;8413;8565;8652;8799;8897;8961;9037;9137;9241;9332;9478;9583;9661;9786;9786;9976;10062;10156;10237;10284;10331;10384;10423;10450];
bmin=[0.0000005,0.005,8000,0.001,0]';
bmax=[0.000050,0.5,15000,50,0]';
ff=optimset; 
ff.LargeScale='off';
ff.TolFun=1e-30; 
ff.Tolx=1e-20; 
ff.TolCon=1e-30; 

b_0=[0.0000046,0.02,10000,10,0]';%需要拟合的参数，依次为 lambda/N miu S0 I0 
%有约束极小优化参数
[btarget,y]=fmincon(@parameter,b_0,[],[],[],[],bmin,bmax,[],ff);
%优化后参数计算数据
date=0:120;
dateRaw=0:length(cases)-1;
[date,X]=ode45(@SIR,date,btarget(3:5),[],btarget(1:2));
S=X(:,1);
I=X(:,2);
R=X(:,3);
Rsquare=corrcoef(R(1:length(cases)),cases);

figure
scatter(dateRaw,cases,'SizeData',50)
hold on 
plot(date,S,'LineWidth',2)
plot(date,I,'LineWidth',2)
plot(date,R,'LineWidth',2)

text_shown={};%元胞数组传入多行图注
text_shown{1}=['\beta = ',num2str(btarget(1))];
text_shown{2}=['\mu = ',num2str(btarget(2))];
text_shown{3}=['S_0 = ',num2str(btarget(3))];
text_shown{4}=['I_0 = ',num2str(btarget(4))];
text_shown{5}=['r^2 = ',num2str(Rsquare(1,2))];

h=legend('# of cases','S_{pre}','I_{pre}','R_{pre}');

g=text(0.5,0.5,text_shown,'FontSize',14);
g.Units='normalized';
g.Position=[0.8,0.5];

xlabel('Date from 20190131 (day)');
ylabel('Cases in Korea');
set(gca,'fontsize',16);
set(h,'Box','off');
pbaspect([1 0.618 1])

clear global

function errx=parameter(b)
%为计算误差平方和函数
global date cases
datax=outbreak(b,date);
errx=sum((cases-datax).^2);
end
function X=outbreak(b,date)
%为单独计算隔离者人数函数
[~,X]=ode45(@SIR,date,b(3:5),[],b(1:2));
X=X(:,3);
end
function dX=SIR(t,X,c)
%为SIR模型函数
dX=[-c(1)*X(1)*X(2);
c(1)*X(1)*X(2)-c(2)*X(2);
c(2)*X(2)];
end