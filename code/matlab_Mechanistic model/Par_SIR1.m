%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%下载工作区ShanghaiMCMC_history.mat，然后运行此程序
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
Data_case=csvread('newly_infected.csv');
Dat_hospital=Data_case';
Dat_hospital_ini=Dat_hospital;


Num1=1;Num2=length(Dat_hospital);
data1=zeros(4,length(Dat_hospital(Num1:Num2)));
data1(1,:)=0:length(Dat_hospital(Num1:Num2))-1;
data1(2,:)=Dat_hospital(Num1:Num2);


Sr0_total=24750000;%总人口
R0=0;%初始恢复者
par2=[Sr0_total R0];
%%
%求解最优参数
%        %   beta    γ     k      I0    Se0
% par1guess=[0.5         0.3    10    3000];   
%        lb=[0.1         0.01   1     2000]; 
%        ub=[1           1      100   400000];
% 
params = {
%     name,  init,       min,     max,      mu,          sig, target?, local?
    {'beta',  0.5,       0.01,     1,       0.5,        Inf}
    {'k',     0.3,        0.01,    1,       0.3,        Inf}
    {'I0',    2,          0.01,    10,      2,        Inf}
    {'Se0',   200000,           100000,  410000,     390000,        Inf}
 %   {'beta',  0.85,       0.4,     1,       0.7,        Inf}
  %  {'k',     0.3,       0.01,    1,       0.3,        Inf}
   % {'I0',    5,         0.01,     12,       1,        Inf}
    %{'Se0',   200000,    100000,  410000,  300000,        Inf}

    };

ydata=[(1:129)',data1(2,:)'];ylabels={'time','newly_infected'};
zdata=par2';zlabels={'par2'};
data=struct('ydata',ydata,'zdata',zdata);

model.ssfun = @kuzdist;%似然函数


options.nsimu = 300000;%生成的参数链个数
[results, chain, s2chain]= mcmcrun(model,data,params,options);
options.nsimu = 50000;
[results, chain, s2chain] = mcmcrun(model,data,params,options, results);
options.nsimu = 50000;
[results, chain, s2chain] = mcmcrun(model,data,params,options, results);

options.nsimu = 50000;
[results, chain, s2chain] = mcmcrun(model,data,params,options, results);

%参数后验分布
figure
mcmcplot(chain,[],results,'pairs');
figure
mcmcplot(chain,[],results,'denspanel',2);



%%
%ps:第一列为参数估计结果,pCI:参数置信区间
ps=chainstats(chain,results)
%置信区间：均值+-95%分位数*标准差。小样本（<30）:用t分布求置信区间；大样本，用正态分布95%分位数求置信区间即1.96
pCI=[ps(1,1)-1.96*ps(1,2) ps(1,1)+1.96*ps(1,2);ps(2,1)-1.96*ps(2,2) ps(2,1)+1.96*ps(2,2);
    ps(3,1)-1.96*ps(3,2) ps(3,1)+1.96*ps(3,2);ps(4,1)-1.96*ps(4,2) ps(4,1)+1.96*ps(4,2)];

%参数估计结果及其置信区间
bestpar1=[ps(:,1),pCI ]



%%
%Sr,Se,I,R,NI的置信区间
dataSe=zeros(129,49300);dataI=zeros(129,49300);dataR=zeros(129,49300);dataNI=zeros(128,49300);%dataSr=zeros(129,49300);
for i=700:50000
par1=chain(i,:);
X0=[par1(4),par1(3),par2(2),5,par2(1)-par1(4)];
[t,x]=ode45(@SIRmodel,1:129,X0,[],par1,par2);

dataSe(:,i)=x(:,1);dataI(:,i)=x(:,2);dataR(:,i)=x(:,3);dataNI(:,i)=diff(x(:,4));
%dataSr(:,i)=x(:,5);
end
dataSr=zeros(129,49999);
for i=2:50000
par1=chain(i,:);
X0=[par1(4),par1(3),par2(2),5,par2(1)-par1(4)];
[t,x]=ode45(@SIRmodel,1:129,X0,[],par1,par2);

dataSr(:,i)=x(:,5);
end

mSe=mean(dataSe,2);sSe=std(dataSe,0,2);UconfSe=mSe+1.96.*sSe;LconfSe=mSe-1.96.*sSe;
mI=mean(dataI,2);sI=std(dataI,0,2);UconfI=mI+1.96.*sI;LconfI=mI-1.96.*sI;
mR=mean(dataR,2);sR=std(dataR,0,2);UconfR=mR+1.96.*sR;LconfR=mR-1.96.*sR;
mNI=mean(dataNI,2);sNI=std(dataNI,0,2);UconfNI=mNI+1.96.*sNI;LconfNI=mNI-1.96.*sNI;
mSr=mean(dataSr,2);sSr=std(dataSr,0,2);UconfSr=mSr+1.96.*sSr;LconfSr=mSr-1.96.*sSr;

%%
%图4
par1=bestpar1(:,1);
X0=[par1(4),par1(3),par2(2),5,par2(1)-par1(4)];
[t,x]=ode45(@SIRmodel,1:0.1:129,X0,[],par1,par2);

fillcolor=[0.58 0.80 1];
figure
subplot(4,2,[1 2])
aaa=1:128;
FNI=[5.0, 1.0916028, 1.4361897, 1.8895411, 2.486003, 3.2707539, 4.303029, 5.661455, 7.447878, 9.799267, 12.890148, 16.959621, 22.308517, 29.344719, 38.605568, 50.763824, 66.77243, 87.81406, 115.4144, 151.74075, 199.43817, 261.87213, 343.86328, 451.34863, 591.47485, 774.4226, 1013.24854, 1322.0635, 1720.6978, 2234.1206, 2888.4365, 3711.1445, 4742.8086, 6010.617, 7534.5723, 9324.801, 11385.852, 13651.969, 16041.859, 18462.312, 20744.68, 22766.742, 24431.453, 25646.594, 26385.734, 26660.672, 26529.61, 26059.812, 25318.562, 24371.219, 23280.656, 22100.375, 20858.344, 19581.094, 18295.469, 17019.562, 15766.5625, 14553.656, 13384.406, 12262.0625, 11190.0, 10171.375, 9209.875, 8307.5, 7469.1875, 6695.5625, 5986.5625, 5341.8125, 4761.25, 4244.8125, 3795.4375, 3402.3125, 3056.375, 2752.25, 2484.5, 2247.9375, 2037.6875, 1853.125, 1685.9375, 1533.875, 1395.9375, 1271.25, 1159.0625, 1058.4375, 968.5625, 888.4375, 817.3125, 751.625, 690.4375, 634.4375, 583.3125, 536.6875, 494.1875, 455.375, 419.875, 387.5, 357.5, 329.9375, 304.5, 281.125, 259.6875, 240.0, 222.0, 205.5625, 190.3125, 176.25, 163.125, 151.0625, 139.875, 129.625, 120.0625, 111.3125, 103.1875, 95.8125, 88.875, 82.5, 76.5, 71.0, 65.8125, 61.125, 56.6875, 52.625, 48.875, 45.4375, 42.1875, 39.1875, 36.4375, 33.8125, 31.4375];
fill([aaa fliplr(aaa)],[UconfNI' fliplr(LconfNI')],fillcolor,'edgealpha',0,'facealpha',0.7);%[]中必须都为行向量
hold on
newly_infected=x(11:end,4)-x(1:1271,4);%diff(x(:,4));
plot(1:0.1:128,newly_infected,'b-','LineWidth',2)
plot(aaa,FNI(2:end),'r--','LineWidth',2)
plot(aaa,data1(2,2:end)','go','LineWidth',2,'MarkerSize',2)
legend('95% CI','Mechanistic model','UDE model','data');ylabel('Daily new infections');
set(gca,'XTick',[ 5 36 66 97 127]);xlim([1 129]);xlabel('Date(days)')
%set(gca,'XTickLabel',{'','2018.01.31','2020.02.31'});
set(gca,'XTickLabel',{'2022/3/1','4/1','5/1','6/1','7/1','FontSize',10});%,'FontSize',14
title('(a)','Position',[-5 40000])
set(gca,'FontName','Times New Roman','FontSize',14)
grid on
grid minor

subplot(423)
 FRt=[2.37174854419768	2.37174002901214	2.37172882618447	2.37171408737793	2.37169469647080	2.37166918624733	2.37163562348400	2.37159147159750	2.37153338221923	2.37145697315584	2.37135644602214	2.37122422121295	2.37105030507187	2.37082152446827	2.37052072896991	2.37012513837943	2.36960499472261	2.36892154855478	2.36802330108429	2.36684323823284	2.36529467440104	2.36326283063871	2.36059855340977	2.35711169704274	2.35255412703324	2.34660427184999	2.33886342416227	2.32882612390418	2.31585676039100	2.29919331916921	2.27795453715801	2.25108804071718	2.21748061491659	2.17603994301975	2.12579804735478	2.06599700124592	1.99650124231186	1.91786517094944	1.83133002948071	1.73906169872574	1.64368812419007	1.54796916299069	1.45460394169009	1.36581041101598	1.28315218620012	1.20750319487770	1.13914608321281	1.07789482796138	1.02322749480205	0.974489029222811	0.930935926069367	0.891862212379694	0.856615504316282	0.824599574132889	0.795294157311928	0.768262418294289	0.743115114053051	0.719530614313893	0.697247983518637	0.676061088230621	0.655814958920868	0.636401097194381	0.617725033277619	0.599706511870423	0.582407836788167	0.566054316485794	0.551038840465108	0.537925704939699	0.527432845767839	0.519519002921219	0.513673275100401	0.509597650950341	0.507023057875969	0.505710906468929	0.505454541762261	0.506081651666821	0.507458626782231	0.509446454543090	0.511918541077694	0.514766169594401	0.517898008742586	0.521239991576958	0.524734892303703	0.528341750142151	0.532035446611113	0.535806045539678	0.539614808856497	0.543417971631073	0.547192514920706	0.550919617961551	0.554584709668348	0.558177321726425	0.561690918066445	0.565122003969999	0.568463863810178	0.571708984063224	0.574851647026393	0.577887995247692	0.580815905765532	0.583634941320281	0.586346440914108	0.588953334491316	0.591460317825014	0.593868847573652	0.596177211316563	0.598386371202053	0.600497689853195	0.602513025278027	0.604434788702095	0.606265812615423	0.608009353904378	0.609669300311923	0.611249608712119	0.612751978528402	0.614178721252724	0.615532136980613	0.6168147  0.6180    0.6192    0.6203    0.6213    0.6223    0.6232    0.6240    0.6249    0.6256    0.6264    0.6270    0.6277];
Rt=5*par1(1).*x(:,1)./(x(:,1)+x(:,2)+x(:,3));
plot(1:0.1:129,Rt,'b-','LineWidth',2)
hold on
plot(1:129,FRt,'r--','LineWidth',2)
%legend({'Mechanism model'});
ylabel('R_t');
line([1 129],[1 1],'linestyle','--','color',[0,0,0])
set(gca,'XTick',[ 6 37 67 98 128]);xlim([1 129]);xlabel('Date(days)')
%set(gca,'XTickLabel',{'','2018.01.31','2020.02.31'});
set(gca,'XTickLabel',{'2022/3/1','4/1','5/1','6/1','7/1'});
title('(b)','Position',[-13 2.5]);legend({'Mechanistic model','UDE model','R_t=1'})
set(gca,'FontName','Times New Roman','FontSize',14)
grid on
grid minor

subplot(424)
rate=par1(2).*par1(1).*x(:,2).*x(:,5)./par2(1);
Frate=[0	-6	-8	-6	-8	-8	-8	-10	-8	-12	-12	-16	-18	-24	-28	-38	-48	-62	-80	-104	-136	-178	-230	-302	-392	-510	-662	-858	-1110	-1432	-1834	-2346	-2982	-3756	-4688	-5798	-7072	-8498	-10044	-11652	-13254	-14770	-16134	-17270	-18114	-18662	-18912	-18854	-18502	-17908	-17114	-16148	-15052	-13856	-12602	-11314	-10024	-8744	-7490	-6284	-5140	-4072	-3086	-2184	-1426	-872	-586	-628	-1056	-1548	-1900	-2188	-2424	-2606	-2750	-2856	-2942	-2994	-3016	-3010	-2984	-2942	-2890	-2834	-2782	-2734	-2680	-2606	-2534	-2456	-2376	-2296	-2218	-2142	-2066	-1988	-1912	-1834	-1758	-1684	-1610	-1542	-1478	-1412	-1350	-1288	-1228	-1168	-1110	-1058	-1004	-954	-906	-860	-816	-772	-732	-692	-654	-616	-584	-550	-520	-490	-462	-436	-408	-386	-362];
plot(1:0.1:129,rate,'b-','LineWidth',2)
hold on
plot(1:129,-Frate,'r--','LineWidth',2)
legend({'Mechanistic model','UDE model'});ylabel('f(t)');xlabel('Date(days)')
xlim([1 129]);title('(c)','Position',[-15 18800])
set(gca,'FontName','Times New Roman','FontSize',14)
grid on
grid minor

subplot(425)
FSr=[2.435e7, 2.4349994e7, 2.4349988e7, 2.434998e7, 2.4349974e7, 2.4349966e7, 2.4349958e7, 2.434995e7, 2.434994e7, 2.4349932e7, 2.434992e7, 2.4349908e7, 2.4349892e7, 2.4349874e7, 2.434985e7, 2.4349822e7, 2.4349784e7, 2.4349736e7, 2.4349674e7, 2.4349594e7, 2.434949e7, 2.4349354e7, 2.4349176e7, 2.4348946e7, 2.4348644e7, 2.4348252e7, 2.4347742e7, 2.434708e7, 2.4346222e7, 2.4345112e7, 2.434368e7, 2.4341846e7, 2.43395e7, 2.4336518e7, 2.4332762e7, 2.4328074e7, 2.4322276e7, 2.4315204e7, 2.4306706e7, 2.4296662e7, 2.428501e7, 2.4271756e7, 2.4256986e7, 2.4240852e7, 2.4223582e7, 2.4205468e7, 2.4186806e7, 2.4167894e7, 2.414904e7, 2.4130538e7, 2.411263e7, 2.4095516e7, 2.4079368e7, 2.4064316e7, 2.405046e7, 2.4037858e7, 2.4026544e7, 2.401652e7, 2.4007776e7, 2.4000286e7, 2.3994002e7, 2.3988862e7, 2.398479e7, 2.3981704e7, 2.397952e7, 2.3978094e7, 2.3977222e7, 2.3976636e7, 2.3976008e7, 2.3974952e7, 2.3973404e7, 2.3971504e7, 2.3969316e7, 2.3966892e7, 2.3964286e7, 2.3961536e7, 2.395868e7, 2.3955738e7, 2.3952744e7, 2.3949728e7, 2.3946718e7, 2.3943734e7, 2.3940792e7, 2.3937902e7, 2.3935068e7, 2.3932286e7, 2.3929552e7, 2.3926872e7, 2.3924266e7, 2.3921732e7, 2.3919276e7, 2.39169e7, 2.3914604e7, 2.3912386e7, 2.3910244e7, 2.3908178e7, 2.390619e7, 2.3904278e7, 2.3902444e7, 2.3900686e7, 2.3899002e7, 2.3897392e7, 2.389585e7, 2.3894372e7, 2.389296e7, 2.389161e7, 2.3890322e7, 2.3889094e7, 2.3887926e7, 2.3886816e7, 2.3885758e7, 2.3884754e7, 2.38838e7, 2.3882894e7, 2.3882034e7, 2.3881218e7, 2.3880446e7, 2.3879714e7, 2.3879022e7, 2.3878368e7, 2.3877752e7, 2.3877168e7, 2.3876618e7, 2.3876098e7, 2.3875608e7, 2.3875146e7, 2.387471e7, 2.3874302e7, 2.3873916e7, 2.3873554e7];
aa=1:129;
fill([aa fliplr(aa)],[UconfSr' fliplr(LconfSr')],fillcolor,'edgealpha',0,'facealpha',0.7);%[]中必须都为行向量
hold on
plot(1:0.1:129,x(:,5),'b-',1:129,FSr(2:130),'r--','LineWidth',2)
legend({'95% CI','Mechanistic model','UDE model'});ylabel('Sr(t)');xlabel('Date(days)')
title('(d)','Position',[-12 24500000])
set(gca,'XTick',[ 6 37 67 98 128]);xlim([1 129]);
%set(gca,'XTickLabel',{'','2018.01.31','2020.02.31'});
set(gca,'XTickLabel',{'2022/3/1','4/1','5/1','6/1','7/1'});
set(gca,'FontName','Times New Roman','FontSize',14)
grid on
grid minor

subplot(426)
aa=1:129;
FSe=[400000.0, 400004.25, 400009.2, 400014.28, 400019.2, 400023.62, 400027.3, 400030.0, 400031.47, 400031.4, 400029.44, 400025.2, 400018.03, 400007.22, 399991.7, 399970.06, 399940.5, 399900.6, 399847.1, 399775.72, 399680.84, 399555.06, 399388.6, 399168.47, 398878.12, 398495.7, 397992.53, 397332.66, 396469.88, 395345.47, 393887.8, 392012.34, 389615.47, 386585.4, 382807.28, 378171.5, 372583.5, 366004.28, 358459.16, 350039.84, 340947.84, 331434.53, 321774.47, 312261.53, 303144.53, 294598.66, 286732.03, 279582.6, 273118.25, 267248.5, 261877.39, 256890.36, 252180.3, 247649.95, 243211.9, 238793.95, 234342.7, 229812.58, 225171.38, 220400.42, 215494.48, 210461.97, 205324.44, 200102.56, 194818.16, 189548.81, 184434.64, 179678.25, 175544.8, 172355.28, 170108.39, 168605.72, 167738.28, 167408.72, 167531.38, 168032.27, 168850.02, 169938.89, 171247.89, 172729.8, 174343.92, 176056.1, 177838.67, 179670.53, 181537.06, 183430.17, 185348.1, 187275.05, 189192.45, 191090.95, 192962.84, 194802.14, 196604.55, 198367.45, 200089.56, 201767.94, 203399.11, 204980.42, 206510.0, 207986.8, 209410.53, 210781.75, 212101.81, 213372.88, 214595.58, 215768.95, 216893.33, 217969.27, 218997.56, 219979.28, 220915.72, 221808.45, 222659.25, 223470.02, 224241.61, 224975.02, 225671.34, 226331.8, 226957.62, 227550.16, 228110.81, 228641.1, 229142.55, 229616.81, 230064.86, 230487.67, 230886.4, 231262.22, 231616.23, 231949.56];

fill([aa fliplr(aa)],[UconfSe' fliplr(LconfSe')],fillcolor,'edgealpha',0,'facealpha',0.7);%[]中必须都为行向量
hold on 
plot(1:0.1:129,x(:,1),'b-','LineWidth',2)
plot(aa,FSe(2:end),'r--','LineWidth',2)

legend({'95% CI','Mechanistic model','UDE model'});ylabel('Se(t)');xlabel('Date(days)')
title('(e)','Position',[-15 478750])
set(gca,'XTick',[6 37 67 98 128]);xlim([1 129]);
%set(gca,'XTickLabel',{'','2018.01.31','2020.02.31'});
set(gca,'XTickLabel',{'2022/3/1','4/1','5/1','6/1','7/1'});
set(gca,'FontName','Times New Roman','FontSize',14)
grid on
grid minor

subplot(427)

FI=[2.0, 2.631351, 3.4619992, 4.5548487, 5.992666, 7.8843465, 10.37304, 13.647365, 17.954815, 23.622093, 31.076788, 40.8847, 53.785473, 70.75438, 93.07708, 122.427635, 161.02994, 211.78938, 278.49066, 366.1647, 481.36005, 632.55396, 830.978, 1091.2323, 1431.9546, 1877.5251, 2459.5264, 3217.2227, 4200.6685, 5472.859, 7109.494, 9198.778, 11847.554, 15168.209, 19272.496, 24261.812, 30215.229, 37145.297, 44988.445, 53602.707, 62719.664, 72012.27, 81125.14, 89678.25, 97339.26, 103862.664, 109083.875, 112922.38, 115381.7, 116539.016, 116498.38, 115393.9, 113362.09, 110539.74, 107063.81, 103066.85, 98662.66, 93954.63, 89035.75, 83988.5, 78884.84, 73786.37, 68743.95, 63800.73, 58994.37, 54357.68, 49918.77, 45701.145, 41723.637, 38000.383, 34545.96, 31361.777, 28441.496, 25775.61, 23351.469, 21153.29, 19162.375, 17360.953, 15733.14, 14263.752, 12938.45, 11743.748, 10666.997, 9696.397, 8821.007, 8030.7134, 7316.2725, 6670.0005, 6085.411, 5556.375, 5077.273, 4642.9966, 4248.9443, 3891.0337, 3565.618, 3269.4875, 2999.9377, 2754.4685, 2530.7722, 2326.7327, 2140.4246, 1970.1189, 1814.2728, 1671.5374, 1540.7214, 1420.7815, 1310.7849, 1209.8634, 1117.2123, 1032.091, 953.8239, 881.7989, 815.4686, 754.3417, 697.9822, 646.0167, 598.0972, 553.8971, 513.1118, 475.45807, 440.67484, 408.52176, 378.78137, 351.2572, 325.78485, 302.2108, 280.38843, 260.18164, 241.46524, 224.12381];

fill([aa fliplr(aa)],[UconfI' fliplr(LconfI')],fillcolor,'edgealpha',0,'facealpha',0.7);
hold on 
plot(1:0.1:129,x(:,2),'b-','LineWidth',2)
plot(aa,FI(2:end),'r--','LineWidth',2)

legend({'95% CI','Mechanistic model','UDE model'});ylabel('I(t)');xlabel('Date(days)')
title('(f)','Position',[-12 150000])
set(gca,'XTick',[ 6 37 67 98 128]);xlim([1 129]);
%set(gca,'XTickLabel',{'','2018.01.31','2020.02.31'});
set(gca,'XTickLabel',{'2022/3/1','4/1','5/1','6/1','7/1'});
set(gca,'FontName','Times New Roman','FontSize',14)
grid on
grid minor

subplot(428)

FR=[0.0, 0.46025187, 1.0657938, 1.862485, 2.9106703, 4.2897425, 6.1040783, 8.491217, 11.631641, 15.763627, 21.199078, 28.350796, 37.75855, 50.134346, 66.4172, 87.830475, 116.00061, 153.05511, 201.76823, 265.83502, 350.07794, 460.75616, 606.1954, 797.28955, 1048.0424, 1376.8942, 1808.1417, 2372.5076, 3109.759, 4071.691, 5323.492, 6945.352, 9039.386, 11729.35, 15159.636, 19495.12, 24927.547, 31649.465, 39848.188, 49696.22, 61323.93, 74798.086, 90116.69, 107210.19, 125934.92, 146072.17, 167380.56, 189601.88, 212461.11, 235675.06, 258996.33, 282201.22, 305091.34, 327494.78, 349266.16, 370282.72, 390453.5, 409715.2, 428018.44, 445327.72, 461621.4, 476891.25, 491143.53, 504394.25, 516669.78, 528002.06, 538427.5, 547987.0, 556725.7, 564693.8, 571943.7, 578530.2, 584506.9, 589924.94, 594833.6, 599279.75, 603308.3, 606962.8, 610276.6, 613279.9, 616001.06, 618467.06, 620702.9, 622731.94, 624575.8, 626254.56, 627786.4, 629184.25, 630459.25, 631622.75, 632685.2, 633656.2, 634544.4, 635357.6, 636102.94, 636786.56, 637413.6, 637989.0, 638517.2, 639002.4, 639448.4, 639858.7, 640236.56, 640584.8, 640905.94, 641202.1, 641475.25, 641727.25, 641959.8, 642174.5, 642372.8, 642556.2, 642725.75, 642882.6, 643027.94, 643162.4, 643286.8, 643401.94, 643508.56, 643607.3, 643698.8, 643783.56, 643862.2, 643935.2, 644002.8, 644065.6, 644123.8, 644177.9, 644228.0, 644274.56];
fill([aa fliplr(aa)],[UconfR' fliplr(LconfR')],fillcolor,'edgealpha',0,'facealpha',0.7);
hold on 
plot(1:0.1:129,x(:,3),'b-','LineWidth',2)
plot(aa,FR(2:end),'r--','LineWidth',2)

legend({'95% CI','Mechanistic model','UDE model'},'Location','northwest');ylabel('R(t)');xlabel('Date(days)')
set(gca,'XTick',[ 6 37 67 98 128]);xlim([1 129]);
%set(gca,'XTickLabel',{'','2018.01.31','2020.02.31'});
set(gca,'XTickLabel',{'2022/3/1','4/1','5/1','6/1','7/1'});
title('(g)','Position',[-15 788800])
set(gca,'FontName','Times New Roman','FontSize',14)
grid on
grid minor



%%
%交叉验证图5，表2
r0=[ 2.2 5 8.2 ];%8.2
%r0=[2 2.445 3.5];%r0=beta/gamma=0.489/0.2,[2,5]
k=[ 0.15 0.2812  0.4 ];

n=200;m=length(1:0.1:n);%解方程长度
nn1=120;mm1=length(1:0.1:nn1);%第一列画图长度
nn2=100;mm2=length(1:0.1:nn2);%第二列画图长度
nn3=80;mm3=length(1:0.1:nn3);%第三列画图长度
nnd=[nn1 nn2 nn3];mmd=[mm1 mm2 mm3];
dataSr=zeros(m,9);dataSe=zeros(m,9);dataI=zeros(m,9);dataR=zeros(m,9);dataRt=zeros(m,9);dataNI=zeros(m-10,9);
name={'(a) Wild type (R_0=2.2)','(b) Delta variant (R_0=5)','(c) Omicron variant (R_0=8.2)','(d)','(e)','(f)','(g)','(h)',...
    '(i)','(j)','(k)','(l)','(m)','(n)','(o)'};

name1={'k=0.15','k=0.0.2812','k=0.4','Model (1)','R_t=1'};
Model1=zeros(m,1);Model1NI=zeros(m-10,1);MmdatacumIrate=zeros(1,9);M1datacumIrate=zeros(1,3);M1dataRt=zeros(m,3);

figure  %交叉验证图5
for i=1:3
   gamma=1/5;par1=bestpar1(:,1);
    par1(1)=r0(i)*gamma;%par1： beta        k      I0    Se0
   
    for j=1:3
        par1(2)=k(j);
        X0=[par1(4),par1(3),par2(2),5,par2(1)-par1(4)];
        [t,x]=ode45(@SIRmodel,1:0.1:n,X0,[],par1,par2);
        dataSr(:,j+3*(i-1))=x(:,5);dataSe(:,j+3*(i-1))=x(:,1);dataI(:,j+3*(i-1))=x(:,2);
        dataR(:,j+3*(i-1))=x(:,3);dataRt(:,j+3*(i-1))=5*par1(1).*x(:,1)./(x(:,1)+x(:,2)+x(:,3));dataNI(:,j+3*(i-1))=x(11:end,4)-x(1:m-10,4);
        for z=1:(m-1)
            if floor(dataRt(z,j+3*(i-1)))==1 && floor(dataRt(z+1,j+3*(i-1)))==0
                MmdatacumIrate(1,j+3*(i-1))=(x(z,4)/Sr0_total)*100;%Rt=1时Mechanism model的累计感染率%
            end
        end
    end
   %Model 1
     XX0=[par2(1),par1(3),par2(2),5];
     [t,xx]=ode45(@SIRmodel1,1:0.1:n,XX0,[],par1,par2);%Model 1 未用到par1(2)，所以这里直接放par1无影响
     Model1=[Model1 xx(:,1:3)];Model1NI=[Model1NI (xx(11:end,4)-xx(1:m-10,4))];M1dataRt(:,i)=5*par1(1).*xx(:,1)./(xx(:,1)+xx(:,2)+xx(:,3));
     for z=1:(m-1)
            if floor(M1dataRt(z,i))==1 && floor(M1dataRt(z+1,i))==0
                M1datacumIrate(1,i)=(xx(z,4)/Sr0_total)*100;%Rt=1时Model1的累计感染率%
            end
     end

    nn=nnd(i);mm=mmd(i);
subplot(5,3,1+i-1)
    plot(1:0.1:nn,dataSr(1:mm,(1+3*(i-1)):(3+3*(i-1))),1:0.1:nn,xx(1:mm,1),'k-','LineWidth',1.5)
%legend(name1,'Location','northwest');
ylabel('Sr(t)');xlim([1 nn]);%xlabel('t')
title(name{1+i-1})
set(gca,'FontName','Times New Roman','FontSize',14)
 grid on
grid minor

subplot(5,3,4+i-1)
    plot(1:0.1:nn,dataSe(1:mm,(1+3*(i-1)):(3+3*(i-1))),'LineWidth',1.5)
%legend(name1,'Location','northwest');
ylabel('Se(t)');xlim([1 nn]);%xlabel('t')
title(name{4+i-1})
set(gca,'FontName','Times New Roman','FontSize',14)
 grid on
grid minor

subplot(5,3,7+i-1)
    plot(1:0.1:nn,dataI(1:mm,(1+3*(i-1)):(3+3*(i-1))),1:0.1:nn,xx(1:mm,2),'k-','LineWidth',1.5)
%legend(name1,'Location','northwest');
ylabel('I(t)');xlim([1 nn]);%xlabel('t')
title(name{7+i-1})
set(gca,'FontName','Times New Roman','FontSize',14)
grid on
grid minor


subplot(5,3,10+i-1)
    plot(1:0.1:nn,dataR(1:mm,(1+3*(i-1)):(3+3*(i-1))),1:0.1:nn,xx(1:mm,3),'k-','LineWidth',1.5)
%legend(name1,'Location','northwest');
ylabel('R(t)');xlim([1 nn]);%xlabel('t')
title(name{10+i-1})
set(gca,'FontName','Times New Roman','FontSize',14)
grid on
grid minor

subplot(5,3,13+i-1)
    plot(1:0.1:nn,dataRt(1:mm,(1+3*(i-1)):(3+3*(i-1))),1:0.1:nn,5*par1(1).*xx(1:mm,1)./(xx(1:mm,1)+xx(1:mm,2)+xx(1:mm,3)),'k-','LineWidth',1.5)
%legend(name1,'Location','northwest');
ylabel('R_t');xlim([1 nn]);xlabel('Time')
title(name{13+i-1})
set(gca,'FontName','Times New Roman','FontSize',14)
line([1 nn],[1 1],'linestyle','--','LineWidth',1,'color',[0.49,0.18,0.56])
if 13+i-1==15
    legend(name1);
end
grid on
grid minor
end
%放大图
axes('Position',[0.15 0.48 0.07 0.07]);
plot(1:0.1:nnd(1),dataI(1:mmd(1),1:3),'LineWidth',1.5);
grid on
grid minor
axes('Position',[0.15 0.305 0.07 0.07]);
plot(1:0.1:nnd(1),dataR(1:mmd(1),1:3),'LineWidth',1.5);
grid on
grid minor

%表2
Model1(:,1)=[];Model1NI(:,1)=[];
%Mechanistic model
final_size=dataR(end,:);attack_rate=final_size./Sr0_total;
peak_size=max(dataI,[],1);peak_size_rate=max(dataI,[],1)./Sr0_total;
peak_sizeNI=max(dataNI,[],1);peak_sizeNI_rate=max(dataNI,[],1)./Sr0_total;%返回每列最大值
Sr_inf=dataSr(end,:);

Mmfinal_size=[0 r0; k' final_size(1:3)' final_size(4:6)' final_size(7:9)' ]'./10^6
Mmattack_rate=[0 r0; k' attack_rate(1:3)' attack_rate(4:6)' attack_rate(7:9)' ]'.*10^2
Mmpeak_size=[0 r0; k' peak_size(1:3)' peak_size(4:6)' peak_size(7:9)' ]'./10^6
Mmpeak_size_rate=[0 r0; k' peak_size_rate(1:3)' peak_size_rate(4:6)' peak_size_rate(7:9)' ]'.*10^2
Mmpeak_sizeNI=[0 r0; k' peak_sizeNI(1:3)' peak_sizeNI(4:6)' peak_sizeNI(7:9)' ]'./10^6
Mmpeak_sizeNI_rate=[0 r0; k' peak_sizeNI_rate(1:3)' peak_sizeNI_rate(4:6)' peak_sizeNI_rate(7:9)' ]'.*10^2
MmSr_inf=[0 r0; k' Sr_inf(1:3)' Sr_inf(4:6)' Sr_inf(7:9)' ]'./10^6
MmdatacumIrate=[0 r0; k' MmdatacumIrate(1:3)' MmdatacumIrate(4:6)' MmdatacumIrate(7:9)' ]'
%Model1
M1final_size=Model1(end,[3 6 9])./10^6
M1attack_rate=(Model1(end,[3 6 9])./Sr0_total).*10^2
M1peak_size=max(Model1(:,[2 5 8]),[],1)./10^6%返回每列最大值
M1peak_size_rate=(max(Model1(:,[2 5 8]),[],1)./Sr0_total).*10^2
M1peak_sizeNI=max(Model1NI,[],1)./10^6.
M1peak_sizeNI_rate=(max(Model1NI,[],1)./Sr0_total).*10^2
M1Sr_inf=Model1(end,[1 4 7])./10^6



%%
%图6  对比图
r0_6=round(1.5:0.1:10,2);k=[ 0.15 0.2812  0.4 ];
%Mechanistic model
dataI_6=zeros(m,length(r0_6)*3);dataR_6=zeros(m,length(r0_6)*3);
dataNI_6=zeros(m-10,length(r0_6)*3);dataRt_6=zeros(m,length(r0_6)*3);
MmdatacumIrate_6=zeros(1,length(r0_6)*3);
%Model 1
M1dataI_6=zeros(m,length(r0_6));M1dataR_6=zeros(m,length(r0_6));
M1dataNI_6=zeros(m-10,length(r0_6));M1dataRt_6=zeros(m,length(r0_6));
M1datacumIrate_6=zeros(1,length(r0_6));

for i=1:length(r0_6)


    gamma=1/5;par1=bestpar1(:,1);
    par1(1)=r0_6(i)*gamma;

    %Mechanistic model
    for j=1:3
        par1(2)=k(j);
        X0=[par1(4),par1(3),par2(2),5,par2(1)-par1(4)];
        [t,x]=ode45(@SIRmodel,1:0.1:n,X0,[],par1,par2);
        dataI_6(:,i+length(r0_6)*(j-1))=x(:,2);dataR_6(:,i+length(r0_6)*(j-1))=x(:,3);
        dataRt_6(:,i+length(r0_6)*(j-1))=5*par1(1).*x(:,1)./(x(:,1)+x(:,2)+x(:,3));
        dataNI_6(:,i+length(r0_6)*(j-1))=x(11:end,4)-x(1:m-10,4);
        for z=1:(m-1)
            if floor(dataRt_6(z,i+length(r0_6)*(j-1)))==1 && floor(dataRt_6(z+1,i+length(r0_6)*(j-1)))==0
                MmdatacumIrate_6(1,i+length(r0_6)*(j-1))=(x(z,4)/Sr0_total)*100;%Rt=1时Mechanism model的累计感染率%
            end
        end
    end
   %Model 1
     XX0=[par2(1),par1(3),par2(2),5];
     [t,xx]=ode45(@SIRmodel1,1:0.1:n,XX0,[],par1,par2);%Model 1 未用到par1(2)，所以这里直接放par1无影响
     M1dataI_6(:,i)=xx(:,2);M1dataR_6(:,i)=xx(:,3);
     M1dataRt_6(:,i)=5*par1(1).*xx(:,1)./(xx(:,1)+xx(:,2)+xx(:,3));
     M1dataNI_6(:,i)=xx(11:end,4)-xx(1:m-10,4);
     
     for z=1:(m-1)
            if floor(M1dataRt_6(z,i))==1 && floor(M1dataRt_6(z+1,i))==0
                M1datacumIrate_6(1,i)=(xx(z,4)/Sr0_total)*100;%Rt=1时Model1的累计感染率%
            end
     end

end
Mmattack_rate_6=dataR_6(end,:)./Sr0_total.*100;
Mmpeak_size_rate_6=max(dataI_6,[],1)./Sr0_total.*100;
Mmpeak_sizeNI_rate_6=max(dataNI_6,[],1)./Sr0_total.*100;%返回每列最大值

M1attack_rate_6=M1dataR_6(end,:)./Sr0_total.*100;
M1peak_size_rate_6=max(M1dataI_6,[],1)./Sr0_total.*100;
M1peak_sizeNI_rate_6=max(M1dataNI_6,[],1)./Sr0_total.*100;%返回每列最大值


figure
subplot(221)
plot([8.2*0.9 8.2*1.1],[82.4 82.4],'Color',[0.72,0.27,1],'LineWidth',4)
hold on 
plot([2.43 5.11],[6.6 6.6],'Color',[0.39,0.83,0.07],'LineWidth',4)
plot([1.4 3.9],[0.8 0.8],'Color',[1 0.07 0.65],'LineWidth',4)

plot(r0_6,Mmattack_rate_6(1:length(r0_6)),'Color',[0 0.45 0.74],'LineWidth',1.5)
plot(r0_6, Mmattack_rate_6(length(r0_6)+1:length(r0_6)*2),'Color',[0.85 0.33 0.1],'LineWidth',1.5)
plot(r0_6,Mmattack_rate_6(length(r0_6)*2+1:length(r0_6)*3)','Color',[0.93 0.69 0.13],'LineWidth',1.5)
plot(r0_6,M1attack_rate_6','k-','LineWidth',1.5)
xlim([1.5 10]);ylim([0 100]);xlabel('R_0');ylabel('Rate/%');title('(a) Attack rate')
set(gca,'XTick',[ 1.5 3 5 7  8.5 10],'FontName','Times New Roman','FontSize',14);
grid on
grid minor
hold off

subplot(222)
r0_6_0=round(1.5:0.5:10,2);pos=zeros(1,length(r0_6_0));%以0.5为步长画图
for i=1:length(r0_6_0)
pos(i)=find(r0_6==r0_6_0(i));%
end
plot(r0_6_0,[MmdatacumIrate_6(pos)' MmdatacumIrate_6(pos+length(r0_6))' MmdatacumIrate_6(pos+length(r0_6)*2)'],r0_6_0,M1datacumIrate_6(pos),'k-','LineWidth',1.5)
xlim([1.5 10]);ylim([0 100]);xlabel('R_0');ylabel('Rate/%');title('(b) Herd immunity')
set(gca,'XTick',[ 1.5 3 5 7  8.5 10],'FontName','Times New Roman','FontSize',14);
grid on
grid minor

subplot(223)
plot(r0_6,[Mmpeak_size_rate_6(1:length(r0_6))' Mmpeak_size_rate_6(length(r0_6)+1:length(r0_6)*2)' Mmpeak_size_rate_6(length(r0_6)*2+1:length(r0_6)*3)'],r0_6,M1peak_size_rate_6,'k-','LineWidth',1.5)
xlim([1.5 10]);ylim([0 70]);xlabel('R_0');ylabel('Rate/%');title('(c) Peak of infection class')
set(gca,'XTick',[ 1.5 3 5 7  8.5 10],'FontName','Times New Roman','FontSize',14);
grid on
grid minor


subplot(224)
plot(r0_6,[Mmpeak_sizeNI_rate_6(1:length(r0_6))' Mmpeak_sizeNI_rate_6(length(r0_6)+1:length(r0_6)*2)' Mmpeak_sizeNI_rate_6(length(r0_6)*2+1:length(r0_6)*3)'],r0_6,M1peak_sizeNI_rate_6,'k-','LineWidth',1.5)
xlim([1.5 10]);ylim([0 50]);xlabel('R_0');ylabel('Rate/%');title('(d) Peak daily new infections')
legend('k=0.15','k=0.0.2812','k=0.4','Model (1)','Location','northwest')
set(gca,'XTick',[ 1.5 3 5 7  8.5 10],'FontName','Times New Roman','FontSize',14);
grid on
grid minor

set(gcf,'unit','centimeters','position',[10 10 20 15])%设置窗口的大小，（10,10）起点坐标，20宽，15高





%%
%图3
%新发报道-感染间隔分布拟合图和新发感染数据图
interval=[0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	3	3	3	4	4	5	5	5	5	5	5	5	6	6	6	6	6	6	6	6	6];
interval_x=[0.5,1.5,2.5,3.5,4.5,5.5,6.5];
interval_p=[9,34,29,3,2,7,9]./sum([9,34,29,3,2,7,9])
figure
subplot(121)
bar(interval_x,interval_p,1,'FaceColor',[.73 .80 1])
hold on
gam=[0.150335471	0.178445207	0.203323613	0.224838769	0.243011002	0.257957895	0.269858075	0.278926632	0.285398122	0.28951472	0.291517978	0.29164314	0.29011531	0.28714697	0.282936495	0.277667379	0.271508005	0.264611789	0.257117595	0.249150344	0.24082174	0.232231086	0.223466125	0.214603908	0.205711651	0.196847569	0.188061686	0.179396602	0.170888217	0.162566415	0.154455695	0.14657576	0.13894206	0.131566288	0.124456835	0.117619207	0.111056398	0.104769232	0.09875667	0.093016083	0.087543504	0.082333844	0.077381089	0.072678477	0.068218647	0.063993778	0.059995706	0.056216027	0.052646188	0.049277564	0.046101524	0.043109491	0.040292984	0.037643664	0.035153366	0.032814122	0.030618191	0.028558066	0.026626497	0.024816494	0.023121332];
weib=[0.175349671	0.191960951	0.206166624	0.218251596	0.228429772	0.236872214	0.243722678	0.249106691	0.2531371	0.255917554	0.257544739	0.258109831	0.257699427	0.256396157	0.254279062	0.25142383	0.247902931	0.243785685	0.239138298	0.234023865	0.228502372	0.222630681	0.216462532	0.210048541	0.203436214	0.196669966	0.189791155	0.182838127	0.175846272	0.16884809	0.161873266	0.154948763	0.148098908	0.141345498	0.134707907	0.128203199	0.121846241	0.115649828	0.109624796	0.103780151	0.098123187	0.092659606	0.087393641	0.082328168	0.077464824	0.072804118	0.068345534	0.064087639	0.060028178	0.056164168	0.05249199	0.049007471	0.045705959	0.042582405	0.039631424	0.036847364	0.034224362	0.031756396	0.029437342	0.027261007	0.025221181];
ln=[0.124680153	0.176605865	0.224178608	0.264461744	0.296363621	0.319945284	0.335900508	0.345212634	0.348946238	0.34813023	0.343698903	0.336467921	0.327130291	0.316263013	0.304338793	0.291739543	0.278769864	0.265669579	0.252624892	0.23977806	0.227235625	0.21507532	0.203351803	0.192101389	0.181345923	0.171095932	0.161353188	0.152112769	0.143364725	0.135095403	0.127288506	0.119925936	0.112988452	0.106456188	0.100309055	0.094527045	0.089090471	0.08398013	0.079177437	0.074664506	0.070424208	0.066440207	0.062696974	0.059179788	0.055874726	0.05276865	0.049849176	0.047104656	0.044524138	0.042097344	0.039814628	0.037666952	0.035645844	0.033743373	0.031952115	0.030265122	0.028675892	0.027178348	0.0257668	0.02443593	0.023180766];
p1=plot(0.5:0.1:6.5,gam,'--','LineWidth',2,'Color',[.929 .694 .125])
p2=plot(0.5:0.1:6.5,weib,'-.','LineWidth',2,'Color',[.47 .67 .19])
p3=plot(0.5:0.1:6.5,ln,'-','LineWidth',2,'Color',[.72 .27 1])
%Weibull AIC=340.9162", "gamma AIC=336.3352", "lognormal AIC=335.9116
legend([p1 p2 p3],{'Gamma AIC=336.3352','Weibull AIC=340.9162','Lognormal AIC=335.9116'})
set(gca,'XTick',[ 0.5 1.5 2.5 3.5 4.5 5.5 6.5]);%xlim([1 129]);
set(gca,'XTickLabel',{'0.5','1','2','3','4','5','6'});
xlabel('Interval length (days)');ylabel('Density');title('(b)');
set(gca,'FontName','Times New Roman','FontSize',14)
grid on
grid minor

%
reported_case=xlsread('reported_case.xlsx');%新发报道数据
%data1(2,:)';%新发感染数据

subplot(122)
plot(6:129,reported_case,'o','LineWidth',2,'MarkerSize',2,'Color',[.47 .67 .19])
hold on
plot(1:129,data1(2,:)','-','LineWidth',2,'Color',[.72 .27 1])
set(gca,'XTick',[ 6 37 67 98 128]);xlim([1 129]);xlabel('Date(days)');ylabel('Cases')
%set(gca,'XTickLabel',{'','2018.01.31','2020.02.31'});
set(gca,'XTickLabel',{'2022/3/1','4/1','5/1','6/1','7/1'});
title('(c)');legend({'Daily reported cases','Daily new infections'})
set(gca,'FontName','Times New Roman','FontSize',14)
grid on
grid minor
set(gcf,'unit','centimeters','position',[10 10 30 11])




