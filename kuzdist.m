function d=kuzdist(par1,data)
%Sr0_total=24750000;%总人口
%R0=0;%初始恢复者
%par2=[Sr0_total R0];
% beta=par1(1);
% gam=par1(2);
% K=par1(3);
% I0=par1(4);
% Se0=par1(5);
% ydata=[(1:129)',data1(2,2:end)'];ylabels={'time','newly_infected'};
% zdata=par2';zlabels={'par2'};
% data=struct('ydata',ydata,'zdata',zdata);

par2=data.zdata(:,1);datanew=data.ydata(2:end,2);
X0=[par1(4),par1(3),par2(2),5,par2(1)-par1(4)];
[t,x]=ode45(@SIRmodel,1:129,X0,[],par1,par2);

newly_infected=diff(x(:,4));
d=sqrt(sum((newly_infected-datanew).^2));

