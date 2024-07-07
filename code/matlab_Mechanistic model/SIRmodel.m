function dy=SIRmodel(i,y,par1,par2)
dy=zeros(5,1);
beta=par1(1);gamma=1/5;k=par1(2);%I0=par1(4);Se0=par1(5);
N=par2(1);
Ne=y(1)+y(2)+y(3);
dy=[
   k*beta*y(2)*(N-Ne)/N-beta*y(1)*y(2)/Ne;
    beta*y(1)*y(2)/Ne-gamma*y(2);
    gamma*y(2);
    beta*y(1)*y(2)/Ne;
    -k*beta*y(2)*(N-Ne)/N;
    ];