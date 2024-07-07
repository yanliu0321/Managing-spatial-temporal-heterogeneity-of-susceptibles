
function dy=SIRmodel1(i,y,par1,par2)
dy=zeros(4,1);
beta=par1(1);gamma=1/5;
N=par2(1);

dy=[
   -beta*y(1)*y(2)/N;
    beta*y(1)*y(2)/N-gamma*y(2);
    gamma*y(2);
    beta*y(1)*y(2)/N
        ];