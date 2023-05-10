%% Lieb lattice edge
clear all; clc;
global Nx;global N0; 
N0 = 1;%N0=0,1表示第一行是C还是AB子晶格
Nx = 41;
sign1 = 1;%1:代表swave, -1代表staggered swave
a = 3.14;

t = 1;
ilaso = 1i*0.1;
Delta0 = 0.1;
mu = -0.0;

Bx =  0.0;
By = -0.0; %Delta=0.5,mu=-0.2

str = strcat('.\Figure\OBC_M5_type1','.jpg');

for kx=0:0.01:3.14
    matrix=[];
    for i=1:1:Nx
       if mod(N0+i,2)==0
           for k=1:1:2
              matrix(posi(i,1)+k,posi(i,2)+k)=2*t*cos(kx);
              matrix(posi(i,2)+k,posi(i,1)+k)=2*t*cos(kx);
              matrix(posi(i,1)+k+2,posi(i,2)+k+2)=-2*t*cos(kx);
              matrix(posi(i,2)+k+2,posi(i,1)+k+2)=-2*t*cos(kx);
           end
          if i>1
              for k=1:1:2
                 matrix(posi(i,1)+k,posi(i-1,1)+k)=t;
                 matrix(posi(i-1,1)+k,posi(i,1)+k)=t;
                 matrix(posi(i,1)+k+2,posi(i-1,1)+k+2)=-t;
                 matrix(posi(i-1,1)+k+2,posi(i,1)+k+2)=-t;
              end
          end
          if i<Nx
              for k=1:1:2
                 matrix(posi(i,1)+k,posi(i+1,1)+k)=t;
                 matrix(posi(i+1,1)+k,posi(i,1)+k)=t;
                 matrix(posi(i,1)+k+2,posi(i+1,1)+k+2)=-t;
                 matrix(posi(i+1,1)+k+2,posi(i,1)+k+2)=-t;
              end
          end
          %s波超导
          matrix(posi(i,1)+1,posi(i,1)+4)=Delta0;
          matrix(posi(i,1)+4,posi(i,1)+1)=Delta0;
          matrix(posi(i,1)+2,posi(i,1)+3)=-Delta0;
          matrix(posi(i,1)+3,posi(i,1)+2)=-Delta0;
          matrix(posi(i,2)+1,posi(i,2)+4)=Delta0;
          matrix(posi(i,2)+4,posi(i,2)+1)=Delta0;
          matrix(posi(i,2)+2,posi(i,2)+3)=-Delta0;
          matrix(posi(i,2)+3,posi(i,2)+2)=-Delta0;
          
          %塞曼场(未完成)
          matrix(posi(i,1)+1,posi(i,1)+2)=matrix(posi(i,1)+1,posi(i,1)+2)+Bx;
          matrix(posi(i,1)+2,posi(i,1)+1)=matrix(posi(i,1)+2,posi(i,1)+1)+Bx;
          matrix(posi(i,1)+3,posi(i,1)+4)=matrix(posi(i,1)+3,posi(i,1)+4)-Bx;
          matrix(posi(i,1)+4,posi(i,1)+3)=matrix(posi(i,1)+4,posi(i,1)+3)-Bx;
          matrix(posi(i,1)+1,posi(i,1)+2)=matrix(posi(i,1)+1,posi(i,1)+2)-1i*By;
          matrix(posi(i,1)+2,posi(i,1)+1)=matrix(posi(i,1)+2,posi(i,1)+1)+1i*By;
          matrix(posi(i,1)+3,posi(i,1)+4)=matrix(posi(i,1)+3,posi(i,1)+4)+1i*By;
          matrix(posi(i,1)+4,posi(i,1)+3)=matrix(posi(i,1)+4,posi(i,1)+3)-1i*By;
          
          %化学势
          matrix(posi(i,1)+1,posi(i,1)+1)=matrix(posi(i,1)+1,posi(i,1)+1)-mu;
          matrix(posi(i,1)+2,posi(i,1)+2)=matrix(posi(i,1)+2,posi(i,1)+2)-mu;
          matrix(posi(i,1)+3,posi(i,1)+3)=matrix(posi(i,1)+3,posi(i,1)+3)+mu;
          matrix(posi(i,1)+4,posi(i,1)+4)=matrix(posi(i,1)+4,posi(i,1)+4)+mu;
          
          matrix(posi(i,2)+1,posi(i,2)+1)=matrix(posi(i,2)+1,posi(i,2)+1)-mu;
          matrix(posi(i,2)+2,posi(i,2)+2)=matrix(posi(i,2)+2,posi(i,2)+2)-mu;
          matrix(posi(i,2)+3,posi(i,2)+3)=matrix(posi(i,2)+3,posi(i,2)+3)+mu;
          matrix(posi(i,2)+4,posi(i,2)+4)=matrix(posi(i,2)+4,posi(i,2)+4)+mu;
       end
       if mod(N0+i,2)==1
               if i>1
                   matrix(posi(i,1)+1,posi(i-1,2)+1)=-ilaso*2*1i*sin(kx);
                   matrix(posi(i,1)+2,posi(i-1,2)+2)=ilaso*2*1i*sin(kx);
                   matrix(posi(i-1,2)+1,posi(i,1)+1)=-ilaso*2*1i*sin(kx);
                   matrix(posi(i-1,2)+2,posi(i,1)+2)=ilaso*2*1i*sin(kx);
                   
                   matrix(posi(i,1)+3,posi(i-1,2)+3)=ilaso*2*1i*sin(-kx);
                   matrix(posi(i,1)+4,posi(i-1,2)+4)=-ilaso*2*1i*sin(-kx);
                   matrix(posi(i-1,2)+3,posi(i,1)+3)=ilaso*2*1i*sin(-kx);
                   matrix(posi(i-1,2)+4,posi(i,1)+4)=-ilaso*2*1i*sin(-kx);
               end
               if i<Nx
                   matrix(posi(i,1)+1,posi(i+1,2)+1)=ilaso*2*1i*sin(kx);
                   matrix(posi(i,1)+2,posi(i+1,2)+2)=-ilaso*2*1i*sin(kx);
                   matrix(posi(i+1,2)+1,posi(i,1)+1)=ilaso*2*1i*sin(kx);
                   matrix(posi(i+1,2)+2,posi(i,1)+2)=-ilaso*2*1i*sin(kx);
                   
                   matrix(posi(i,1)+3,posi(i+1,2)+3)=-ilaso*2*1i*sin(-kx);
                   matrix(posi(i,1)+4,posi(i+1,2)+4)=ilaso*2*1i*sin(-kx);
                   matrix(posi(i+1,2)+3,posi(i,1)+3)=-ilaso*2*1i*sin(-kx);
                   matrix(posi(i+1,2)+4,posi(i,1)+4)=ilaso*2*1i*sin(-kx);
               end
             matrix(posi(i,1)+1,posi(i,1)+4)=Delta0*sign1;
             matrix(posi(i,1)+4,posi(i,1)+1)=Delta0*sign1;
             matrix(posi(i,1)+2,posi(i,1)+3)=-Delta0*sign1;
             matrix(posi(i,1)+3,posi(i,1)+2)=-Delta0*sign1;
             
             %化学势
             matrix(posi(i,1)+1,posi(i,1)+1)=matrix(posi(i,1)+1,posi(i,1)+1)-mu;
             matrix(posi(i,1)+2,posi(i,1)+2)=matrix(posi(i,1)+2,posi(i,1)+2)-mu;
             matrix(posi(i,1)+3,posi(i,1)+3)=matrix(posi(i,1)+3,posi(i,1)+3)+mu;
             matrix(posi(i,1)+4,posi(i,1)+4)=matrix(posi(i,1)+4,posi(i,1)+4)+mu;
             
             %塞曼场
             matrix(posi(i,1)+1,posi(i,1)+2)=matrix(posi(i,1)+1,posi(i,1)+2)+Bx;
             matrix(posi(i,1)+2,posi(i,1)+1)=matrix(posi(i,1)+2,posi(i,1)+1)+Bx;
             matrix(posi(i,1)+3,posi(i,1)+4)=matrix(posi(i,1)+3,posi(i,1)+4)-Bx;
             matrix(posi(i,1)+4,posi(i,1)+3)=matrix(posi(i,1)+4,posi(i,1)+3)-Bx;
             matrix(posi(i,1)+1,posi(i,1)+2)=matrix(posi(i,1)+1,posi(i,1)+2)-1i*By;
             matrix(posi(i,1)+2,posi(i,1)+1)=matrix(posi(i,1)+2,posi(i,1)+1)+1i*By;
             matrix(posi(i,1)+3,posi(i,1)+4)=matrix(posi(i,1)+3,posi(i,1)+4)+1i*By;
             matrix(posi(i,1)+4,posi(i,1)+3)=matrix(posi(i,1)+4,posi(i,1)+3)-1i*By;
       end
    end
    ans(:,round((kx+0)/0.01)+1)=eig(matrix);
end



for (j=1:1:4*(fix(Nx/2)*3+mod(Nx,2)*(N0+1)))
   if (j>=2*(fix(Nx/2)*3+mod(Nx,2)*(N0+1))-1&&j<=2*(fix(Nx/2)*3+mod(Nx,2)*(N0+1))+2) 
       plot(0:0.01:a,ans(j,:),'k','linewidth',1); 
   else
       plot(0:0.01:a,ans(j,:),'k','linewidth',1); 
   end

   hold on;
end
set(gca,'FontName','Times New Roman','FontSize',20);
xlabel('\it{k_x}');
ylabel('\it{E}');
xlim([0,pi]);
ylim([-3,3]);
xticks([0 pi/2 pi]);

ax = gca;
ax.TickLabelInterpreter = 'latex';
xticklabels({'0','$\pi/2$','$\pi$'});
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.5]);
saveas(gcf,str);

function P=posi(i,j)
   global Nx;global N0;
   if mod(N0,2)==0
       if mod(i,2)==1
          P=4*(3*1*fix((i-1)/2)+(j-1));
       elseif mod(i,2)==0
          P=4*(3*1*fix((i-1)/2)+1+(j-1));
       end
   elseif mod(N0,2)==1
       if mod(i,2)==1
          P=4*(3*1*fix((i-1)/2)+(j-1)); 
       elseif mod(i,2)==0
           P=4*(3*1*fix((i-1)/2)+2*1+(j-1));
       end
   end
end