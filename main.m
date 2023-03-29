clear ;
rand('state',sum(100*clock));
Dimension=[10 30 50 100];
Xmin=-100;
Xmax=100;
runs=5;     
fhd=str2func('cec17_func');
Findex= [1,3:30];

for idea=2.4 %%% The parameter c in Alg.7

tic
for dimension=1:1 
  
    D=Dimension(dimension);
    Max_FEs=10000*D;
for i=1:29
    func_num=Findex(i);

    for jj=1:runs
        FEs=0;
        total=200;
         mu=100;
         vara=[func_num,-100,100,zeros(1,100),0,0];
         VRmin=repmat(Xmin,total,D);
         VRmax=repmat(Xmax,total,D);
         posinitial=VRmin+(VRmax-VRmin).*rand(total,D);
         e=feval(fhd,posinitial',vara(:));
         FEs=FEs+total;
         weights=log(mu+1/2)-log(1:mu)';
         weights=weights/sum(weights);
         [a1,a2]=sort(e);
         bestval=a1(1);
          bestvec=posinitial(a2(1),:);
         pos=posinitial(a2(1:mu),:);
         meanval=mean(pos);
         stdval=std(pos);
         for k=1:total
          pos(k,:)=meanval+stdval.*randn(1,D);
         end
         for k=1:total
           for j=1:D
             if pos(k,j)>100
                 pos(k,j)=meanval(j)+stdval(j).*randn;
             elseif pos(k,j)<-100
                 pos(k,j)=meanval(j)+stdval(j).*randn;
             end
           end
         end
         cc1=0;
         kk = 1;

xiaomei = [];
xiaomei(1)=bestval;
        while USE_DQN(kk-1,xiaomei,idea) == 0
             e=feval(fhd,pos',vara(:));
             FEs=FEs+total;

    if FEs-0.01*Max_FEs<total
        processvalue(1)=a1(1);
    elseif FEs-0.02*Max_FEs<total
         processvalue(2)=a1(1);
    elseif FEs-0.03*Max_FEs<total
         processvalue(3)=a1(1);
         elseif FEs-0.05*Max_FEs<total
         processvalue(4)=a1(1);
         elseif FEs-0.1*Max_FEs<total
         processvalue(5)=a1(1);
         elseif FEs-0.2*Max_FEs<total
         processvalue(6)=a1(1);
         elseif FEs-0.3*Max_FEs<total
         processvalue(7)=a1(1);
         elseif FEs-0.4*Max_FEs<total
         processvalue(8)=a1(1);
    elseif FEs-0.5*Max_FEs<total
         processvalue(9)=a1(1);
         elseif FEs-0.6*Max_FEs<total
         processvalue(10)=a1(1);
         elseif FEs-0.7*Max_FEs<total
         processvalue(11)=a1(1);
         elseif FEs-0.8*Max_FEs<total
         processvalue(12)=a1(1);
         elseif FEs-0.9*Max_FEs<total
         processvalue(13)=a1(1);
         elseif FEs-1*Max_FEs<total
         processvalue(14)=a1(1);
    end
             
             
             [a1,a2]=sort(e);
             if a1(1)<bestval
                 bestval=a1(1);
                 bestvec=pos(a2(1),:);
             end
             xiaomei(kk+1)=bestval;
             newpos=pos(a2(1:mu),:);
             meanval=(newpos(:,1:D)'*weights)';
             stdval=1*std(newpos);
             FV(kk)=a1(1);
             if kk>30 
               if mod(kk,20)==0
                  [aa1,aa2]=min(FV);
                  if aa2<kk-20
                     cc1=1;
                  end
               end
             end
             for k=1:total
                 if cc1==1      %kk>300
                    a=0.96*randn(1,D);
                 else
                     a=randn(1,D);
                 end
                 pos(k,:)=meanval+stdval.*a;
             end
              for k=1:total
                 for j=1:D
                    if pos(k,j)>100
                       pos(k,j)=mod(pos(k,j),100);
                    elseif pos(k,j)<-100
                       pos(k,j)=mod(pos(k,j),-100);
                    end
                 end
              end
              kk = kk+1;
        end
%         disp(kk-1)
% zht (i,jj) = kk;
        previousbest=a1(1);
              if D<=30
                  Times=2;
             else
                  Times=1;
              end
             arfitnessbest=bestval.*ones(1,Times);
             xvalbest=repmat(bestvec',1,Times);
              N=D;
              index=1;
             for kkk=1:Times
             sigma=0.2;
             stopfitness=1e-8;
             if D<=30
                   stopeval=Max_FEs/4;
             else
                  stopeval=Max_FEs/2;
             end
             lambda=floor(3*log(N))+80;
             mu=lambda/2;
             weights=log(mu+1/2)-log(1:mu)';
             mu=floor(mu);
             weights=weights/sum(weights);
             mueff=sum(weights)^2/sum(weights.^2);
             cc=(4+mueff/N)/(N+4+2*mueff/N);
             cs=(mueff+2)/(N+mueff+5);
             c1=2/((N+1.3)^2+mueff);
             cmu=2*(mueff-2+1/mueff)/((N+2)^2+2*mueff/2);
             damps=1+2*max(0,sqrt((mueff-1)/(N+1))-1)+cs;           
              pc=zeros(N,1);
              ps=zeros(N,1);
              B=eye(N);
              DD=eye(N);
              C=B*DD*(B*DD)';
              eigenval=0;
             chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));                 %unknown
             counteval=0;
             xmean=bestvec';
             while counteval<stopeval
                 for k=1:lambda
                     arz(:,k)=randn(N,1);
                     arxx(:,k)=xmean+1*sigma*B*DD*arz(:,k);
                  for jjj=1:N
                      if real(arxx(jjj,k))>100
                         arxx(jjj,k)=mod(real(arxx(jjj,k)),100);
                    elseif real(arxx(jjj,k))<-100
                         arxx(jjj,k)=mod(real(arxx(jjj,k)),-100);
                      end
                  end
                     arfitness(k)=feval(fhd,arxx(:,k),vara(:));
                     counteval=counteval+1;
                      FEs=FEs+1;
                 end
                 
                  if FEs-0.01*Max_FEs<total
        processvalue(1)=a1(1);
    elseif FEs-0.02*Max_FEs<total
         processvalue(2)=a1(1);
    elseif FEs-0.03*Max_FEs<total
         processvalue(3)=a1(1);
         elseif FEs-0.05*Max_FEs<total
         processvalue(4)=a1(1);
         elseif FEs-0.1*Max_FEs<total
         processvalue(5)=a1(1);
         elseif FEs-0.2*Max_FEs<total
         processvalue(6)=a1(1);
         elseif FEs-0.3*Max_FEs<total
         processvalue(7)=a1(1);
         elseif FEs-0.4*Max_FEs<total
         processvalue(8)=a1(1);
    elseif FEs-0.5*Max_FEs<total
         processvalue(9)=a1(1);
         elseif FEs-0.6*Max_FEs<total
         processvalue(10)=a1(1);
         elseif FEs-0.7*Max_FEs<total
         processvalue(11)=a1(1);
         elseif FEs-0.8*Max_FEs<total
         processvalue(12)=a1(1);
         elseif FEs-0.9*Max_FEs<total
         processvalue(13)=a1(1);
         elseif FEs-1*Max_FEs<total
         processvalue(14)=a1(1);
                  end
    
                  
                     [arfitness, arindex]=sort(arfitness);
                     xval=arxx(:,arindex(1)); 
                     
                     if abs(arfitness(1)-previousbest)<1*10^(-11)
                         break;
                     else
                         previousbest=arfitness(1);
                     end
                     
%                      if arfitnessbest(kkk)>arfitnessbest(1)
%                          arfitnessbest(kkk)=arfitnessbest(1);
%                      end                 
                     if arfitnessbest(kkk)>arfitness(1)
                        arfitnessbest(kkk)=arfitness(1);
                        xvalbest(:,kkk)=arxx(:,arindex(1));
                     end
                     xmean=arxx(:,arindex(1:mu))*weights;
                     zmean=arz(:,arindex(1:mu))*weights;
                     ps=(1-cs)*ps+(sqrt(cs*(2-cs)*mueff))*(B*zmean);
                     hsig=norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN<1.4+2/(N+1);    
                     pc=(1-cc)*pc+hsig*sqrt(cc*(2-cc)*mueff)*(B*DD*zmean);
                     C=(1-c1-cmu)*C+c1*(pc*pc'+(1-hsig)*cc*(2-cc)*C)+cmu*(B*DD*arz(:,arindex(1:mu)))*diag(weights)*(B*DD*arz(:,arindex(1:mu)))';
                     sigma=sigma*exp((cs/damps)*(norm(ps)/chiN-1));
                     xx(counteval/lambda)=sigma;
                     if counteval-eigenval>lambda/(cmu)/N/10
                         eigenval=counteval;
                         C=triu(C)+triu(C,1)';
                         [B,DD]=eig(C);
                          DD=diag(sqrt(diag(DD)));
                     end
                     
                     if arfitness(1)==arfitness(ceil(0.7*lambda))
                        sigma=sigma*exp(0.2+cs/damps);
                        xx(counteval/lambda)=sigma;
                     end  

             end
             end
             xiuxiuxiu=index+kk;

             
             if D<=30
               total=200;
               mu=160;
             elseif D==50
               total=450;
               mu=360;
             else
                  total=600;
               mu=480;
             end
            if D>=50
                 if FEs<=0.3*Max_FEs
                     total=total+200;
                     mu=mu+160;
                 end
             end
             weights=log(mu+1/2)-log(1:mu)';
             weights=weights/sum(weights);
             if D<=30
                 ppp1=std(xvalbest');
              ppp2=sort(std(xvalbest'));
                if ppp2(1)>0.2
                   dividevalue=0;
                elseif max(ppp2)<0.01
                    dividevalue=1;
                else
                   for dd=2:D
                       indicatorppp(dd)=(ppp2(dd)-ppp2(dd-1))/ppp2(dd-1);
                   end
                       indicatorppp(1)=min(indicatorppp)-0.001;
                      [value1,value2]=sort(indicatorppp,'descend');
                   for dd=1:D
                      if ppp2(value2(dd))<10
                      if ppp2(value2(dd))>0.1
                      dividevalue=ppp2(value2(dd))-0.001;
                      break;
                      end
                     elseif ppp2(value2(dd)-1)<0.01
                     dividevalue=ppp2(value2(dd))-0.001;
                    break;
                      end
                      if dd==D
                           dividevalue=ppp2(value2(dd))-0.001;
                      end
                    end     
                end
             else
                  for kkk2=1:total/5
                      spos(kkk2,:)=xvalbest(:);
                 end
             for d=1:D
                for k=1:total/5
                      spos(k,d)=xvalbest(d)-0.1*total+1*k;
                 end
                 e=feval(fhd,spos',vara(:));
                 FEs=FEs+total/5;
                 bbpbbp(d)=abs(max(e)/arfitnessbest);
                 for k=1:total/5
                     spos(k,d)=xvalbest(d);
                 end
             end
             
              if max(bbpbbp)<3.1
                 for d=1:D
                     bbpbb(d)=1;
                 end
              else
                  [aaa1,aaa2]=sort(bbpbbp);
                  for d=1:D-1
                      diaaa1(d)=aaa1(d+1)/aaa1(d);
                  end
                      [aab1,aab2]=sort(diaaa1,'descend');
                  if aaa1(D/2)<=2
                    for d=1:D-1
                       if aaa1(aab2(d))<1.8
                          division=aaa1(aab2(d))+0.01;
                          break;
                       end
                    end
                    for d=1:D
                       if bbpbbp(d)<=division
                          bbpbb(d)=1;
                       else
                           bbpbb(d)=0;
                       end
                    end
                  else
                      for d=1:D-1
                        if aaa1(aab2(d))<4
                          division=aaa1(aab2(d))+0.01;
                          break;
                        else division=0;
                       end
                      end
                  for d=1:D
                    if bbpbbp(d)<=division
                       bbpbb(d)=1;
                    else
                        bbpbb(d)=0;
                    end
                  end
                  end
              end
             end
              kk=1;
              cc2=0;
             
               VRmin=repmat(Xmin,total,D);
               VRmax=repmat(Xmax,total,D);
               pos=VRmin+(VRmax-VRmin).*rand(total,D);
               FEs=FEs+total;
               
            while FEs<Max_FEs-total
            e1=feval(fhd,pos',vara(:));
            FEs=FEs+total;
            [a1,a2]=sort(e1);
         % a1(1)
          xmin(kk)=a1(1);
            if kk==1
            [arfitnessbest,seq]=min(arfitnessbest);
            end
            if a1(1)<arfitnessbest
                xy(kk)=a1(1);
                xyvector=pos(a2(1),:);
            else
                xy(kk)=arfitnessbest;
                 xyvector=xvalbest(:,seq(1));
            end
            
             if FEs-0.01*Max_FEs<total
        processvalue(1)=min(xy);
    elseif FEs-0.02*Max_FEs<total
         processvalue(2)=min(xy);
    elseif FEs-0.03*Max_FEs<total
         processvalue(3)=a1(1);
         elseif FEs-0.05*Max_FEs<total
         processvalue(4)=min(xy);
         elseif FEs-0.1*Max_FEs<total
         processvalue(5)=min(xy);
         elseif FEs-0.2*Max_FEs<total
         processvalue(6)=min(xy);
         elseif FEs-0.3*Max_FEs<total
         processvalue(7)=min(xy);
         elseif FEs-0.4*Max_FEs<total
         processvalue(8)=min(xy);
    elseif FEs-0.5*Max_FEs<total
         processvalue(9)=min(xy);
         elseif FEs-0.6*Max_FEs<total
         processvalue(10)=min(xy);
         elseif FEs-0.7*Max_FEs<total
         processvalue(11)=min(xy);
         elseif FEs-0.8*Max_FEs<total
         processvalue(12)=min(xy);
         elseif FEs-0.9*Max_FEs<total
         processvalue(13)=min(xy);
%              elseif FEs-0.99*Max_FEs<1*total
%          processvalue(14)=min(xy);
    end
            
            newpos=pos(a2(1:mu),:);
            meanval=(newpos(:,1:D)'*weights)';
             stdval=1*std(newpos);
%              if max(stdval)<0.0001
%                  break;
%              end
            if kk==1
                if D>=50
                   for jjj=1:D
                     if bbpbb(jjj)==0
                     stdval(jjj)=0.001;
                     meanval(jjj)=xvalbest(jjj);
                     end
                   end
                else
                    for jjj=1:D
                      if ppp1(jjj)<dividevalue
                          stdval(jjj)=0.001;
                     meanval(jjj)=xvalbest(jjj,seq(1));
                      end
                    end
                end
            end
            kk=kk+1;
            if kk>30 
               if mod(kk,20)==0
                   [aaa,bbb]=min(xmin);
               if bbb<kk-20
                   cc2=1;
               else
                   cc2=0;
                end
               end
            end
            for k=1:total
            if cc2==1      
                   a=0.96*randn(1,D);
            else 
                   a=1*randn(1,D);
            end
             pos(k,:)=meanval+stdval.*a;
            end           
             for k=1:total
                 for j=1:D
                    if pos(k,j)>100
                        pos(k,j)=meanval(j)+stdval(j).*randn;
                   elseif pos(k,j)<-100
                        pos(k,j)=meanval(j)+stdval(j).*randn;
                    end
                 end
              end
%              pos=pos(1:total,:);
%              pos(1,:)=xyvector(:);
            xiaomei(jj,kk+xiuxiuxiu-2)=min(xy);
            end
            
            
            cri_fes(i,jj) =FEs<Max_FEs;
               processvalue(14)=min(xy);
         finalvalue(jj)=min(xy);  
         process(jj,:)=processvalue;
         clear xy;
         clear spos;
         meanval;
    end 
    process=process'-100*func_num;
    finalvalue=finalvalue-100*func_num;
    for pr=1:14
        for pc=1:runs
            if process(pr,pc)<10^(-8)
                process(pr,pc)=0;
            end
        end
    end 
    wwww=0;

     clear process;
     clear xiaomei;

       disp(mean(finalvalue))
end
clear arz;
clear arxx;
end
toc
end

    

