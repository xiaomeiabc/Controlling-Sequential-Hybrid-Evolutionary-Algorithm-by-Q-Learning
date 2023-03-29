function [criterion] = USE_DQN(kk,xiaomei,idea)
    xiaomei = log(xiaomei);
    if mod(kk,10)==0 && kk>=30
        state1 = (xiaomei(kk+1-20)-xiaomei(kk+1))/xiaomei(kk+1-20);
        state2 = (xiaomei(1)-xiaomei(kk+1))/xiaomei(1);
        

        meta=[1,10,15,18,29,2,6,7,8,9,10];
        qflag1=0;
        qflag2=0;
    for k=1:11
        if k>5

            state = [state1;state2];
            state=state*10;
            load(['network\cec13\D10_new_gai_',num2str(meta(k)),'_1'])
            load(['network\cec13\D10_new_gai_',num2str(meta(k)),'_2'])

            Q1=Deep_Q_net_generatedata( net_w1,state(:,1)' );
            Q2=Deep_Q_net_generatedata( net_w2,state(:,1)' );
            load(['data\cec13\D10_new_gai_',num2str(meta(k))])
            dis=sqrt((state(1,1)-state11(1,3))^2+(state(2,1)-state11(2,3))^2);

        for k1=4:19

            dis1=sqrt((state(1,1)-state11(1,k1))^2+(state(2,1)-state11(2,k1))^2);
            if dis1<dis
                dis=dis1;
            end
        end
        else

            state = [state1;state2];
            state=state*10;
            load(['network\cec17\D10_gai_',num2str(meta(k)),'_1'])
            load(['network\cec17\D10_gai_',num2str(meta(k)),'_2'])

            Q1=Deep_Q_net_generatedata( net_w1,state(:,1)' );
            Q2=Deep_Q_net_generatedata( net_w2,state(:,1)' );
            load(['data\cec17\D10_gai_',num2str(meta(k))])
            dis=sqrt((state(1,1)-state11(1,3))^2+(state(2,1)-state11(2,3))^2);

        for k1=4:19

            dis1=sqrt((state(1,1)-state11(1,k1))^2+(state(2,1)-state11(2,k1))^2);
            if dis1<dis
                dis=dis1;
            end
        end
        end
        
        
        if dis>idea
            Q1=0;
            Q2=0;
        end
        
        if Q1>Q2
            qflag1=qflag1+1;
        elseif Q1<Q2
            qflag2=qflag2+1;
        else
            qflag1=qflag1+1;
            qflag2=qflag2+1;
        end
    end
        
        if qflag1< qflag2
            criterion = 1;
        elseif qflag1> qflag2
            criterion = 0;
        else
            if rand(1,1)>=0.5
                criterion = 1;
            else
                criterion = 0;
            end
        end
        
        if kk>=200
           criterion = 1;
        end
    else
        criterion = 0;%%%continue
    end

end

