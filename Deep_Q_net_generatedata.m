function [ output ] = Deep_Q_net_generatedata( net_w,data )
%Deep_Q_net_generatedata �����Q������������ݺ���������data,��ǰ���磬���Q������Ŀ��

        %�����µ�label
             hidnet=data*net_w.layer{1}+net_w.bias{1};
             hid=(1+exp(-hidnet)).^(-1)-0.5;
             net=hid*net_w.layer{2}+net_w.bias{2};
             output=((1+exp(-net)).^(-1)-0.5)*5;
         

end