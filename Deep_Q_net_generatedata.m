function [ output ] = Deep_Q_net_generatedata( net_w,data )
%Deep_Q_net_generatedata 是深度Q网络的生成数据函数，输入data,当前网络，输出Q函数的目标

        %生成新的label
             hidnet=data*net_w.layer{1}+net_w.bias{1};
             hid=(1+exp(-hidnet)).^(-1)-0.5;
             net=hid*net_w.layer{2}+net_w.bias{2};
             output=((1+exp(-net)).^(-1)-0.5)*5;
         

end