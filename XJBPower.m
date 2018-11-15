function [ power_xjb,v_xjb,q1out_xjb,qspill_xjb ] =Power(v1_xjb,v2_xjb,qin_xjb,qout_xjb)
%% 输入：v1_xjb初始库容，v2_xjb末库容，qin_xjb，入库，qout_xjb出库
%  输出：powerout_xjb出力，v_xjb末库容；q1out_xjb实际出库；qspill_xjb弃水量
%% 只有出力大于预想出力时才有弃水
%% KN：效率系数；KF:效率系数；Qelse:它用流量;hmax_xjb查XZ_xjb变量；
    KN=8.8; KF=0.933; QELSE=120; hmax_xjb=113.6; hmin_xjb=76; npre_xjb=2009;
  global zv_xjb ZQDown_xjb lose_xjb XZ_xjb 
  global dt vmin_xjb vsheji_xjb  qoutmin_xjb
    %**********************计算出力****************************
    vup_xjb=(v1_xjb+v2_xjb)/2;
    zup_xjb=insert(zv_xjb(:,2),zv_xjb(:,1),vup_xjb);
    zdown_xjb=insert(ZQDown_xjb(:,1),ZQDown_xjb(:,2),qout_xjb);
    hlost_xjb=insert(lose_xjb(:,1),lose_xjb(:,2),qout_xjb);
    hpower_xjb=zup_xjb-zdown_xjb-hlost_xjb;
    %当qout<Qelse，N1为0
    if(qout_xjb<QELSE)
        N1=0;
    else
        N1=KN*(qout_xjb-QELSE)*hpower_xjb/1.0E3;
    end
    %hpower小于50，大于150均不作发电
    if((hpower_xjb<hmin_xjb)||(hpower_xjb>hmax_xjb))
        N2=0;
    else
        N2=insert(XZ_xjb(:,1),XZ_xjb(:,2),hpower_xjb);
    end
    %确定发电量power，下一阶段库容v与末库容v2实质一样，
    %放水q1out实质与qout一样，弃水量qspill初步为0；
    
    %%白鹤滩XZ_xjb
    power_xjb=min(N1,KF*N2);
    v_xjb=v2_xjb;
    q1out_xjb=qout_xjb;
    qspill_xjb=0.0;
    
    %% 讨论发电小于保证出力，大于保证出力的情况
    isout=1;
    %% 发电小于npre,为使power=npre,此时qspill仍为0，计算v,q1out
    if(power_xjb<npre_xjb)
        while(isout)
            v_xjb=v1_xjb+(qin_xjb-q1out_xjb)*dt/1.0E8;
            vuptemp_xjb=(v1_xjb+v_xjb)/2;
            zuptemp_xjb=insert(zv_xjb(:,2),zv_xjb(:,1),vuptemp_xjb);
            zdowntemp_xjb=insert(ZQDown_xjb(:,1),ZQDown_xjb(:,2),q1out_xjb);
            hlosttemp_xjb=insert(lose_xjb(:,1),lose_xjb(:,2),q1out_xjb);
            hpowertemp_xjb=zuptemp_xjb-zdowntemp_xjb-hlosttemp_xjb;
            qqtemp_xjb=(npre_xjb*1.0E3)/(KN*hpowertemp_xjb)+QELSE;%保证出力计算的实际流量
            if(abs(q1out_xjb-qqtemp_xjb)<0.1)
                power_xjb=npre_xjb;
                v_xjb=v1_xjb+(qin_xjb-q1out_xjb)*dt/1.0E8;
                isout=0;
            else
                q1out_xjb=(q1out_xjb+qqtemp_xjb)/2;
            end
        end
     %% 计算发电大于保证出力时的弃水量,v=v2, q1out=qout,  power=min(N1,N2*KF);
    else
        vuptemp_xjb=(v1_xjb + v_xjb)/2;
        zuptemp_xjb=insert(zv_xjb(:,2),zv_xjb(:,1),vuptemp_xjb);
        zdowntemp_xjb=insert(ZQDown_xjb(:,1),ZQDown_xjb(:,2),q1out_xjb);
        hlosttemp_xjb=insert(lose_xjb(:,1),lose_xjb(:,2),q1out_xjb);
        hpowertemp_xjb=zuptemp_xjb-zdowntemp_xjb-hlosttemp_xjb;
        qyy_xjb=(power_xjb*1.0E3)/(KN*hpowertemp_xjb)+QELSE;            
        qspill_xjb = q1out_xjb - qyy_xjb ;
        if(abs(qspill_xjb-0)<0.01)
            qspill_xjb=0;
        end
    end
    
    %% 当出水比qoutmin小的情况,q1out=qmin,计算v,power,qspill=0
    if(q1out_xjb<qoutmin_xjb)
        q1out_xjb=qoutmin_xjb;
        v_xjb=v1_xjb+(qin_xjb-q1out_xjb)*dt/1.0E8;
        vuptemp_xjb=(v1_xjb + v_xjb)/2;
        zup_xjb=insert(zv_xjb(:,2),zv_xjb(:,1),vuptemp_xjb);
        zdown_xjb=insert(ZQDown_xjb(:,1),ZQDown_xjb(:,2),q1out_xjb);
        hlost_xjb=insert(lose_xjb(:,1),lose_xjb(:,2),q1out_xjb);
        hpower_xjb=zup_xjb - zdown_xjb -hlost_xjb;
        power_xjb=KN*(q1out_xjb - QELSE)*hpower_xjb/1.0E3;
    end
    
    %% 当库容小于vmin的情况,v=vmin,计算q1out,power,qspill=0,
    if(v_xjb<vmin_xjb)
        v_xjb=vmin_xjb;
        q1out_xjb=qin_xjb-(v_xjb-v1_xjb)*1.0E8/dt;
        vuptemp_xjb=(v1_xjb + v_xjb)/2.0;
        zup_xjb=insert(zv_xjb(:,2),zv_xjb(:,1),vuptemp_xjb);
        zdown_xjb=insert(ZQDown_xjb(:,1),ZQDown_xjb(:,2),q1out_xjb);
        hlost_xjb=insert(lose_xjb(:,1),lose_xjb(:,2),q1out_xjb);
        hpower_xjb=zup_xjb - zdown_xjb -hlost_xjb;
        power_xjb=KN*(q1out_xjb - QELSE)*hpower_xjb/1.0E3;
  %% 计算发电大于保证出力时的弃水量,v=v2, q1out=qout,  power=N1,qspill在之前已经计算过;;  
    elseif(v_xjb>vsheji_xjb)% 当库容大于vshejisx
        q1out_xjb=qin_xjb-(v_xjb-v1_xjb)*1.0E8/dt;
        vuptemp_xjb=(v1_xjb + v_xjb)/2.0;
        zup_xjb=insert(zv_xjb(:,2),zv_xjb(:,1),vuptemp_xjb);
        zdown_xjb=insert(ZQDown_xjb(:,1),ZQDown_xjb(:,2),q1out_xjb);
        hlost_xjb=insert(lose_xjb(:,1),lose_xjb(:,2),q1out_xjb);
        hpower_xjb=zup_xjb - zdown_xjb -hlost_xjb;
        power_xjb=KN*(q1out_xjb - QELSE)*hpower_xjb/1.0E3;
        %% qspill并不是没算，而是在之前已经计算过，实际也为0
    end
    
end