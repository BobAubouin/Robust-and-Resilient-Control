function [dU,BIU_out] = attacker_U(U,dU_1,Time,BIU_in,choix_scenario,Data,Ts)
%Choix du scénario selon la case coché dans le subsystem
BIU_out=BIU_in;

switch choix_scenario
    case 2 %Deny of service
        
        if (Time >= Data.DoS.Tstart) && (Time <= Data.DoS.Tend)
                dU(1:2) = 0*dU_1(1:2);
                %dU(3:5) = 0*dU_1(3:5);
        else
            dU = U;
        end
        
    case 3 % Bias injection
        
        if (Time >= Data.BI.Tstart) && (Time <= Data.BI.Tend)
            BIU_out = BIU_in*Data.BI.transientU + (1-Data.BI.transientU)*Data.BI.FinalBiasU;
            dU = U + BIU_out;
        else
            dU = U;
        end
        
    case 4 % Upper saturation
        
        if (Time >= Data.US.Tstart) && (Time <= Data.US.Tend)
            dU = Data.US.U;
        else
            dU = U;
        end
        
    case 5 % Replay attack
        if (Time >= Data.RA.Tstart_atk) && (Time <= Data.RA.Tend_atk)%phase d'attack
            dU = Data.RA.U;
        else 
            dU = U;
        end
        
%     case 6 % False disturbance prediction
%         dU = U;
%         if (Time >= Data.FdP.Tstart_save) && (Time <= Data.FdP.Tend_save) %Phase d'attaque
%            dY = data.FdP.Y(Data.FdP.i);
%            Data1.FdP.i = Data.FdP.i + 1;
%         else
%            dY = Y;
%         end
        
    otherwise %%nothing 
        
        dU = U;
        
end
    
end
