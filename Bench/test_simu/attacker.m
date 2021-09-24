function [dU,dY,BIU_out,BIY_out,Replay_data_out] = attacker(U,dU_1,Y,dY_1,Time,BIU_in,BIY_in,Replay_data_in,choix_scenario,  Data,Ts)
%Choix du scénario selon la case coché dans le subsystem
BIU_out=BIU_in;
BIY_out=BIY_in;
Replay_data_out=Replay_data_in;
switch choix_scenario
    case 2 %Deny of service
        
        if (Time >= Data.DoS.Tstart) && (Time <= Data.DoS.Tend)
            if rand()>Data.DoS.pu
                dU = dU_1;
            else
                dU = U;
            end
            if rand()>Data.DoS.py
                dY = dY_1;
            else
                dY = Y;
            end
        else
            dU = U;
            dY = Y;
        end
        
    case 3 % Bias injection
        
        if (Time >= Data.BI.Tstart) && (Time <= Data.BI.Tend)
            BIU_out = BIU_in*Data.BI.transientU + (1-Data.BI.transientU)*Data.BI.FinalBiasU;
            dU = U + BIU_out;
            BIY_out = BIY_in * Data.BI.transientY + (1-Data.BI.transientY)*Data.BI.FinalBiasY;
            dY = Y + BIY_out;
        else
            dU = U;
            dY = Y;
        end
        
    case 4 % Upper saturation
        
        if (Time >= Data.US.Tstart) && (Time <= Data.US.Tend)
            dU = Data.US.U;
            dY = Y;
        else
            dU = U;
            dY = Y;
        end
        
    case 5 % Replay attack
        if (Time >= Data.RA.Tstart_save) && (Time <= Data.RA.Tend_save) %Phase d'enregistrement
            Replay_data_out(:,(Time - Data.RA.Tstart_save)/Ts+1) = Y;
            dY = Y;
            dU = U;
        elseif (Time >= Data.RA.Tstart_atk) && (Time <= Data.RA.Tend_atk)%phase d'attack
            dU = Data.RA.U;
            dY = Replay_data_in(:,(Time - Data.RA.Tstart_atk)/Ts+1);
        else 
            dU = U;
            dY = Y;
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
        dY = Y;
        
end
    
end
