function [dY,BIY_out,Replay_data_out,pred_out] = attacker_Y(Y,pred,Time,BIY_in,Replay_data_in,choix_scenario,Data,Ts)
%Choix du scénario selon la case coché dans le subsystem
BIY_out=BIY_in;
Replay_data_out=Replay_data_in;
pred_out = pred;
switch choix_scenario
    case 2 %Deny of service
        
        if (Time >= Data.DoS.Tstart) && (Time <= Data.DoS.Tend)
            dY = Y*0;
        else
            dY = Y;
        end
        
    case 3 % Bias injection
        
        if (Time >= Data.BI.Tstart) && (Time <= Data.BI.Tend)
            BIY_out = BIY_in * Data.BI.transientY + (1-Data.BI.transientY)*Data.BI.FinalBiasY;
            dY = Y + BIY_out;
        else
            dY = Y;
        end
        
    case 4 % Upper saturation
        
        dY = Y;
        
    case 5 % Replay attack
        if (Time >= Data.RA.Tstart_save) && (Time <= Data.RA.Tend_save) %Phase d'enregistrement
            Replay_data_out(:,(Time - Data.RA.Tstart_save)/Ts+1) = Y;
            dY = Y;
        elseif (Time >= Data.RA.Tstart_atk) && (Time <= Data.RA.Tend_atk)%phase d'attack
            dY = Replay_data_in(:,(Time - Data.RA.Tstart_atk)/Ts+1);
        else 
            dY = Y;
        end
        
    case 6 % False disturbance prediction
        if (Time >= Data.FDP.Tstart_save) && (Time <= Data.FDP.Tend_save) %Phase d'attaque
            pred_out = pred;
            for i=1:size(pred_out,1)
                pred_out(i,:,:) = 0*pred_out(i,:,:) + Data.FDP.Bias(i);
            end
        end
        dY = Y;
        
    otherwise %%nothing 
        
        dY = Y;
        
end
    
end
