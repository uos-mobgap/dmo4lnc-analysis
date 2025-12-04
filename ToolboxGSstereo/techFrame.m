function newTraj = techFrame(Traj,mrks, mrkToRec, gap)
    N = length(Traj.(mrks{1}));
    
    % LRF
     for i = 1:N
        x(i,:) = (Traj.(mrks{2})(i,:) - Traj.(mrks{1})(i,:))/...
            norm(Traj.(mrks{2})(i,:) - Traj.(mrks{1})(i,:));     
        v(i,:) = (Traj.(mrks{3})(i,:) - Traj.(mrks{1})(i,:))/...
            norm(Traj.(mrks{3})(i,:) - Traj.(mrks{1})(i,:)); 
        z(i,:) = cross(x(i,:),v(i,:))/norm(cross(x(i,:),v(i,:)));                 

        y(i,:) = cross(z(i,:), x(i,:))/norm(cross(z(i,:), x(i,:)));            

        R(:,:,i) = [x(i,:)' y(i,:)' z(i,:)'];  % LRF Rotation Matix           
     end
     
     t = Traj.(mrks{1});   % Origin of the LRF
     
     newTrajL = Traj.(mrkToRec);
     newTraj = Traj.(mrkToRec);
     %% Rigid body recontruction
     GAPS = find(diff(gap)>1);
     Gaps_start =[gap(1); gap(GAPS+1)];
     Gaps_end = [gap(GAPS); gap(end)];
     for i = 1:size(Gaps_start,1)
         if Gaps_start(i)==1 % Gap at the beginning of the trial
             n = Gaps_end(i)+1;
             if n(end)>N
                 n = gap(end);
                 fprintf('Marker missing on the whole trial - RB reconstruction still not possible. \n')
             end
             
             for k = 1:n
                 Lmrk = (R(:,:,k)'*(Traj.(mrkToRec)(n,:) - t(k,:))')';   % from GRF to LRF
                 newTrajL(k,:) = Lmrk;   
             end
         elseif Gaps_end(i)==N                                         % Gap at the end of the trial
             n = Gaps_start(end)-1;
             for k = Gaps_start(i):Gaps_end(i)
                 Lmrk = (R(:,:,k)'*(Traj.(mrkToRec)(n,:) - t(k,:))')';   % from GRF to LRF
                 newTrajL(k,:) = Lmrk;   
             end
         else
             n = Gaps_end(i)+1;
             for k = Gaps_start(i):n
                 Lmrk = (R(:,:,k)'*(Traj.(mrkToRec)(n,:) - t(k,:))')';   % from GRF to LRF
                 newTrajL(k,:) = Lmrk;
             end
         end
         
         for k = Gaps_start(i):Gaps_end(i)
             newTraj(k,:) = (R(:,:,k)*newTrajL(k,:)'+ t(k,:)')';
         end
     end     
end