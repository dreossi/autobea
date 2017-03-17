clear all
close all


addpath ./models/
system_model = 'laub_loomis';
load_system(system_model);

delta = 0.2;
traj.x = [ 1 ; 1; 1; 1; 1; 1; 0.1]; %[x7,x1,x2,x3,x4,x5,x6]
traj.y = [  0; traj.x(1,end) ];
traj.t = 0;
P = [0.9 2.5 1.5 0.6 0.8 1.0 1.3 0.3 0.8 0.7 4.9 23 4.5];


q = 0;

hold on
plot(traj.t(end),traj.y(2,end),'.b')
axis([0 100 0.5 3.5]);
drawnow

global ts;
ts = 0;

alpha = 0;
epsP = 0.02;

set(gca,'FontSize',18)

for i=1:500
       
    [traj.x(:,end+1),traj.y(:,end+1),traj.t(end+1)] = simModel(traj.x(:,end),traj.t(end),delta);
    q = automaton(q,[traj.y(:,end)' traj.t(end)]);
    
    sampleState(q,P,epsP,alpha);
    
    plotTrace(traj,q);   
        
end
% for i=1:125
%     
%     traj.y(:,end+1) = [0 ; traj.y(2,end)];
%     traj.t(end+1) = traj.t(end) + delta;   
%     
%     q = automaton(q,[traj.y(:,end)' traj.t(end)])
%     
%     plotTrace(traj,q);
%     
% end
% 
% for i=1:175
%      
%     [traj.x(:,end+1),traj.y(:,end+1),traj.t(end+1)] = simModel(traj.x(:,end),traj.t(end),delta);
%     
%     q = automaton(q,[traj.y(:,end)' traj.t(end)])
%     
%     plotTrace(traj,q);
%     
%         
% end
% 
% for i=1:75
%     
%     traj.y(:,end+1) = [0 ; traj.y(2,end)];
%     traj.t(end+1) = traj.t(end) + delta;   
%     
%     q = automaton(q,[traj.y(:,end)' traj.t(end)])
%     
%     plotTrace(traj,q);
%     
% end
        
    


function [x1,y1,t1] = simModel(x0, t0, delta)

    system_model = 'laub_loomis.slx';
    t1 = t0 + delta;
    
    out = sim(system_model,'StartTime',num2str(t0),'StopTime',num2str(t1),...
            'SolverType','variable-step','MaxStep','auto',...
            'InitInArrayFormatMsg', 'None','LoadInitialState','on','InitialState',mat2str(x0),...
            'SaveFinalState','on','SaveFormat','Structure','SaveFinalState','on',...
            'FinalStateName','xFinal','StateSaveName','xoutNew','TimeSaveName','tout',...
            'srcworkspace','current',...
            'AlgebraicLoopSolver','LineSearch');

    x1 = zeros(1,7);
    for i = 1:7
        x1(i) = out.get('xFinal').signals(i).values;
    end
    
    y1 = zeros(1,2);
    for i = 1:2
        y1(i) = out.get('yout').signals(i).values(end);
    end

end

function P = alterParams(P,eps)

    for i=1:length(P)
        P(i) = P(i) + 2*eps - eps;
    end

    set_param('laub_loomis/Constant2','Value',num2str(P(1)));
    set_param('laub_loomis/Constant3','Value',num2str(P(2)));
    set_param('laub_loomis/Constant4','Value',num2str(P(3)));
    set_param('laub_loomis/Constant5','Value',num2str(P(4)));
    set_param('laub_loomis/Constant6','Value',num2str(P(5)));
    set_param('laub_loomis/Constant7','Value',num2str(P(6)));
    set_param('laub_loomis/Constant8','Value',num2str(P(7)));
    set_param('laub_loomis/Constant9','Value',num2str(P(8)));
    set_param('laub_loomis/Constant10','Value',num2str(P(9)));
    set_param('laub_loomis/Constant11','Value',num2str(P(10)));
    set_param('laub_loomis/Constant12','Value',num2str(P(11)));
    set_param('laub_loomis/Constant13','Value',num2str(P(12)));
    set_param('laub_loomis/Constant14','Value',num2str(P(13)));

end

function qp = automaton(q, a)

    global c
    EPSxdot = 0.125;
    EPSx = 0.15;
    EPSc = 0.1;    
    delta = 1;
    
    T = 7.4;
    
    xdot = a(1);
    x = a(2);
    t = a(3);
    ticc(t);
    
    
    switch q
        case 0
            if (c > T) && ixdot(a(1),EPSxdot)
                wxp(x);
                resetc(t);
                qp = 1;
                return                
            end
            ticc(t);
            qp = 0;
            return
            
        case 1  
            
            if period(T,EPSc) && ix(x,EPSx) && ixdot(xdot,EPSxdot)
                wp();
                wxp(x);
                resetc(t);                
                qp = 2;
                
                plot(t,x,'*r')
                drawnow
                
                return
            end         
            if (c < T) && (xdot == 0)
                wp();
                wxp(x);
                resetc(t);                
                qp = 3;
                return;
            end               
            ticc(t);
            qp = 1;
            return
            
        case 2
            
            if ic(EPSc) && ix(x,EPSx) && ixdot(xdot,EPSxdot)                
                wxp(x);
                resetc(t);               
                qp = 2;   
                
                plot(t,x,'*k')
                return
            end
            if ic(EPSc) && (~(ix(x,EPSx)) || ~(ixdot(xdot,EPSxdot)))
                resetc(t);               
                qp = 0;
                return;
            end           
            ticc(t);
            qp = 2; 
            return; 
        case 3
            if ic(EPSc) && ix(x,EPSx)
               resetc(t);               
               qp = 3;
               return;
            end
            if ic(EPSc) && ~(ix(x,EPSx))
                resetc(t);               
               qp = 0;
               return;
            end
            ticc(t);
            qp = 3; 
            return;             
    end
end


% Aux guards
function b = ix(x, EPS)
    global xp
    b = abs(x-xp) < EPS;
end

function b = ixdot(xdot, EPS)
    b = abs(xdot) < EPS;
end

function b = ic(EPS)
    global c p;
    b = abs(c-p) < EPS;
end

function wp()
    global c p;
    p = c;
end

function wxp(x)
    global xp;
    xp = x;
end

function resetc(t)
    global c ts;
    c = 0;
    ts = t;    
end

function ticc(t)
    global c ts;
    c = t - ts; 
end

function b = period(T,EPS)
    global c
    b = abs(c - T) < EPS;
end

function plotTrace(traj,q)
    switch q
        case 0
            color = '-b';            
        case 1
            color = '-c';
        case 2
            color = '-g';
        case 3
            color = '-r';            
    end
    hold on
    plot([traj.t(end-1) traj.t(end)],[traj.y(2,end-1) traj.y(2,end)],color);
    axis([0 100 0.5 3.5]);
    drawnow
end

function P = sampleState(q, P, eps, alpha)
    
    adj = [0 1 0 0 ; 0 0 0.5-alpha 0.5+alpha ; 0.5+alpha 0 0.5-alpha 0 ; 0.5+alpha 0 0 0.5-alpha];
    alt = [0 0 0 0 ; 0 0 0 1 ; 1 0 0 0 ; 1 0 0 0];
    
    
    probs = adj(q+1,:);
    toss = rand;
    acc = 0;
    
    for i=1:length(probs)
        if probs(i) > 0
            if toss <= (probs(i) + acc)
                if alt(q+1,i)
                    alterParams(P,eps);
                end
                return
            else
                acc = acc + probs(i);                
            end
        end
    end
    

end


