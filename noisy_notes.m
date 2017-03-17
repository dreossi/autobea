clear all
close all


fs = 8000;
T = 3; 

t = 0:(1/fs):T;
t = t(1:end-19);
t = [t t+t(end)];

t = t(1:end-10);

f = 110;
a = 0.5;
c = 1;
eps = 0.001;

delta = t(2)-t(1);

traj.x = [0 345];
traj.t = t;

q = 0;

global ts;
ts = 0;

global xp;
xp = a;


figure(1)
hold on

alpha = 0;
P = 0;
epsP = 0.001;

for i=2:size(t,2)
    
    traj.x(i,1) = simModel(f,a,c,traj.t(i),delta,P);
    traj.x(i-1,2) = (traj.x(i,1) - traj.x(i-1,1))/delta;
    
    qa = [traj.x(i-1,2), traj.x(i-1,1), traj.t(i-1)];
    q = automaton(q,qa);
    
    P = sampleState(q,P,epsP,alpha);
    
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
    plot([traj.t(i-1) traj.t(i)],[traj.x(i-1) traj.x(i)],color);
    %drawnow
    
    %plotTrace(traj,q);   
        
end
        
    
function [y] = simModel(f, a, c, t, delta, eps)
%WAVEGEN Generate a dumped sound wave

    if(t <= 2.9976)
        y = a*sin(2*pi*f*t)*c*(exp(-t)+eps);
    else
        y = simModel(f,a,c,5.9940-t,delta,eps);
    end
    
end

function P = alterParams(P,eps)

end

function qp = automaton(q, a)

    global c
    global xp
    global p;
    
    p = 0.0091;
    
    EPSxdot = 15;
    EPSx = 0.15;
    EPSc = 0.0001;    
    delta = 1;
    
    T = 0.05;
    
    xdot = a(1);
    x = a(2);
    t = a(3);
    ticc(t);  
    
    switch q
        
        case 0
            
            if (c > T) && ixdot(xdot,EPSxdot)
                wxp(x);
                resetc(t);                
                qp = 1;
                plot(t,xp,'*k'); 
                return                
            end
            ticc(t);
            qp = 0;
            return     
        
        case 1
            
            if ixdot(xdot,EPSxdot)
                if(( xp > 0) && (x < xp)) || (( xp < 0) && (x > xp))
                    wxp(x);
                    resetc(t);                
                    qp = 2;
                    plot(t,xp,'*g');
                    return
                else
                    wxp(x);
                    resetc(t);                
                    qp = 3;
                    plot(t,xp,'*r');
                    return
                end
            end
            ticc(t);
            qp = 1;
            return
            
            
            
        case 2            
                      
            if ixdot(xdot,EPSxdot) && (c > 0.009)
                if (( xp > 0) && (x < xp)) || (( xp < 0) && (x > xp))
                    wxp(x);
                    resetc(t);                
                    qp = 2;   
                    plot(t,x,'*b');                
                    return
                end
                
                resetc(t);
                qp = 0;
                return                         
            end         
                  
            ticc(t);
            qp = 2;
            return 
            
        case 3            
                      
            if ixdot(xdot,EPSxdot) && (c > 0.009)
                if (( xp > 0) && (x > xp)) || (( xp < 0) && (x < xp))
                    wxp(x);
                    resetc(t);                
                    qp = 3;   
                    plot(t,x,'*b');                
                    return
                end
                resetc(t);
                qp = 0;
                return                         
            end         
                  
            ticc(t);
            qp = 3;
            return  
    end
end


% Aux guards
function b = ix(x, EPS)
    global xp
    %b = abs(x-xp) < EPS;
    if(xp > 0)
        b = x < xp;
    else
        b = x > xp;
    end
        
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
    %axis([0 100 0.5 3.5]);
    %drawnow
end

function P = sampleState(q, P, eps, alpha)
    
    adj = [0.5+alpha 0.5-alpha 0 0 ; 0 0 0.5 0.5 ; 0.5+alpha 0 0.5-alpha 0 ; 0.5+alpha 0 0 0.5-alpha];
    alt = [1 0 0 0 ; 0 0 0 0 ; 1 0 0 0 ; 1 0 0 0];
    
    
    probs = adj(q+1,:);
    toss = rand;
    acc = 0;
    
    P = 0;
    
    for i=1:length(probs)
        if probs(i) > 0
            if toss <= (probs(i) + acc)
                if alt(q+1,i)
                    alterParams(P,eps);
                    P = rand*eps - (eps/2);
                end
                return
            else
                acc = acc + probs(i);                
            end
        end
    end
    

end


