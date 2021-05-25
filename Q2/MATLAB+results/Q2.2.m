% define system parameters:

c =     50 ;    % kg/s
k =     10000;     % N/m
m =       331 ;  % kg

Y=.05;
v=(65.*1.60934.*1000)./3600;
wb=v.*2.*pi;

% define the time interval and the number of integration steps
n    = 100000; 
t    = zeros(n,1);
tf   = 10;
t(1) = 0;
delt = (tf-t(1))/n;

% initializing variables
x    = zeros(n,1);
v    = zeros(n,1);

% assign initial conditions
x(1) = 0;
v(1) = 0;


% start the integration directly from the second step (with i=2)
for i=2:1:n
    % update the velocity and acceleration for x at step i 
    a=(-c./m).*v(i-1)+(-k./m).*x(i-1)+((c.*Y.*wb)./m).*cos(wb.*t(i-1))+((k.*Y)./m).*sin(wb.*t(i-1));
    
    % integrate x and v
    x(i) = x(i-1)+v(i-1).*delt;
    v(i)= v(i-1)+a.*delt;
    
    % integrate time
    t(i)=t(i-1)+delt;
    
end % end of each integration step  

% plot results 
figure(20)
clf; 
% plot vertical position
subplot(2,1,1) 
box on; grid on; hold on;
plot(t,x,'linewidth',2);
ylabel('Position [m]', 'fontsize', 14, 'fontname', 'times');
% plot vertical velocity
subplot(2,1,2) 
box on; grid on; hold on;
plot(t,v,'linewidth',2);
xlabel('t [s]', 'fontsize', 14, 'fontname', 'times');
ylabel('Velocity [m/s]', 'fontsize', 14, 'fontname', 'times'); 
