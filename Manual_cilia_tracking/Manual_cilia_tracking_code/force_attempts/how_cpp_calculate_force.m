clear
close all

% define angles
phi = linspace(0,2*pi,1e4 + 1);
phi(end) = [];

% ellipse's half axes
Rx = 1;
Ry = 0.5;

% angle of points on ellipse
phi0 = pi/6;

% distance of the trap from the point on the ellpse. This influences the
% magnitude of the force applied, and is in general a function of the angle
epsilon = 0.5; 

% bead's intended trajectory
Bx = Rx.*cos(phi0);
By = Ry.*sin(phi0);

% trap position with small angles approximation
phi1 = epsilon .* sqrt( Bx^2 + By^2 );
Tx = Rx.*cos(phi1);
Ty = Ry.*sin(phi1);




%% How BeadsTraps.cpp works:

% bead's position
x = Rx.*cos(phi(1000)) + 0.2*rand(1)-0.1;
y = Ry.*sin(phi(1000)) + 0.2*rand(1)-0.1;

% bead's position on the ellipse
theta0 = atan2(y/Ry, x/Rx); % Phase of the bead on the cycle (this is not an actual angle, but theta for the ellipse parametrization)
x0 = Rx*cos(theta0);
y0 = Ry*sin(theta0);

% distance of (x0,y0) from centre of ellipse
norm = sqrt(Rx*Rx*sin(theta0)*sin(theta0) + Ry*Ry*cos(theta0)*cos(theta0)); %"radius"

% tangent versor
tx = -Rx*sin(theta0)/norm; %from derivative of trajectory
ty = Ry*cos(theta0)/norm;

% trap position, epsilon away 
xtrap = x0 + epsilon * tx;
ytrap = y0 + epsilon * ty;

% plot
figure;
hold on
box on;
plot(Rx.*cos(phi), Ry.*sin(phi),'.');   % trajectory
% plot(Bx, By, '.r');
% plot(Tx, Ty, '.g');
axis image
plot(x,y,'cx');                         % bead's position
plot([0,x],[0,y],'c-.')                 % vector to centre
plot(x0,y0,'c+');                       % corresponding point on ellipse
plot(xtrap, ytrap, 'r*');               % trap's position
quiver(x0,y0,tx,ty,'g', 'AutoScale','off');     % tangent versor



%% How BeadsSolving.cpp works:

% x,y still bead's position, and xtrap, ytrap trap position
% same norm definition
% same tx,ty definition

% normal versor (also works with altternative definition ie swapping signs)
nx = -ty;
ny = tx;

% Vector bead-trap in t,n coordinates
xe = tx*(x-xtrap) + ty*(y-ytrap);
ye = nx*(x-xtrap) + ny*(y-ytrap);

% stiffness
kt = 0.5;
kn = 0.5;
alpha = 2;
beta = 2;

% Force in the (t, n) coordinates
% the alpha and beta in front of everything come from the potential being defined as V = k*x^alpha, so then force is F = -alpha*k*x^(alpha-1)
ForceT = - beta  * kt * xe * abs(xe)^(beta -2);
ForceN = - alpha * kn * ye * abs(ye)^(alpha-2);

% And force in the (x, y) coordinates
% signs seem to check out despite n being defined the opposite way of what I would have done
ForceX = ForceT*tx + ForceN*nx;
ForceY = ForceT*ty + ForceN*ny;

% plot
quiver(x0,y0,nx,ny,'r', 'AutoScale','off');     % tangent versor
quiver(x,y,ForceX, ForceY, 'b', 'AutoScale','off'); 
        
        
        