function [ S, varargout ] = Bilinear_Newmark_v031218( T, z, dy, alpha, ag, Dtg, Dt )
% Author:  Nicolas Luco
% Revised by: Jawad Fayaz
% Reference:  "Dynamics of Structures" (1995) by A.K. Chopra 
%
% INPUT:
% ======
% T     = periods                           ( 1 x num_oscillators, or scalar )
% z     = damping ratios                    ( " )
% dy    = yield displacements               ( " )
% alpha = strain hardening ratios           ( " )
% ag    = ground acceleration time history  ( length(ag) x 1 )
% Dtg   = time step of ag                   ( scalar )
% Dt    = analyis time step                 ( scalar )
%
% OUTPUT:
% =======
% S.d  = relative displacement spectrum               ( 1 x num_oscillators ) 
%  .v  = relative velocity spectrum                   ( " )
%  .a  = acceleration spectrum                        ( " )
% varargout{1} = E, where ...
% E.K  = kinetic energy (at end)                      ( 1 x num_oscillators )
%  .D  = energy dissipated by damping (at end)        ( " )
%  .S  = strain energy (at end)                       ( " )
%  .Y  = energy dissipated by yielding (at end)       ( " )
%  .I  = total input energy (at end)                  ( " )
% varargout{2} = H, where ...
% H.d  = relative displacement response time history  ( length(ag) x num_oscillators ) 
%  .v  = relative velocity response time history      ( " )
%  .a  = acceleration response time history           ( " )
%  .fs = force response time history                  ( " )
%


% Newmark parameters definition
% -----------------------------
GAMMA = 1/2;
BETA = 1/6;   % linear acceleration (stable if Dt/T<=0.551)
%BETA = 1/4;   % average acceleration (unconditionally stable)
    
% Number of oscillators
% ---------------------
num_oscillators = max( [ length(T), length(z), length(dy), length(alpha) ] );

% m*a + c*v + fs(k,fy,kalpha) = p
% -------------------------------
m = 1;
w = 2*pi ./ T;
c = z .* (2*m*w);
k = (w.^2 * m) .* ones(1,num_oscillators);
fy = k .* dy;
kalpha = k .* alpha;

% Interpolate p=-ag*m (linearly)
% ------------------------------
tg = [ 0:(length(ag)-1) ]' * Dtg;
t = [ 0:Dt:tg(end) ]';
p = interp1( tg, -ag*m, t );

% Memory allocation & initial conditions
% --------------------------------------
lp = length( p );
d = zeros( lp, num_oscillators );
v = zeros( lp, num_oscillators );
a = zeros( lp, num_oscillators );
fs = zeros( lp, num_oscillators );

% Initial calculations
% --------------------
a(1,:) = ( p(1) - c.*v(1,:) - fs(1,:) ) / m;
A = 1/(BETA*Dt)*m + GAMMA/BETA*c;
B = 1/(2*BETA)*m + Dt*(GAMMA/(2*BETA)-1)*c;

% Time stepping
% -------------
for i = 1:(lp-1)

    DPi = p(i+1)-p(i) + A.*v(i,:) + B.*a(i,:);
    ki = k;
    jj = find( ( DPi>0 & fs(i,:)>=+fy+kalpha.*(d(i,:)-dy) ) | ...
               ( DPi<0 & fs(i,:)<=-fy+kalpha.*(d(i,:)+dy) ) );
    ki(jj) = kalpha(jj);
    Ki = ki + A/Dt;  % = ki + GAMMA/(BETA*Dt)*c + 1/(BETA*Dt^2)*m;
    
    Ddi = DPi ./ Ki;
    fs(i+1,:) = fs(i,:) + ki.*Ddi;
    d(i+1,:) = d(i,:) + Ddi;
    fsmax =  fy + kalpha.*(d(i+1,:)-dy);
    fsmin = -fy + kalpha.*(d(i+1,:)+dy);
    jjabove = find( fs(i+1,:) > fsmax );
    jjbelow = find( fs(i+1,:) < fsmin );
    if ~isempty( [ jjabove jjbelow ] )
        fs(i+1,jjabove) = fsmax(jjabove);
        fs(i+1,jjbelow) = fsmin(jjbelow);
        Df = fs(i+1,:) - fs(i,:) + (Ki-ki).*Ddi;
        DR = DPi - Df;
        Ddi = Ddi + DR ./ Ki;
        d(i+1,:) = d(i,:) + Ddi;
    end

    Dvi = GAMMA/(BETA*Dt)*Ddi - GAMMA/BETA*v(i,:) + Dt*(1-GAMMA/(2*BETA))*a(i,:);
    v(i+1,:) = v(i,:) + Dvi;
    a(i+1,:) = ( p(i+1) - c.*v(i+1,:) - fs(i+1,:) ) / m;
%   Dai = 1/(BETA*Dt^2)*Ddi - 1/(BETA*Dt)*v(i,:) - 1/(2*BETA)*a(i,:);  % alternative
%   a(i+1,:) = a(i,:) + Dai;  % alternative
   
end

% Spectral values
% ---------------
S.d = max( abs(d) );
S.v = max( abs(v) );
S.a = max( abs( a + (-p/m)*ones(1,num_oscillators) ) );  % Note: ag itself was not interpolated

% Energies at end
% ---------------
if nargout > 1
    E.K = m * v(end,:).^2 / 2;
    E.D = c .* trapz( v.^2 )*Dt;
    E.S = fs(end,:).^2 ./ (2*k);
    E.Y = trapz( fs.*v )*Dt - E.S;
    E.I = trapz( (p*ones(1,num_oscillators)).*v )*Dt;
    tmp = (E.K + E.D + E.S + E.Y) ./ E.I;
%     disp( [ 'Energy Balance:  ' num2str(min(tmp)) ...
%             ' < (E.K + E.D + E.S + E.Y) ./ E.I < ' num2str(max(tmp)) ] )
    varargout{1} = E;
end
    
% Histories at tg
% ---------------
if nargout > 2
    H.d = interp1( t, d, tg );
    H.v = interp1( t, v, tg );
    H.a = interp1( t, a, tg );
    H.fs = interp1( t, fs, tg );
    varargout{2} = H;
end
