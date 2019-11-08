function [Ypred1,Y,t,x,expvar] = test_ceof(iopt);

% a simple test using traveling wave solutions as data
%
%  y = cos(k_1 x - omega_1 t) + 0.7 cos(k_2 x + omega_2 t)
%
%
%  set the range of x and t
x = (0 : pi/32 : 4 * pi);
t = (0 : 0.02 : 4);

% pick wavenumber and frequency
k1 = 1;
k2 = 4;
omega1 = 3;
omega2 = 6;

% construct data matrix Y(time, space)
[X,T]=meshgrid(x,t);
Y = cos(k1 * X - omega1 * T) +  0.7 * cos(k2*X + omega2*T) + 0.1*randn(size(X));

% complex eof
[e,pc,expvar]=calCeof(Y,3,4);

% reconstruct first and second mode
Ypred1 = pc(1,:).' * conj(e(1,:));
Ypred2 = pc(2,:).' * conj(e(2,:));
expvar

% plot
if (iopt ==1)
subplot(3,1,1)
contourf(Y);
xlabel('x')
ylabel('time')
colorbar
subplot(3,1,2)
contourf(real(Ypred1));
xlabel('x')
ylabel('time')
colorbar
subplot(3,1,3)
contourf(real(Ypred2));
xlabel('x')
ylabel('time')
colorbar
end

% spatial phase
if (iopt==2)
    kx1 = atan2(imag(Ypred1(50,:)),real(Ypred1(50,:)));
    kx2 = atan2(imag(Ypred2(50,:)),real(Ypred2(50,:)));
    plot(x,kx1,x,kx2,'--')
    xlabel 'x'
    ylabel 'phase kx'
 
end

if (iopt ==3)
    wt1 = atan2(imag(Ypred1(:,32)),real(Ypred1(:,32)));
    wt2 = atan2(imag(Ypred2(:,32)),real(Ypred2(:,32)));
    plot(t,wt1,t,wt2,'--')
    xlabel 't'
    ylabel 'phase wt'
end






