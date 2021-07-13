global a lam_min lam_max f_Al a1 a2 a3 A
a = 4.050; % angstrom
lam_min = 0.6; % angstrom
lam_max = 1.1; % angstrom
f_Al = 13; % -
a1 = [1, -1, 1];
a2 = [1, 1, 0];
a3 = [-1, 1, 2];
A = [a1',a2',a3'];
 
L = 3; % camera length, cm
n = 8; % initial guess on the max of Miller indices 
KK = [n, n, n];
while (n - max(abs(max(KK)))) <= 1e-6
    n = n + 1;
    clear DATA t
    DATA = zeros(1,12); % data to be exported
    t = 0; % number of diffraction conditions
    for h = -n:1:n
        for k = -n:1:n
            for l = -n:1:n
                K = [h,k,l]; % K vector
                F = f_Al*(1+exp(1i*pi*(h+k))+exp(1i*pi*(k+l))+exp(1i*pi*(h+l)));
                FF = abs(F);
                if (h == 0) && (k == 0) && (l == 0)
                    t = t + 1;
                    KK(t,1:3) = K;
                    theta = 0.5*pi;
                    d = Inf;
                    deg = theta*180/pi; % degree
                    lam = 2*d*sin(theta);
                    rho = 0;
                    phi = 0;
                    m = 1;
                    [x, y] = pol2cart(phi, rho);
                    DATA(t,1:3) = K; 
                    DATA(t,4:12) = ...
                        [m, FF, d, deg, lam, rho, phi*180/pi, x, y];
                elseif FF >= 1e-6
                    S = A\K';
                    if S(1) > 0
                        theta = 0.5*2*(acos(sum((-a1).*K)/norm(-a1)/norm(K))-0.5*pi);
                        if 2*theta > 0.5*pi
                            d = a/sqrt(h^2+k^2+l^2);
                            LAM = 2*d*sin(theta);
                            m = 1;
                            lam = LAM/m;
                            while (lam <= lam_max) && (lam >= lam_min)
                                rho = L*tan(pi-(2*theta));
                                K = [h, k, l];
                                K0 = (-a1/norm(-a1))*1/lam;
                                K1 = K0 + K;
                                P = K1 - ((sum(K1.*(a1))/norm(a1))*(a1)/norm(a1));
                                P2 = P(1, 1:2);
                                A2 = A(1:2, 2:3);
                                S2 = A2\P2';
                                ph = acos(sum(P.*(-a2))/norm(P)/norm(-a2));
                                if S2(2) < 0
                                    phi = 2*pi-ph;
                                else
                                    phi = ph;
                                end
                                [x, y] = pol2cart(phi, rho);
                                if (abs(x) <= 5) && (abs(y) <= 5)
                                    t = t + 1;
                                    KK(t,1:3) = K;
                                    deg = theta*180/pi; % degree
                                    DATA(t,1:3) = K;
                                    DATA(t,4:12) = ...
                                        [m, FF, d, deg, lam, rho, phi*180/pi, x, y];
                                end
                                m = m + 1;
                                lam = LAM/m;
                            end
                        end
                    end
                end
            end
        end
    end
end
n = n - 1;
X = DATA(:,11);
Y = DATA(:,12);
scatter(X,Y,25,'ko','filled')
xlabel('position, x (cm)'), ylabel('position, y (cm)')