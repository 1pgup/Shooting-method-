clc, clear, close all
datetime('now')

%Planck's constant in J * sec
hbar = 1.0546e-34;
%Mass of electron in kg
m0 = 9.1094e-31;
%Conversion constant from joules to electronvolts
J2eV = 6.2415e18;

% Pit length in nm
L = 1;
%Potential in eV
U = 0;

%Cell array for energy comparison
C = cell(3, 4); 
C(1, 1) = {'n'};
C(2, 1) = {'E_shoot, eV'};
C(3, 1) = {'E_an, eV'};
C(1, 2) = {1};
C(1, 3) = {2};
C(1, 4) = {3};
%% 
% Shooting method

%Energies close to the analytically calculated energies of the first three energy levels
E_apr = [0.37, 1.504, 3.38]; %eV

%Dependence of the wave function on the coordinate
figure(1);
subplot(2, 3, 1);

%First energy level
[xp, psi_sh_1, E_sh_1] = shooting(E_apr(1), U, L, 0.001, 0.001);
plot(xp, psi_sh_1, 'LineWidth', 2);
hold on;

%Second energy level
[xp, psi_sh_2, E_sh_2] = shooting(E_apr(2), U, L, 0.001, 0.001);
plot(xp, psi_sh_2, 'LineWidth', 2);

%Third energy level
[xp, psi_sh_3, E_sh_3] = shooting(E_apr(3), U, L, 0.001, 0.001);
plot(xp, psi_sh_3, 'LineWidth', 2);

%Filling cell array with energies calculated by shooting method
C(2, 3) = num2cell(E_sh_2);
C(2, 2) = num2cell(E_sh_1);
C(2, 4) = num2cell(E_sh_3);

%Graphics customization
legend('\psi_{n = 1}', '\psi_{n = 2}', '\psi_{n = 3}', 'Location', 'southwest');
title('Shooting method');
grid on;
xlabel('x, nm');
ylabel('\psi, nm^{-0.5}');
xlim([0, L]);

%Dependence of the squared modulus of the wave function on the coordinate
subplot(2, 3, 4);

%First energy level
plot(xp, abs(psi_sh_1) .^2, 'LineWidth', 2);
hold on;

%Second energy level
plot(xp, abs(psi_sh_2) .^2, 'LineWidth', 2);

%Third energy level
plot(xp, abs(psi_sh_3) .^2, 'LineWidth', 2);

%Graphics customization
title('Shooting method');
grid on;
xlabel('x, nm');
ylabel('|\psi|^2, nm^{-1}');
xlim([0, L]);
%% 
% Analytical method

%Principal quantum numbers
n = [1, 2, 3];

%Energies of the first three energy levels
E_an = pi ^ 2 * n .^ 2 * hbar ^ 2 / (2 * m0 * (L * 1e-9) ^2) * J2eV;

%Filling cell array with analytically calculated energies
C(3, 2) = num2cell(E_an(1));
C(3, 3) = num2cell(E_an(2));
C(3, 4) = num2cell(E_an(3));

%Dependence of the wave function on the coordinate
subplot(2, 3, 2);

%First energy level
psi_an_1 = sqrt(2 / L) * sin(pi * n(1) * xp / L);
plot(xp, psi_an_1, 'LineWidth', 2);
hold on;

%Second energy level
psi_an_2 = sqrt(2 / L) * sin(pi * n(2) * xp / L);
plot(xp, psi_an_2, 'LineWidth', 2);

%Third energy level
psi_an_3 = sqrt(2 / L) * sin(pi * n(3) * xp / L);
plot(xp, psi_an_3, 'LineWidth', 2);

%Graphics customization
legend('\psi_{n = 1}', '\psi_{n = 2}', '\psi_{n = 3}', 'Location', 'southwest');
title('Analytical method');
grid on;
xlabel('x, nm');
ylabel('\psi, nm^{-0.5}');
xlim([0, L]);

%Dependence of the squared modulus of the wave function on the coordinate
subplot(2, 3, 5);

%First energy level
plot(xp, abs(psi_an_1) .^2, 'LineWidth', 2);
hold on;

%Second energy level
plot(xp, abs(psi_an_2) .^2, 'LineWidth', 2);

%Third energy level
plot(xp, abs(psi_an_3) .^2, 'LineWidth', 2);

%Graphics customization
title('Analytical method');
grid on;
xlabel('x, nm');
ylabel('|\psi|^2, nm^{-1}');
xlim([0, L]);
sgtitle('Electron in 1d-PW');
%% 
% Difference

%Abs difference between psi_sh and psi_an
subplot(2, 3, 3)

%psi_sh_1 and psi_an_1
plot(xp, abs(psi_an_1 - psi_sh_1), 'LineWidth', 2);
hold on;
%psi_sh_2 and psi_an_2
plot(xp, abs(psi_an_2 - psi_sh_2), 'LineWidth', 2);

%Graphics customization
title('Difference');
legend('|\psi_{an}^{n = 1} - \psi_{shoot}^{n = 1}|', ...
    '|\psi_{an}^{n = 2} - \psi_{shoot}^{n = 2}|', ...
    'Location', 'northwest');
grid on;
xlabel('x, nm');
ylabel('|\psi_{an} - \psi_{shoot}|, nm^{-0.5}');
xlim([0, L]);
sgtitle('Electron in 1d-PW');

subplot(2, 3, 6)
%psi_sh_3 and psi_an_3
plot(xp, abs(psi_an_3 - psi_sh_3), 'LineWidth', 2);
hold on;

%Graphics customization
title('Difference');
grid on;
xlabel('x, nm');
ylabel('|\psi_{an}^{n = 3} - \psi_{shoot}^{n = 3}|, nm^{-0.5}');
xlim([0, L]);
sgtitle('Electron in 1d-PW');

%Outputting filled cell array
C
%% 
% Функция находит $\psi$, удовлетворяющую второму граничному условию $\psi \left(L\right)=0$.

function [xp, psip, E] = shooting(E, U, L, dx, tolerance)
    %hbar * c [eV * nm]
    hbc = 1.0546e-34 * 6.2415e18 * 3e17;
    %mass of electron in mc^2
    m = 510998.95;
    %Coefficient in the Schrodinger equation
    k = 2 * m / hbc ^2;
  
    %Increasing Energy to Satisfy the Second Boundary Condition
    dE = 0.001; %eV

    %Vector for storing the coordinate
    xp = dx : dx : L;

    %First boundary condition
    psi = 0;
    %The last point of the wave function should be approximately equal to zero
    %We assume that it is not equal to 0
    psifinal = 1;

    %Loop until psifinal is close to zero
    while abs(psifinal) > tolerance
        %First boundary condition
        dpsi = 1;

        %Vector for storing value of the wave function
        psip = zeros(1, numel(xp));
        
        %Loop until we reach the end of the pit
        for i = 1 : L / dx
            %Schrodinger equation
            ddpsi = -k * psi * (E - U);

            %Representing the derivative as a limit
            dpsi = dpsi + ddpsi * dx;
            psi = psi + dpsi * dx;

            %Filling vector
            psip(i) = psi;
        end
        %Assigning the last value of the wave function to psifinal
        psifinal = psi;

        %Increasing the energy to the desired value only 
        % if the second boundary condition is not met
        if (abs(psifinal) > tolerance)
            E = E + dE;
        end
    end

    %Normalization
    area = sum(psip .^ 2)  * dx;

    psip = psip / sqrt(area);
end