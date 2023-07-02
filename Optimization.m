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

%Principal quantum number
n = 5;

%Energy of the fifth energy level
E_an = pi ^ 2 * n ^ 2 * hbar ^ 2 / (2 * m0 * (L * 1e-9) ^2) * J2eV;

%Energy close to the analytically calculated energy
E_apr = 0.8 * E_an; %eV

%Dependence of the wave function on the coordinate
subplot(2, 1, 1);

%Fifth energy level
tic
[xp, psi_sh_1, E_sh] = shooting(E_apr, U, L, 0.001, 0.0001);
toc
plot(xp, psi_sh_1, 'LineWidth', 2);

%Graphics customization
grid on;
xlabel('x, nm');
ylabel('\psi, nm^{-0.5}');
xlim([0, L]);

%Dependence of the wave function on the coordinate
subplot(2, 1, 2);

%Fifth energy level
tic
%Using modified function
[xp, psi_sh_1, E_sh] = shooting_modified(0.8 * E_apr, U, L, 0.001, 0.0001);
toc
plot(xp, psi_sh_1, 'LineWidth', 2);

%Graphics customization
grid on;
xlabel('x, nm');
ylabel('\psi, nm^{-0.5}');
xlim([0, L]);
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

function [xp, psip, E] = shooting_modified(E, U, L, dx, tolerance)
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
        %temporary variable to hold the old value of psifinal
        tmp = psifinal;

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

        %If psifinal does not change significantly, then exit from the loop
        if (abs(psifinal - tmp)) < 1e-5
            break;
        end

        %Increasing the energy to the desired value only 
        % if the second boundary condition is not met
        if (abs(psifinal) / tolerance >= 5)
            E = E + 4 * dE;
        elseif (abs(psifinal) / tolerance >= 2)
            E = E + 2 * dE;
        elseif (abs(psifinal) / tolerance >= 1.2)
            E = E + dE / 2;
        end
    end

    %Normalization
    area = sum(psip .^ 2)  * dx;

    psip = psip / sqrt(area);
end