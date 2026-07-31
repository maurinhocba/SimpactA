%**************************************************************************
% Debug del point con regla de masing
%**************************************************************************
clear all;close all;clc
%**************************************************************************
% desplazamientos 

tini = 0;
tfin = 2;
proj = 'point5';
node = 2;
nEstrats = 1;
fr=1;

% ABRIR ARCHIVOS
correr=false;
for i=1:nEstrats
        str1=['(echo ' proj '.@' num2str(i) ' ' num2str(fr) ' 0 && echo d &&'];
        str2=[' echo ' num2str(node) ' ' num2str(0) ' u' num2str(i) '.dat && echo s) | curvas'];
        [status,cmdout]=dos([str1 str2]);
        
        if status==0
            [aux, msgU] = fopen(['u' num2str(i) '.dat'],'r');
            fidU(i) = aux;
            correr(i)=isempty(msgU) ;
        else
            correr(i)=0;
        end
end

correr=false;
for i=1:nEstrats
        str1=['(echo ' proj '.@' num2str(i) ' ' num2str(fr) ' 0 && echo l &&'];
        str2=[' echo ' num2str(node) ' ' num2str(0) ' f' num2str(i) '.dat && echo s) | curvas'];
        [status,cmdout]=dos([str1 str2]);
        
        if status==0
            [aux, msgU] = fopen(['f' num2str(i) '.dat'],'r');
            fidU(i) = aux;
            correr(i)=isempty(msgU) ;
        else
            correr(i)=0;
        end
end
fclose all;

uSimpact = [];
fSimpact = [];
for i=1:nEstrats
      uSimpact = [uSimpact;load(['u' num2str(i) '.dat'])];
      fSimpact = [fSimpact;load(['f' num2str(i) '.dat'])];
end

time   = uSimpact(:,1);

% rotm = eul2rotm(datD01(1,5:end))
% Trabaja con los angulos de euler para calcular las rotaciones
for i=1:length(time)
    if i==1
        MR_0 = rotmat(uSimpact(1,5:7));% Matriz rotacion inicial, a tiempo 0
    end
    MR_i = rotmat(uSimpact(i,5:7));
        
    MR_02i = MR_i*MR_0';
    VG_02i = rotm2axang(MR_02i);
    uSimpact(i,5:7)= VG_02i(1:3)*VG_02i(4);
end

figure(1);
subplot(2,1,1);hold on;grid on;
plot(time,uSimpact(:,2))
plot(time,uSimpact(:,3))
plot(time,uSimpact(:,4))
subplot(2,1,2);hold on;grid on;
plot(time,fSimpact(:,2))
plot(time,fSimpact(:,3))
plot(time,fSimpact(:,4))

figure(11);
subplot(2,1,1);hold on;grid on;
plot(time,uSimpact(:,5))
plot(time,uSimpact(:,6))
plot(time,uSimpact(:,7))
subplot(2,1,2);hold on;grid on;
plot(time,fSimpact(:,5))
plot(time,fSimpact(:,6))
plot(time,fSimpact(:,7))


figure(2);hold on;grid on;
plot(uSimpact(:,2),fSimpact(:,2))
plot(uSimpact(:,3),fSimpact(:,3))
plot(uSimpact(:,4),fSimpact(:,4))

figure(22);hold on;grid on;
plot(uSimpact(:,5),fSimpact(:,5))
plot(uSimpact(:,6),fSimpact(:,6))
plot(uSimpact(:,7),fSimpact(:,7))



% FUNCIONES AUXILIARES  

function R=rotmat(AE)
psi=AE(1);%*pi/180;
theta=AE(2);%*pi/180;
phi=AE(3);%*pi/180;

Rz1=[cos(psi) sin(psi) 0;-sin(psi) cos(psi) 0;0 0 1];
Rx2=[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)];
Rz3=[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0;0 0 1];
R=Rz3*Rx2*Rz1;
end



function axang = rotm2axang(rotm)
    if ( (size(rotm,1) ~= 3) || (size(rotm,2) ~= 3) )
        error('rotm2axang: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end
    axang   = zeros(4,1);
    epsilon = 1e-12; % min. value to treat a number as zero ...

    % Translate a given rotation matrix R into the corresponding axis-angle representation (u, theta).
    % For further details about the computation, see:
    %   [1] Technical Concepts: Orientation, Rotation, Velocity and Acceleration and the SRM, P. Berner, Version 2.0, 2008,
    %       <http://sedris.org/wg8home/Documents/WG80485.pdf>, pp. 32-33.
    %   [2] A Mathematical Introduction to Robotic Manipulation, Murray & Li & Sastry, CRC Press, 1994, p. 30, eq. (2.17) & (2.18).
    %   [3] Modelling and Control of Robot Manipulators, L. Sciavicco & B. Siciliano, 2nd Edition, Springer, 2008,
    %       p. 35, formula (2.25).
    %   [4] Introduction to Robotics: Mechanics and Control, John J. Craig, 3rd Edition, Pearson/Prentice Hall, 2005,
    %       pp. 47-48, eq. (2.81) & (2.82).
    tr = rotm(1,1) + rotm(2,2) + rotm(3,3);
    if (abs(tr - 3) <= epsilon) % tr = 3 --> theta = 0:
        % Null rotation --> singularity: The rotation matrix R is the identity matrix and the axis of rotation u is undefined.
        % By convention, set u to the default value (0, 0, 1) according to the ISO/IEC IS 19775-1:2013 standard of the Web3D Consortium.
        % See: <http://www.web3d.org/documents/specifications/19775-1/V3.3/Part01/fieldsDef.html#SFRotationAndMFRotation>
        axang(3,1) = 1;
    elseif (abs(tr + 1) <= epsilon) % tr = -1 --> theta = pi:
        if ( (rotm(1,1) > rotm(2,2)) && (rotm(1,1) > rotm(3,3)) )
            u = vertcat(rotm(1,1)+1, rotm(1,2), rotm(1,3));
        elseif (rotm(2,2) > rotm(3,3))
            u = vertcat(rotm(2,1), rotm(2,2)+1, rotm(2,3));
        else
            u = vertcat(rotm(3,1), rotm(3,2), rotm(3,3)+1);
        end
        n = u.'*u;
        axang(1:3,1) = u./sqrt(n); % normalize
        axang(4,1)   = pi;
    else % general case, tr ~= 3 and tr ~= -1:
        axang(4,1) = acos((tr - 1)*0.5); % rotation angle theta within the range (0, pi).
        n_inv = 1/(2*sin(axang(4,1)));
        % unit vector u:
        axang(1,1) = (rotm(3,2) - rotm(2,3))*n_inv;
        axang(2,1) = (rotm(1,3) - rotm(3,1))*n_inv;
        axang(3,1) = (rotm(2,1) - rotm(1,2))*n_inv;
    end
end
