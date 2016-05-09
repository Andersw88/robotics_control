function ha=plot_plane(A,b,r,col)
%plot_plane(A,b,r,col)
%
%Plots a plane Ax+b=0 in R^3 visualized as a circular disc
%
%A   ... plane unit normal vector in R^3
%b   ... scalar plane offset  
%r   ... scalar radius of the disc visualization
%col ... color 
%
% Robert Krug, 2016/04/20
% robert.krug@oru.se
%

shade=0.25; %opacity of the plane

if nargin==3
col='g'; %color of the plane
end

tol=1e-14;
A=A(:);
% make sure that the normal vector is of unit magnitude
lA=norm(A); A=A/lA; b=b/lA;

% if (b < 0)
%     error('b has to be >= 0');
% end
% if norm(A) < 1-tol || norm(A) > 1+tol
%     error('normal vector A is not normalized'); 
% end

if numel(A) ~=3, error('A has the wrong dimension'); end
% local z axis
z_ax = [0;0;1];
N=100;

a = linspace(-pi,pi,N+1);
a = a(1:end-1); % exclude the last because it is the same as the first

x = r*cos(a);
y = r*sin(a);
z = zeros(1,N);
v = [x;y;z];
rot_ax = cross(z_ax,A);
rot_an = acos(z_ax'*A); %norm(z_ax) = norm(n) = 1 is assumed

if norm(rot_ax) > 1e-15
    rot_ax = rot_ax/norm(rot_ax);
    R = aa2R(rot_an,rot_ax); 
else
    if abs(rot_an) < 1e-15 
        R = eye(3);
    else % rot_an == pi
        R = rx(pi);
    end
end

v = R*v;
n_test = R*z_ax;

if norm(A-n_test) > 1e-14    
    %norm(n-n_test)
    disp('WARNING in rotation test')
end

% normalize
for i=1:N
    % v(:,i) = v(:,i)/norm(v(:,i));
end
v=v';
offset=-b*A;
v=v+repmat(offset',N,1);

% plot related
% ------------------------------------------------------------------
ha(1)=patch(v(:,1),v(:,2),v(:,3),col,'FaceAlpha',shade);

%plot normals
if 0
	hold on;
    s=1; %scale for the normals
    ha(2)=plot3([offset(1) offset(1)+A(1)*s],[offset(2) offset(2)+A(2)*s],[offset(3) offset(3)+A(3)*s],col);
    ha(3)=plot3(offset(1)+A(1)*s,offset(2)+A(2)*s,offset(3)+A(3)*s,strcat(col,'<'));
	hold off;
end



% EOF

