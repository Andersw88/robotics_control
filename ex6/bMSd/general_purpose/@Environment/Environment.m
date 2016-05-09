%
% Robert Krug, 2016/04/20
% robert.krug@oru.se
%
% Generates a planar environment which reacts with a force
% f_e=-k_e*dz-c_e*ddz
%
% The environment plane is defnied by the x/y plane of cf_e_
%
classdef Environment<handle
    
    properties (Access = protected)
        vis_scale_
        g_

    end
    
    properties  (Access = public)
        h_
        k_e_ %environment stiffness
        c_e_ %environment damping
        cf_e_ %environment frame (represented as a hgtransform)
        fn_
    end
    
    methods
        function environment = Environment(cf_e,h)
%             if nargin == 1
%                 environment.h_=figure;
%             else
%                 environment.h_=h;
%             end
            scale=1; %define visualization scaling factor
            
            %guesstimate environment stiffness/damping parameters (k=c^2/4 yields critical damping)
            environment.c_e_=150;
            environment.k_e_=environment.c_e_^2/4;
            environment.fn_=zeros(3,1);
            
            environment.vis_scale_=scale;
            
            phi = linspace(-pi,pi,100);
            x = zeros(1,100); y = scale*cos(phi); z = scale*sin(phi); 
            
            %plot plane with normal along x
            h(1)=patch(x,y,z,'c','FaceAlpha',0.25);
            h(2:7)=draw_frame(eye(3),zeros(3,1),'rgb',scale/3);
            
            environment.cf_e_=hgtransform('parent',gca);
            set(h,'parent',environment.cf_e_);
            environment.cf_e_.Matrix=cf_e;
            
            environment.h_=[environment.h_; h(:)];
        end
        
        function [fn dfn]=getContactForce(environment,p_e,v_e,dt)
            %generate plane H/K representation (a'x+b=0)
            a=environment.cf_e_.Matrix(1:3,1);
            b=-a'*environment.cf_e_.Matrix(1:3,4);
            
            %project end-effector positon & velolcity on the environment normal
            dx=a'*p_e+b;
            ddx=a'*v_e;
            
            %frictionless environment -> torque parts are always zero
            if dx > 0
                fn=a*(-dx*environment.k_e_-ddx*environment.c_e_);
            else
                fn=zeros(3,1);
            end
            
            %compute force derivative
            dfn=(environment.fn_-fn)/dt;
            environment.fn_=fn;
            
            
        end
    end
end
