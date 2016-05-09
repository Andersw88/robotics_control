%
% Robert Krug, 2016/04/20
% robert.krug@oru.se
%
classdef MBSVisualizer
	
	properties (Access = protected)
		links_
		vis_scale_
	end
	
	properties  (Access = public)
		h_
	end
	
	methods
		function mbs_visualizer = MBSVisualizer(SP,SV,h)
			if nargin == 2
				mbs_visualizer.h_=figure;
			else
				mbs_visualizer.h_=h;
			end
			
			assert(SP.n >=1);
			
			%Appearence parameter
			mbs_visualizer.vis_scale_=1;
			
			hold on;
			C=SP.C;
			scale=mbs_visualizer.vis_scale_;
			
			%visualize links, joints and CoM's
			for(i=1:SP.n+1)
				R=SV.L(i).R; %rotation of link i
				if i==1
					p=SV.L(i).p; %link 1 is positioned in it's CoM by definition
					c=zeros(3,1); %CoM of link 1 expressed in CF1 is zero
					
				else
					p=SV.L(i).p-R*SP.J(i-1).f; %position of CFi in the world frame
					c=SP.J(i-1).f; %CoM position of link i expressed in link i
				end
				
				g=[];
				%plot CoM
				g(1)=plot3(c(1),c(2),c(3),'gs','MarkerSize',scale*10,'MarkerFaceColor','g','MarkerEdgeColor','k');
				%plot vector IJ to CoM
				g(2)=plot3([0 c(1)],[0 c(2)],[0 c(3)],'b','LineWidth',scale*4);
				
				if i>1
					if 0
						%plot joint
						g(3)=DrawCylinder([0;0;0], [0;0;1], scale/15, scale/2, [1 0 0]);
					else
						%plot joint axis
						a=[0;0;1]*scale/12;
						g(3)=plot3([-a(1) a(1)],[-a(2) a(2)],[-a(3) a(3)],'r','LineWidth',scale*2);
					end
				end
				links(i)=hgtransform('parent',gca);
				set(g,'parent',links(i));
				links(i).Matrix=[R p; zeros(1,3) 1];
				
				if i>SP.n
					continue;
				end
				%visualize vectors from CoM of link C(i+1) to joint i
				g=[];
				if C(i+1)==1
					c=zeros(3,1); %CoM of link 1 expressed in CF1 is zero
				else
					c=SP.J(C(i)).f; %CoM position of link C(i+1) expressed in link C(i+1)
				end
				t=SP.J(i).t; %vector from CoM of link C(i+1) expressed in link C(i+1)
				g(1)=plot3([c(1) c(1)+t(1)],[c(2) c(2)+t(2)],[c(3) c(3)+t(3)],'b','LineWidth',scale*4);
				set(g,'parent',links(C(i+1)));
			end
			
			%visualize end-effectors
			bN=SP.bN;
			bP=SP.bP;
			for i=1:length(bN)
				g=[];
				if bN(i)==1
					c=zeros(3,1); %CoM of link 1 expressed in CF1 is zero
				else
					c=SP.J(bN(i)-1).f; %CoM position of link bN(i) expressed in CF of link bN(i)
					g(1)=plot3(c(1)+bP(1,i),c(2)+bP(2,i),c(3)+bP(3,i),'ko','MarkerSize',scale*5,'MarkerFaceColor','r');
					g(2)=plot3([c(1) c(1)+bP(1,i)],[c(2) c(2)+bP(2,i)],[c(3) c(3)+bP(3,i)],'b','LineWidth',scale*4);
				end
				set(g,'parent',links(bN(i)));
			end
			
			if 1
				%draw the world frame
				g=draw_frame(eye(3),zeros(3,1),'rgb',scale/3);
				set(g,'parent',gca);
			end
			
			mbs_visualizer.links_ = links;
			drawnow; hold off; rotate3d on;
		end
		
		function update(visualizer,SP,SV)
			n=numel(visualizer.links_);
			assert(length(SV.L)==n);
			
			for i=1:n
				R=SV.L(i).R; %rotation of link i
				if i==1
					p=SV.L(i).p; %link 1 is positioned in it's CoM by definition
				else
					p=SV.L(i).p-R*SP.J(i-1).f; %position of CFi in the world frame
				end
				visualizer.links_(i).Matrix=[R p; zeros(1,3) 1];
			end
			drawnow;
		end
	end
end
