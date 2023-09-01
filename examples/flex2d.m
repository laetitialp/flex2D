%----------------------------------------------------------------------------
% 
%   4 four-node isoparametric elements of the shear deformable displacement
%   formulation.
%
% Variable descriptions
%%% rheology 
%   kinmtpb = matrix for kinematic equation for bending
%   matmtpb = matrix for material property for bending
%   kinmtps = matrix for kinematic equation for shear
%   matmtps = matrix for material property for shear
%%%% element and global matrices 
%   k = element matrix
%   k  = kb + ks 
%   kb = element matrix for bending stiffness
%   kb = kb1 + kb2 * delta_rho_g_loc
%   ks = element matrix for shear stiffness
%   f = element force vector
%   kk = system stifness matrix
%   ff = system force vector
%   disp = system nodal displacement vector
%   (disp includes vertical(omega)+bending moment thetax and thetay for each nodes)

%%%%  BOOK KEEPING = CONNECTIVITY ARRAYS %%%%%% 
%   gcoord = coordinate values of each node
%   nodes = nodal connectivity of each element
%   index = a vector containing system dofs associated with each element
%   pointb = matrix containing sampling points for bending term
%   weightb = matrix containing weighting coefficients for bending term
%   points = matrix containing sampling points for shear term
%   weights = matrix containing weighting coefficients for shear term
%  
%%%%%  BOUNDARY CONDITIONS : GLOBAL  
%   bcdof = a vector containing dofs associated with boundary conditions
%   bcval = a vector containing boundary condition values associated with
%           the dofs in 'bcdof'
%%%%%  BOUNDARY CONDITIONS : LOCAL  
%  use function bc apply 

%  [BCDOF,BCVAL] = bcapply(nodes_bc,nodes_bc,nodes_bc,0,0,0,BCDOF,BCVAL);
%  fixes all three theta_x; theta_y and omega to 0 on nodes nodes_bc.


%  [BCDOF,BCVAL] = bcapply([],[],nodes_bc,[],[],0,BCDOF,BCVAL);
%  fixes all three omega to 0 on nodes nodes_bc but leave bending moment
%  free (slope is non zero on the side of your model)

%  [BCDOF,BCVAL] = bcapply(nodes_bc,[],nodes_bc,0,[],1000,BCDOF,BCVAL);
%  fixes omega to 1000 on nodes nodes_bc but leave bending moment y 
%  free while fixing bending moment x to 0...  Hope this is clear enough 


% 
%----------------------------------------------------------------------------

clear;
addpath('../src')
fig = 20;
%------------------------------------
%  input data for control parameters Unit SI
%------------------------------------
Lx       = 2000e3;            % size model in x 
Ly       = 1000e3;            % size model in y
ne_dim    = [100,50];        % number of element in x and y 
emodule=7e10;             % elastic modulus  kg*L*s-2/L^2
poisson=0.25;             % Poisson's ratio
% plate thickness 
% put array for a parametric study or a single value if you wanna run only
% one model
TE =[30]*1e3;        
%difference of density between sediment/water/air filling and the mantle or the lower crust

rho_m = 3300; % density mantle
rho_infill = 1600; % density of infill sediments:2000 - water:1000 - air:0
rho_ice    = 2700;
g          = 9.81;
% there is no infill in the area where ice is present.
% here I do a circular Ice sheet 
ice_type = 'halfspace'; 

%ice_type = 'cylinder'; 
%y_ice = 1000e3;
%x_ice = 1000e3;
%Radius_ice = 400e3;

%================================================================



dx = Lx/ne_dim(1);            
dy = Ly/ne_dim(2);
nel       = prod(ne_dim);      % number of elements
nn_dim    = [ne_dim]+1;
nnpe      =4;                  % number of nodes per element
ndof      =3;                  % number of dofs per node
nnode=prod(nn_dim);                 % total number of nodes in system
sdof=nnode*ndof;         % total system dofs
edof=nnpe*ndof;          % degrees of freedom per element
nglxb=2; nglyb=2;        % 2x2 Gauss-Legendre quadrature for bending
nglb=nglxb*nglyb;        % number of sampling points per element for bending
nglxs=1; nglys=1;        % 1x1 Gauss-Legendre quadrature for shear
ngls=nglxs*nglys;        % number of sampling points per element for shear


nnrm    = reshape(1:nnode,nn_dim);
nerm    = reshape(1:nel,ne_dim);
te      = TE(1)*ones(size(nerm));


  
 
%---------------------------------------------
%  gcoord(i,j) nodal coordinate values
%  where i->node no. and j->x or y
%---------------------------------------------
for j = 1:nn_dim(2)
    for i= 1:nn_dim(1)
        gcoord (nnrm(i,j),:) = [0e3+(i-1)*dx,-0e3+(j-1)*dy];
    end
end
Xplot   = reshape(gcoord(:,1),size(nnrm));
Yplot   = reshape(gcoord(:,2),size(nnrm));
milieu_x = (Xplot(1:end-1,1:end-1)+Xplot(2:end,2:end))/2;
milieu_y = (Yplot(1:end-1,1:end-1)+Yplot(2:end,2:end))/2;
te(abs(milieu_y-Ly/2)<75e3)=te(abs(milieu_y-Ly/2)<75e3)*1;
%---------------------------------------------------------
%  nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%  l2g is local2global (l2g(node_number,local_dof)-> global dof number)
%---------------------------------------------------------

for j = 1:ne_dim(2)
    for i= 1:ne_dim(1)
        nodes(nerm(i,j),:)=[nnrm(i,j),nnrm(i+1,j),nnrm(i+1,j+1),nnrm(i,j+1)];
        for k=1:nnpe
            l2g(nerm(i,j),k,1:ndof)=(nodes(nerm(i,j),k)-1)*3+[1:ndof];
        end
    end
end

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------
BCDOF = [];
BCVAL = [];

%bord bottom
nodes_bc = nnrm(:,1)';
%[BCDOF,BCVAL] = bcapply(nodes_bc,nodes_bc,nodes_bc,0,0,0,BCDOF,BCVAL);
[BCDOF,BCVAL] = bcapply([],nodes_bc,[],0,0,0,BCDOF,BCVAL);
%bordtop
nodes_bc = nnrm(:,end)';
%[BCDOF,BCVAL] = bcapply(nodes_bc,nodes_bc,nodes_bc,0,0,0,BCDOF,BCVAL);
[BCDOF,BCVAL] = bcapply([],nodes_bc,[],0,0,0,BCDOF,BCVAL);
%bord left
nodes_bc = nnrm(1,:);
%[BCDOF,BCVAL] = bcapply(nodes_bc,nodes_bc,nodes_bc,0,0,0,BCDOF,BCVAL);
[BCDOF,BCVAL] = bcapply(nodes_bc,[],[],0,0,0,BCDOF,BCVAL);
%bord right
nodes_bc = nnrm(end,:);
[BCDOF,BCVAL] = bcapply(nodes_bc,nodes_bc,nodes_bc,0,0,0,BCDOF,BCVAL);


[bcdof,I,J] = unique(BCDOF);
bcval       = BCVAL(I);

%-------------------------------------
%  input data for ice
%-------------------------------------

bcdofn=[];
bcvaln=[];

switch ice_type 
    case 'cylinder'
        %cylinder ice
% and define the nodes where ice is present
is_there_ice = ((Xplot-x_ice).^2 +  (Yplot-y_ice).^2 ) < Radius_ice^2;
% and assign the constant ice thickness which will only apply in area with ice
thick_ice   = is_there_ice*1500;
case 'halfspace'
        %cylinder ice
% 
is_there_ice = ((Xplot-Lx/4)<0 & abs(Yplot-Ly/2)>75e3)| ((Xplot-Lx/6)<0);
% and assign the constant ice thickness which will only apply in area with ice
thick_ice   = 3000.*is_there_ice;
%thick_ice   = atan((Xplot-Lx/6)/200e3*2*pi)*3000*2/pi.*is_there_ice;
pcolor(Xplot,Yplot,thick_ice);pause(0.5)

    case 'from_file'
   % providing you have a file with x,y,z for your ice
   x = rand(1000,1);
   y = rand(1000,1);
   z = 0.3-sqrt((y-0.5).^2*4 +(x-0.5).^2*2);
   x = x*Lx; y=y*Ly; z=z*5000;
   z(z<0)=0;
   % 
   F = TriScatteredInterp(x,y,z);
   thick_ice = F(Xplot,Yplot); 
   thick_ice(isnan(thick_ice))=0; %replace NaN in the place with no ice interpolated data by zero     
   is_there_ice = thick_ice >0; 
   pcolor(Xplot,Yplot,thick_ice);pause(0.5)
    otherwise
        error('you did not provide a method for ice')
end

forcevalues = -rho_ice*thick_ice*g*dx*dy;
forcevalues(:,1)  = forcevalues(:,1)/2; 
forcevalues(:,end)= forcevalues(:,end)/2;
forcevalues(1,:)  = forcevalues(1,:)/2; 
forcevalues(end,:)= forcevalues(end,:)/2;

deltarhog = ~is_there_ice*((rho_m-rho_infill)*g)+is_there_ice*(rho_m*g);
forcedof = 3+3*(-1+find(forcevalues(:)));
forceval = forcevalues(find(forcevalues(:)));

bcdofn=[bcdofn;forcedof(:)];
bcvaln=[bcvaln;forceval(:)];



for nte = 1:length(TE)
t = TE(nte);
D = fematiso(1,emodule,poisson);
matmtpb=D*t^3/12;  % bending material property
shearm=0.5*emodule/(1.0+poisson);          % shear modulus
shcof=5/6;                                 % shear correction factor
matmtps=shearm*shcof*t^3/12*[1 0; 0 1]; %   % shear material property

%----------------------------------------------
%  initialization of matrices and vectors
%----------------------------------------------

ff=zeros(sdof,1);       % system force vector
kk=sparse(sdof,sdof);    % system matrix
disp=zeros(sdof,1);     % system displacement vector
index=zeros(edof,1);    % index vector
kinmtpb=zeros(3,edof);   % kinematic matrix for bending

kinmtps=zeros(2,edof);   % kinematic matrix for shear


%-----------------------------------------------------------------
%  computation of element matrices and vectors
%-----------------------------------------------------------------
%
%  for bending stiffness
%
[pointb,weightb]=feglqd2(nglxb,nglyb);     % sampling points & weights
%
%  for shear stiffness
%
[points,weights]=feglqd2(nglxs,nglys);     % sampling points & weights



for i=1:nnpe
    %nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
    nd(i)=nodes(1,i);         % extract connected node for (iel)-th element
    xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
    ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
end

k=zeros(edof,edof);         % initialization of element matrix to zero
kb1=zeros(edof,edof);        % initialization of bending matrix to zero
kb2=zeros(edof,edof);        % initialization of bending matrix to zero
ks=zeros(edof,edof);        % initialization of shear matrix to zero
kb1s=zeros(edof,edof);        % initialization of bending matrix to zero
kb2s=zeros(edof,edof);
%------------------------------------------------------
%  numerical integration for bending term
%------------------------------------------------------
ip = 1
for intx=1:nglxb
    x=pointb(intx,1);                  % sampling point in x-axis
    wtx=weightb(intx,1);               % weight in x-axis
    for inty=1:nglyb
        y=pointb(inty,2);                  % sampling point in y-axis
        wty=weightb(inty,2) ;              % weight in y-axis
        
        [shape,dhdr,dhds]=feisoq4(x,y);     % compute shape functions and
        % derivatives at sampling point
        
        jacob2=fejacob2(nnpe,dhdr,dhds,xcoord,ycoord);  % compute Jacobian
        
        detjacob=det(jacob2);                 % determinant of Jacobian
        invjacob=inv(jacob2);                 % inverse of Jacobian matrix
        
        % derivatives w.r.t.
        dhdX = invjacob*[dhdr;dhds]; % physical coordinate
        
        kinmtpb=fekinepb(nnpe,dhdX(1,:),dhdX(2,:));          % bending kinematic matrix
        
        %--------------------------------------------
        %  compute bending element matrix
        %--------------------------------------------
        Load            = zeros(1,12);
        Load(1,3:3:end) = shape;
        
        %
        kb1=kb1+wtx*wty*detjacob*(kinmtpb'*matmtpb*kinmtpb);
        kb2=kb2+wtx*wty*detjacob*(Load'*Load);
        kb1s=kb1s+wtx*wty*detjacob*(kinmtpb'*(matmtpb/t^3)*kinmtpb);
       
        %kb=kb+kinmtpb'*matmtpb*kinmtpb*wtx*wty*detjacob;
        
        curv(ip,:,:)    = D*kinmtpb;
        ip = ip+1;
        
    end
end                      % end of numerical integration loop for bending term
trucb=kb1s-kb1/t^3;
%------------------------------------------------------
%  numerical integration for bending term
%------------------------------------------------------
ks2 = ks;
for intx=1:nglxs
    x=points(intx,1);                  % sampling point in x-axis
    wtx=weights(intx,1);               % weight in x-axis
    for inty=1:nglys
        y=points(inty,2);                  % sampling point in y-axis
        wty=weights(inty,2) ;              % weight in y-axis
        
        [shape,dhdr,dhds]=feisoq4(x,y);     % compute shape functions and
        % derivatives at sampling point
        
        jacob2=fejacob2(nnpe,dhdr,dhds,xcoord,ycoord);  % compute Jacobian
        
        detjacob=det(jacob2);                 % determinant of Jacobian
        invjacob=inv(jacob2);                 % inverse of Jacobian matrix
        
        dhdX = invjacob*[dhdr;dhds]
        dhdx = dhdX(1,:);
        dhdy = dhdX(2,:);
        % shear kinematic matrix
        i = 1:nnpe;
        i1=(i-1)*3+1;
        i2=i1+1;
        i3=i2+1;
        kinmtps(1,i1)=-shape;
        kinmtps(:,i3)=dhdX;
        kinmtps(2,i2)=-shape;
        
        %----------------------------------------
        %  compute shear element matrix
        %----------------------------------------
        
        ks=ks+kinmtps'*matmtps*kinmtps*wtx*wty*detjacob;
        ks2=ks2+kinmtps'*(matmtps/t^3)*kinmtps*wtx*wty*detjacob;
    end
end                      % end of numerical integration loop for shear term
trucs=ks2-ks/t^3;

%--------------------------------
%  element assembly
%--------------------------------

for iel=1:nel
    ks  = (ks2-trucs)*te(iel)^3+trucs;  
    kb1 = (kb1s-trucb)*te(iel)^3+trucb;
    k=kb1+kb2*mean(deltarhog(nodes(iel,:)))+ks;
    index=feeldof(nodes(iel,:),nnpe,ndof);% extract system dofs associated with element
    kk(index,index)=kk(index,index)+k;% assemble element matrices
    
end

%-----------------------------
%   apply boundary conditions
%-----------------------------

ff(bcdofn)=bcvaln;
kk(bcdof,:)=0;
kk(bcdof,bcdof)=eye(length(bcdof));
ff(bcdof)=bcval;
%----------------------------
%  solve the matrix equation
%----------------------------

disp=kk\ff;

%%% compute stresses is slow so if you dont want to plot stress no reason to
%%% make an other loop ! 
% %%%%%%%%%%%%%%%%%%%%%COMPUTE STRESSS%%%%%%%%%%%%
%  nip = nglxb*nglyb;
% for iel=1:prod(ne_dim)                  % loop for the total number of elements
%         index    = feeldof(nodes(iel,:),nnpe,ndof);% extract system dofs associated with element
%         disp_loc = disp(index);
%         SB      = zeros(3,1);
%
%
%         for ip=1:nip
%             stressb = squeeze(curv(ip,:,:))*disp_loc;
%             SB      = SB+stressb;
%         end                      % end of numerical integration loop for bending term
%
%         STRESS(iel,:)=SB/4;
%     end
%
%
%




num=1:1:sdof;
displace=[num' disp];                  % print nodal displacements
thetax  = reshape(disp(1:3:end),size(nnrm));
thetay  = reshape(disp(2:3:end),size(nnrm));
omega   = reshape(disp(3:3:end),size(nnrm));

fh = figure(1+nte*10);
set(fh,'colormap',jet);
contour(Xplot/1000,Yplot/1000,omega);colorbar('horiz');
shading flat;axis equal tight;
fh = figure(2+nte*10);
set(fh,'colormap',jet);
surf(Xplot/1000,Yplot/1000,omega+thick_ice);colorbar('horiz');
shading flat;axis equal tight;



deflect(nte,:)=omega(:,51);
end
figure(1); plot(Xplot(:,51)/1000,deflect');


