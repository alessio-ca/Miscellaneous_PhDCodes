pdem=createpde();

%% Geometry
g = @circleg;
geometryFromEdges(pdem,g);
figure(1) 
pdegplot(g, 'edgeLabels', 'on');
axis equal

%% Problem init
applyBoundaryCondition(pdem,'Edge',1:pdem.Geometry.NumEdges,'u',0);
c = 1;
a = 0;
f = 1;


%% Init mesh
[p,e,t] = initmesh(g,'hmax',1);
figure(2)
pdemesh(p,e,t);
axis equal

%% Refinement
er = Inf;
while er > 0.001
    [p,e,t] = refinemesh(g,p,e,t);
    u = assempde(pdem,p,e,t,c,a,f);
    exact = (1-p(1,:).^2-p(2,:).^2)'/4;
    er = norm(u-exact,'inf');
    fprintf('Error: %e. Number of nodes: %d\n',er,size(p,2));
end
figure(2)
pdemesh(p,e,t);
axis equal

%% Plot solution
figure(3)
pdesurf(p,t,u);
