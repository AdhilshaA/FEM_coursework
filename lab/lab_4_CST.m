%% LAB 4: CONSTANT STRAIN TRIANGLE (CST) ELEMENT AND FEM ASSEMBLY
% Units: N, mm, MPa
% Plane-stress CST formulation.
%
% This script solves:
%   1. A single CST element
%   2. A two-element CST patch test
%   3. A structured-mesh cantilever beam with mesh refinement
%
% The beam dimensions in the PDF are converted from metres to millimetres
% so that they are consistent with MPa (N/mm^2), N, and mm.

clear; clc; close all;
format short g;

%% MATERIAL PROPERTIES
E = 200000;             % Young's modulus (MPa)
nu = 0.30;              % Poisson's ratio
t = 1;                  % Thickness (mm)
D = plane_stress_matrix(E, nu);

%% 1. SINGLE CST ELEMENT
fprintf('\n====================================================\n');
fprintf('1. SINGLE CST ELEMENT\n');
fprintf('====================================================\n');

nodes1 = [0 0; 1 0; 0 1];
element1 = [1 2 3];
[A1, B1, D1, Ke1] = cst_matrices(nodes1(element1,:), E, nu, t);

F1 = [0; 0; 0; 0; 0; -100];
fixedDofs1 = [1 2 4];
U1 = solve_with_constraints(Ke1, F1, fixedDofs1);
R1 = Ke1*U1 - F1;
strain1 = B1*U1;
stress1 = D1*strain1;

fprintf('Area A (mm^2):\n');
disp(A1);
fprintf('B matrix (1/mm):\n');
disp(B1);
fprintf('D matrix (MPa):\n');
disp(D1);
fprintf('Element stiffness Ke (N/mm):\n');
disp(Ke1);
fprintf('Nodal displacement vector [u1 v1 u2 v2 u3 v3] (mm):\n');
disp(U1);
fprintf('Support reactions (N):\n');
disp(R1);
fprintf('Strain [epsilon_x epsilon_y gamma_xy]:\n');
disp(strain1);
fprintf('Stress [sigma_x sigma_y tau_xy] (MPa):\n');
disp(stress1);

figure('Name','Single CST Element');
patch(nodes1(:,1), nodes1(:,2), [0.65 0.82 1.00], ...
    'FaceAlpha', 0.45, 'EdgeColor', 'b', 'LineWidth', 1.5);
hold on;
scale1 = 200;
deformedNodes1 = nodes1 + scale1*reshape(U1, 2, []).';
patch(deformedNodes1(:,1), deformedNodes1(:,2), [1.00 0.70 0.45], ...
    'FaceAlpha', 0.45, 'EdgeColor', 'r', 'LineWidth', 1.5);
plot(nodes1(:,1), nodes1(:,2), 'bo', 'MarkerFaceColor', 'b');
plot(deformedNodes1(:,1), deformedNodes1(:,2), 'ro', 'MarkerFaceColor', 'r');
axis equal; grid on;
xlabel('x (mm)'); ylabel('y (mm)');
title('Single CST Element: Original and Deformed (magnified)');
legend('Original', 'Deformed', 'Location', 'best');

%% 2. CST PATCH TEST
fprintf('\n====================================================\n');
fprintf('2. CST PATCH TEST\n');
fprintf('====================================================\n');

% setting up a 1x1 square patch with two triangular elements
nodes2 = [0 0; 1 0; 1 1; 0 1];
elements2 = [1 2 3; 1 3 4];
nNodes2 = size(nodes2, 1);
nDofs2 = 2*nNodes2;

% Assemble the global stiffness matrix for the patch test
[K2, elementData2] = assemble_global_stiffness(nodes2, elements2, E, nu, t);

% Prescribed displacement field: u = 0.001*x, v = 0.002*y.
U2 = zeros(nDofs2, 1);
U2(1:2:end) = 0.001*nodes2(:,1); % setting u displacements
U2(2:2:end) = 0.002*nodes2(:,2); % setting v displacements
exactStrain2 = [0.001; 0.002; 0];

fprintf('Prescribed nodal displacement vector (mm):\n');
disp(U2);
fprintf('Element results [element area eps_x eps_y gamma_xy sigma_x sigma_y tau_xy]:\n');

patchResults = zeros(size(elements2,1), 8);
for e = 1:size(elements2,1)

    % select the degrees of freedom for this element
    dofs = element_dofs(elements2(e,:));

    % Extract the displacements for this element
    ue = U2(dofs);

    % Compute the strain and stress for this element
    strain = elementData2(e).B*ue;
    stress = elementData2(e).D*strain;

    % Store the results for this element
    patchResults(e,:) = [e, elementData2(e).area, strain.', stress.'];
    fprintf('%d  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f\n', ...
        patchResults(e,:));
end

strainError2 = patchResults(:,3:5) - repmat(exactStrain2.', size(elements2,1), 1);
fprintf('Maximum absolute strain error: %.3e\n', max(abs(strainError2(:))));

figure('Name','CST Patch Test');
triplot(elements2, nodes2(:,1), nodes2(:,2), 'k-', 'LineWidth', 1.5);
hold on;
deformedNodes2 = nodes2 + 200*reshape(U2, 2, []).';
triplot(elements2, deformedNodes2(:,1), deformedNodes2(:,2), ...
    'b-', 'LineWidth', 1.5);
axis equal; grid on;
xlabel('x (mm)'); ylabel('y (mm)');
title('CST Patch Test: Original and Deformed (magnified)');
legend('Original mesh', 'Deformed mesh', 'Location', 'best');

%% 3. CANTILEVER BEAM USING CST ELEMENTS
fprintf('\n====================================================\n');
fprintf('3. CANTILEVER BEAM USING CST ELEMENTS\n');
fprintf('====================================================\n');

% dimensions 
L = 1000;  % mm
H = 100;   % mm
t = 1;  % mm

% load 
Fy = -100;  % N

% Material properties
E = 200000;  % MPa
nu = 0.30; 
I = t*H^3/12; % Moment of inertia for a rectangular cross-section


% different mesh cases for convergence study
meshCases = [10 2; 20 4; 40 8; 80 16];
nCases = size(meshCases, 1); % no. of mesh cases

% initialize results storage
beamResults = zeros(nCases, 8);
beamSolutions = cell(nCases, 1);

% calculate the analytical tip displacement magnitude using Euler-Bernoulli beam theory
analyticalTipMagnitude = abs(Fy)*L^3/(3*E*I);

for iCase = 1:nCases
    % generate the structured mesh
    Nx = meshCases(iCase,1);
    Ny = meshCases(iCase,2);
    [nodes, elements, nodeId] = structured_mesh(Nx, Ny, L, H);

    % assemble the global stiffness matrix
    [K, ~] = assemble_global_stiffness(nodes, elements, E, nu, t);

    % set no. of nodes and initialize the global force vector
    nNodes = size(nodes,1);
    F = zeros(2*nNodes, 1);


    tipNode = nodeId(Nx, floor(Ny/2)); % find node ID of the tip node at the free end
    tipVDof = 2*tipNode; % vertical degree of freedom for the tip node
    F(tipVDof) = Fy; % apply the vertical load at the tip node

    % fix the nodes at the clamped end
    fixedNodes = zeros(Ny+1, 1);
    for iy = 0:Ny
        fixedNodes(iy+1) = nodeId(0, iy);
    end
    fixedDofs = reshape([2*fixedNodes-1, 2*fixedNodes].', [], 1);

    % solve for the displacements with the applied constraints
    U = solve_with_constraints(K, F, fixedDofs);

    % isolate the tip displacement and compute the relative error
    tipDisplacement = U(tipVDof);
    relativeError = abs(abs(tipDisplacement) - analyticalTipMagnitude) ...
        / analyticalTipMagnitude * 100;

    % store the results for this mesh case
    beamResults(iCase,:) = [Nx, Ny, nNodes, size(elements,1), ...
        tipDisplacement, abs(tipDisplacement), ...
        analyticalTipMagnitude, relativeError];
    beamSolutions{iCase} = struct('nodes', nodes, 'elements', elements, ...
        'U', U, 'tipNode', tipNode);
end

fprintf(['Nx  Ny  Nodes  Elements  FEM tip v (mm)  FEM magnitude (mm)  ', ...
    'Analytical (mm)  Error (%%)\n']);
fprintf('%2d  %2d  %5d  %8d  %16.6f  %19.6f  %16.6f  %8.3f\n', ...
    beamResults.');

fprintf('\nSecond moment of area I = %.6f mm^4\n', I);
fprintf('Analytical tip displacement magnitude = %.6f mm\n', ...
    analyticalTipMagnitude);

% Plot the coarsest mesh and its magnified deformation.
coarse = beamSolutions{1};
scaleBeam = 100;
deformedCoarse = coarse.nodes + scaleBeam*reshape(coarse.U, 2, []).';

figure('Name','Cantilever CST Mesh');
triplot(coarse.elements, coarse.nodes(:,1), coarse.nodes(:,2), ...
    'k--', 'LineWidth', 1);
hold on;
triplot(coarse.elements, deformedCoarse(:,1), deformedCoarse(:,2), ...
    'r-', 'LineWidth', 1.2);
axis equal; grid on;
xlabel('x (mm)'); ylabel('y (mm)');
title('Cantilever: Original and Deformed Mesh (magnified)');
legend('Original mesh', 'Deformed mesh', 'Location', 'best');

figure('Name','CST Cantilever Convergence');
plot(beamResults(:,4), beamResults(:,6), 'o-', 'LineWidth', 1.5);
hold on;
yline(analyticalTipMagnitude, 'r--', 'Analytical');
grid on;
xlabel('Number of triangular elements');
ylabel('Tip displacement magnitude (mm)');
title('CST Cantilever Mesh Refinement');
legend('CST FEM', 'Euler-Bernoulli analytical', 'Location', 'best');

% Save the convergence table for the lab report.
resultsTable = array2table(beamResults, ...
    'VariableNames', {'Nx','Ny','Nodes','Elements','FEM_tip_v_mm', ...
    'FEM_magnitude_mm','Analytical_mm','Relative_error_percent'});
writetable(resultsTable, 'Lab4_CST_Cantilever_Convergence.csv');
disp(resultsTable);

fprintf('\nLab 4 CST analysis completed successfully.\n');
fprintf(['The CST model is overly stiff on coarse bending meshes because ', ...
    'strain is constant in each triangle.\n']);

%% LOCAL FUNCTIONS

function D = plane_stress_matrix(E, nu)
    D = E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
end

function [A, B, D, Ke] = cst_matrices(coords, E, nu, thickness)

    % separate x and y coordinates
    x = coords(:,1);
    y = coords(:,2);

    % calculate the area of the triangle
    twoA = x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2));
    A = abs(twoA)/2;

    if A <= 0
        error('CST element has zero area.');
    end

    % calculate the B matrix
    beta = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
    gamma = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];
    B = [beta(1) 0 beta(2) 0 beta(3) 0; ...
         0 gamma(1) 0 gamma(2) 0 gamma(3); ...
         gamma(1) beta(1) gamma(2) beta(2) gamma(3) beta(3)]/(2*A);

    % calculate the D matrix and element stiffness matrix
    D = plane_stress_matrix(E, nu);
    Ke = thickness*A*(B.'*D*B);
end

function dofs = element_dofs(conn)
    % conn is a vector of node indices for the element

    % select the degrees of freedom for the element
    dofs = reshape([2*conn-1; 2*conn], [], 1);
end

function [K, elementData] = assemble_global_stiffness(nodes, elements, E, nu, thickness)
    
    % setting up the number of degrees of freedom for the global stiffness matrix
    nDofs = 2*size(nodes,1);

    % initialize the global stiffness matrix and element data structure
    K = sparse(nDofs, nDofs);
    elementData = repmat(struct('area', [], 'B', [], 'D', [], 'Ke', []), ...
        size(elements,1), 1);

    for e = 1:size(elements,1)
        % select the nodes for this element
        conn = elements(e,:);

        % compute the element matrices using the CST formulation
        [A, B, D, Ke] = cst_matrices(nodes(conn,:), E, nu, thickness);

        % assemble the global stiffness matrix and store element data
        dofs = element_dofs(conn);
        K(dofs,dofs) = K(dofs,dofs) + Ke;

        % set the element data for this element
        elementData(e).area = A;
        elementData(e).B = B;
        elementData(e).D = D;
        elementData(e).Ke = Ke;
    end
end

function U = solve_with_constraints(K, F, fixedDofs)

    allDofs = (1:length(F)).'; % all degrees of freedom
    freeDofs = setdiff(allDofs, fixedDofs); % determine the free degrees of freedom
    U = zeros(size(F)); % initialize the displacement vector
    U(freeDofs) = K(freeDofs,freeDofs)\F(freeDofs); % solve for the free displacements
end

function [nodes, elements, nodeId] = structured_mesh(Nx, Ny, L, H)

    % generate a structured mesh of triangular elements over a rectangular domain
    x = linspace(0, L, Nx+1);
    y = linspace(0, H, Ny+1);
    [X, Y] = meshgrid(x, y);

    % create the node coordinates and a function to map (ix, iy) to node ID
    nodes = [X(:), Y(:)];
    nodeId = @(ix, iy) iy + 1 + ix*(Ny+1);

    % create the triangular elements by splitting each rectangle into two triangles
    elements = zeros(2*Nx*Ny, 3);
    e = 0;
    for ix = 0:Nx-1
        for iy = 0:Ny-1
            n1 = nodeId(ix, iy);
            n2 = nodeId(ix+1, iy);
            n3 = nodeId(ix+1, iy+1);
            n4 = nodeId(ix, iy+1);
            e = e + 1;
            elements(e,:) = [n1 n2 n3];
            e = e + 1;
            elements(e,:) = [n1 n3 n4];
        end
    end
end
