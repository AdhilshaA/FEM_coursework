function fem_periodic_cylinder()
% FEM_PERIODIC_CYLINDER
% 2-D TMz finite-element scattering from a periodic dielectric cylinder array.
% Run this file in MATLAB with no toolbox/PDE add-ons.

    c0 = 299792458;
    mu0 = 4*pi*1e-7;
    eps0 = 1/(mu0*c0^2);

    p.a = 300e-9;
    p.r = 75e-9;
    p.eps_c = 2.5^2;
    p.eps_b = 1;
    p.mu_r = 1;
    p.yb = -750e-9;
    p.yt = 750e-9;
    p.xc = 0;
    p.yc = 0;
    p.loss_c = 0;
    p.wavelengths = (400:10:700)*1e-9;

    outdir = fullfile(fileparts(mfilename('fullpath')), 'figures');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    draw_geometry(p, outdir);

    nx = 61;
    ny = 51;
    [nodes, triangles, edgeInfo] = make_mesh(p, nx, ny);

    draw_mesh(nodes, triangles, p, fullfile(outdir, 'geometry_mesh.png'));
    print_mesh_statistics(nodes, triangles);

    R = zeros(size(p.wavelengths));
    T = zeros(size(p.wavelengths));
    A = zeros(size(p.wavelengths));

    sampleLambda = 550e-9;
    Ez_sample = [];
    Hx_sample = [];
    Hy_sample = [];

    for iw = 1:numel(p.wavelengths)
        lambda0 = p.wavelengths(iw);
        k0 = 2*pi/lambda0;
        omega = 2*pi*c0/lambda0;
        [Aglobal, rhs, area, gradU] = assemble_system(nodes, triangles, edgeInfo, p, k0);
        U = solve_periodic(Aglobal, rhs, nodes);

        [rcoef, tcoef] = extract_spectrum_coefficients(U, nodes, p, k0);
        R(iw) = abs(rcoef)^2;
        T(iw) = abs(tcoef)^2;
        A(iw) = max(0, 1 - R(iw) - T(iw));

        if abs(lambda0-sampleLambda) < eps
            Ez_sample = U;
            [Hx_sample, Hy_sample] = magnetic_field(nodes, triangles, U, omega, mu0*p.mu_r, gradU);
        end
    end

    writetable(table(p.wavelengths(:)*1e9, R(:), T(:), A(:), ...
        'VariableNames', {'Wavelength_nm','Reflection','Transmission','Absorption'}), ...
        fullfile(outdir, 'spectrum.csv'));

    draw_spectrum(p.wavelengths*1e9, R, T, A, fullfile(outdir, 'spectra.png'));

    if isempty(Ez_sample)
        lambda0 = sampleLambda;
        k0 = 2*pi/lambda0;
        omega = 2*pi*c0/lambda0;
        [Aglobal, rhs, ~, gradU] = assemble_system(nodes, triangles, edgeInfo, p, k0);
        Ez_sample = solve_periodic(Aglobal, rhs, nodes);
        [Hx_sample, Hy_sample] = magnetic_field(nodes, triangles, Ez_sample, omega, mu0*p.mu_r, gradU);
    end

    draw_fields(nodes, triangles, Ez_sample, Hx_sample, Hy_sample, p, outdir);

    run_mesh_convergence(p, sampleLambda, outdir);

    energyError = max(abs(R + T + A - 1));
    fprintf('Maximum post-processed R+T+A error: %.3e\n', energyError);
end

function [nodes, triangles, edgeInfo] = make_mesh(p, nx, ny)
    x = linspace(-p.a/2, p.a/2, nx);
    y = linspace(p.yb, p.yt, ny);
    [X, Y] = meshgrid(x, y);
    nodes = [X(:), Y(:)];

    triangles = zeros((nx-1)*(ny-1)*2, 3);
    ir = 0;
    for j = 1:ny-1
        for i = 1:nx-1
            n00 = (j-1)*nx + i;
            n10 = n00 + 1;
            n01 = n00 + nx;
            n11 = n01 + 1;
            ir = ir + 1;
            triangles(ir,:) = [n00 n10 n11];
            ir = ir + 1;
            triangles(ir,:) = [n00 n11 n01];
        end
    end

    tol = max(p.a, p.yt-p.yb)*1e-10;
    left = find(abs(nodes(:,1)+p.a/2) < tol);
    right = find(abs(nodes(:,1)-p.a/2) < tol);
    [~,ordL] = sort(nodes(left,2));
    [~,ordR] = sort(nodes(right,2));
    left = left(ordL);
    right = right(ordR);
    bottom = find(abs(nodes(:,2)-p.yb) < tol);
    top = find(abs(nodes(:,2)-p.yt) < tol);

    edgeInfo.left = left;
    edgeInfo.right = right;
    edgeInfo.bottom = bottom;
    edgeInfo.top = top;
end

function draw_geometry(p, outdir)
    figure('Visible','off','Position',[100 100 900 600]);
    rectangle('Position',[-p.a/2,p.yb,p.a,p.yt-p.yb], ...
        'FaceColor',[0.95 0.98 1],'EdgeColor','k','LineWidth',1.5);
    rectangle('Position',[-p.r,p.r-p.r-p.r,2*p.r,2*p.r], ...
        'Curvature',[1 1],'FaceColor',[0.7 0.85 1],'EdgeColor','b','LineWidth',1.5);
    hold on;
    plot([p.xc p.xc],[p.yb p.yt],'k:');
    quiver(-p.a/3,p.yb-0.18*(p.yt-p.yb),0,0.12*(p.yt-p.yb),0,'b','LineWidth',2,'MaxHeadSize',0.8);
    quiver(-p.a/3,p.yt+0.03*(p.yt-p.yb),0,-0.06*(p.yt-p.yb),0,'r','LineWidth',1.5,'MaxHeadSize',0.8);
    quiver(p.a/3,p.yb-0.03*(p.yt-p.yb),0,0.06*(p.yt-p.yb),0,'g','LineWidth',1.5,'MaxHeadSize',0.8);
    xline(-p.a/2,'--k','Left periodic boundary');
    xline(p.a/2,'--k','Right periodic boundary');
    axis equal;
    xlim([-1.3*p.a/2 1.3*p.a/2]);
    ylim([p.yb-0.15*(p.yt-p.yb) p.yt+0.15*(p.yt-p.yb)]);
    xlabel('x (nm)');
    ylabel('y (nm)');
    xticklabels(xticks*1e9);
    yticklabels(yticks*1e9);
    title('Computational unit cell');
    legend({'Air domain','Dielectric cylinder','Periodic axis','Incident wave','Reflected wave','Transmitted wave'}, ...
        'Location','northoutside','Orientation','horizontal');
    exportgraphics(gcf,fullfile(outdir,'geometry.png'),'Resolution',300);
    close(gcf);
end

function draw_mesh(nodes, triangles, p, filename)
    figure('Visible','off','Position',[100 100 900 650]);
    patch('Faces',triangles,'Vertices',[nodes(:,1)*1e9 nodes(:,2)*1e9], ...
        'FaceColor',[0.95 0.98 1],'EdgeColor',[0.2 0.2 0.2],'LineWidth',0.15);
    view(2); axis equal; grid off;
    hold on;
    th = linspace(0,2*pi,200);
    plot(p.r*cos(th)*1e9,p.r*sin(th)*1e9,'b','LineWidth',2);
    xlabel('x (nm)'); ylabel('y (nm)');
    title('Triangular finite-element mesh');
    exportgraphics(gcf,filename,'Resolution',300);
    close(gcf);
end

function print_mesh_statistics(nodes, triangles)
    nNodes = size(nodes,1);
    nElems = size(triangles,1);
    uniqueEdges = zeros(nElems*3,2);
    k = 0;
    for e = 1:nElems
        te = triangles(e,:);
        es = [te([1 2]); te([2 3]); te([3 1])];
        es = [min(es,[],2) max(es,[],2)];
        uniqueEdges(k+1:k+3,:) = es;
        k = k + 3;
    end
    uniqueEdges = unique(sort(uniqueEdges,2),'rows');
    fprintf('Mesh: %d nodes, %d triangular elements, %d edges\n', nNodes,nElems,size(uniqueEdges,1));
end

function [Aglobal, rhs, area, gradU] = assemble_system(nodes, triangles, edgeInfo, p, k0)
    n = size(nodes,1);
    Aglobal = sparse(n,n);
    rhs = complex(sparse(n,1));

    area = zeros(size(triangles,1),1);
    gradU = zeros(size(triangles,1),2);

    for e = 1:size(triangles,1)
        te = triangles(e,:);
        xy = nodes(te,:);
        x1 = xy(1,1); y1 = xy(1,2);
        x2 = xy(2,1); y2 = xy(2,2);
        x3 = xy(3,1); y3 = xy(3,2);

        detJ = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
        Ae = 0.5*abs(detJ);
        area(e) = Ae;
        if detJ < 0
            te = te([1 3 2]);
            xy = nodes(te,:);
            x1=xy(1,1); y1=xy(1,2);
            x2=xy(2,1); y2=xy(2,2);
            x3=xy(3,1); y3=xy(3,2);
            Ae = 0.5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
        end

        b = [y2-y3, y3-y1, y1-y2];
        c = [x3-x2, x1-x3, x2-x1];
        B = [b/(2*Ae); c/(2*Ae)];

        xc = mean(xy(:,1));
        yc = mean(xy(:,2));
        if (xc-p.xc)^2+(yc-p.yc)^2 <= p.r^2
            epsr = p.eps_c + 1i*p.loss_c;
        else
            epsr = p.eps_b;
        end

        Ke = Ae/(p.mu_r)*(B'*B);
        Me = epsr*Ae/12*[2 1 1;1 2 1;1 1 2];
        Ae_mat = Ke-k0^2*Me;
        for i = 1:3
            for j = 1:3
                Aglobal(te(i),te(j)) = Aglobal(te(i),te(j))+Ae_mat(i,j);
            end
        end
    end

    topSorted = sortrows([nodes(edgeInfo.top,:) edgeInfo.top],1);
    topN = topSorted(:,3);
    for ib = 1:numel(topN)-1
        n1 = topN(ib);
        n2 = topN(ib+1);
        L = norm(nodes(n2,:)-nodes(n1,:));
        Mb = L/6*[2 1;1 2];
        Cb = -1i*k0*Mb;
        te = [n1 n2];
        for i = 1:2
            for j = 1:2
                Aglobal(te(i),te(j)) = Aglobal(te(i),te(j))+Cb(i,j);
            end
        end
    end

    bottomSorted = sortrows([nodes(edgeInfo.bottom,:) edgeInfo.bottom],1);
    bottomN = bottomSorted(:,3);
    for ib = 1:numel(bottomN)-1
        n1 = bottomN(ib);
        n2 = bottomN(ib+1);
        L = norm(nodes(n2,:)-nodes(n1,:));
        Mb = L/6*[2 1;1 2];
        Cb = -1i*k0*Mb;
        te = [n1 n2];
        for i = 1:2
            for j = 1:2
                Aglobal(te(i),te(j)) = Aglobal(te(i),te(j))+Cb(i,j);
            end
        end

        xmid = mean(nodes(te,1));
        uinc = exp(1i*k0*p.yb);
        fb = -1i*k0*L*uinc*[1;1];
        rhs(te) = rhs(te)+fb;
    end

    gradU = zeros(size(triangles,1),2);
end

function U = solve_periodic(Aglobal, rhs, nodes)
    n = size(nodes,1);
    tol = max(max(nodes,[],1)-min(nodes,[],1))*1e-9;
    leftCandidates = find(abs(nodes(:,1)-min(nodes(:,1))) < tol);
    rightCandidates = find(abs(nodes(:,1)-max(nodes(:,1))) < tol);

    [~,ordL] = sort(nodes(leftCandidates,2));
    [~,ordR] = sort(nodes(rightCandidates,2));
    leftCandidates = leftCandidates(ordL);
    rightCandidates = rightCandidates(ordR);

    independent = true(n,1);
    independent(rightCandidates) = false;
    P = sparse(n,sum(independent));
    rnodes = find(independent);
    for r = 1:numel(rnodes)
        P(rnodes(r),r) = 1;
    end
    for j = 1:min(numel(leftCandidates),numel(rightCandidates))
        P(rightCandidates(j), find(rnodes==leftCandidates(j))) = 1;
    end

    Ar = P'*Aglobal*P;
    fr = P'*rhs;
    Ur = Ar\fr;
    U = P*Ur;
end

function [rcoef, tcoef] = extract_spectrum_coefficients(U, nodes, p, k0)
    tol = max(p.a,p.yt-p.yb)*1e-9;
    bottom = find(abs(nodes(:,2)-p.yb)<tol);
    top = find(abs(nodes(:,2)-p.yt)<tol);

    Ebottom = mean(U(bottom));
    Etop = mean(U(top));
    uincBottom = exp(1i*k0*p.yb);
    urefBottom = exp(-1i*k0*p.yb);
    uincTop = exp(1i*k0*p.yt);

    rcoef = (Ebottom-uincBottom)/urefBottom;
    tcoef = Etop/uincTop;
end

function [Hx, Hy] = magnetic_field(nodes, triangles, U, omega, mu, gradU)
    n = size(nodes,1);
    gx = zeros(n,1);
    gy = zeros(n,1);
    count = zeros(n,1);
    for e = 1:size(triangles,1)
        te = triangles(e,:);
        xy = nodes(te,:);
        detJ = (xy(2,1)-xy(1,1))*(xy(3,2)-xy(1,2))-(xy(3,1)-xy(1,1))*(xy(2,2)-xy(1,2));
        Ae = 0.5*abs(detJ);
        if detJ < 0
            te = te([1 3 2]);
            xy = nodes(te,:);
        end
        b = [xy(2,2)-xy(3,2), xy(3,2)-xy(1,2), xy(1,2)-xy(2,2)];
        c = [xy(3,1)-xy(2,1), xy(1,1)-xy(3,1), xy(2,1)-xy(1,1)];
        B = [b/(2*Ae); c/(2*Ae)];
        ue = U(te);
        ge = B*ue;
        gx(te) = gx(te)+real(ge(1));
        gy(te) = gy(te)+real(ge(2));
        count(te) = count(te)+1;
    end
    gx = gx./count;
    gy = gy./count;
    Hx = gx/(1i*omega*mu);
    Hy = -gy/(1i*omega*mu);
end

function draw_spectrum(wavelength_nm, R, T, A, filename)
    figure('Visible','off','Position',[100 100 900 600]);
    plot(wavelength_nm,R,'b-o','LineWidth',1.5,'MarkerSize',3); hold on;
    plot(wavelength_nm,T,'r-s','LineWidth',1.5,'MarkerSize',3);
    plot(wavelength_nm,A,'k-d','LineWidth',1.5,'MarkerSize',3);
    ylim([0 1]); grid on;
    xlabel('Wavelength (nm)');
    ylabel('Normalized power');
    legend('Reflection R','Transmission T','Absorption A','Location','best');
    title('Periodic-cylinder spectral response');
    exportgraphics(gcf,filename,'Resolution',300);
    close(gcf);
end

function draw_fields(nodes, triangles, Ez, Hx, Hy, p, outdir)
    xnm = reshape(nodes(:,1)*1e9,[]);
    ynm = reshape(nodes(:,2)*1e9,[]);

    figure('Visible','off','Position',[100 100 900 650]);
    patch('Faces',triangles,'Vertices',[xnm ynm],'FaceVertexCData',abs(Ez), ...
        'FaceColor','interp','EdgeColor','none');
    axis equal; colorbar; xlabel('x (nm)'); ylabel('y (nm)');
    title('|E_z| at 550 nm');
    exportgraphics(gcf,fullfile(outdir,'electric_field_magnitude.png'),'Resolution',300);
    close(gcf);

    figure('Visible','off','Position',[100 100 900 650]);
    patch('Faces',triangles,'Vertices',[xnm ynm],'FaceVertexCData',real(Ez), ...
        'FaceColor','interp','EdgeColor','none');
    axis equal; colorbar; xlabel('x (nm)'); ylabel('y (nm)');
    title('Re(E_z) at 550 nm');
    exportgraphics(gcf,fullfile(outdir,'electric_field_real.png'),'Resolution',300);
    close(gcf);

    figure('Visible','off','Position',[100 100 900 650]);
    patch('Faces',triangles,'Vertices',[xnm ynm],'FaceVertexCData',sqrt(abs(Hx).^2+abs(Hy).^2), ...
        'FaceColor','interp','EdgeColor','none');
    hold on;
    quiver(nodes(:,1)*1e9,nodes(:,2)*1e9,real(Hx),real(Hy),0.6,'k','LineWidth',0.25);
    axis equal; colorbar; xlabel('x (nm)'); ylabel('y (nm)');
    title('|H| and in-plane magnetic-field vectors at 550 nm');
    exportgraphics(gcf,fullfile(outdir,'magnetic_field.png'),'Resolution',300);
    close(gcf);
end

function run_mesh_convergence(p, lambda0, outdir)
    k0 = 2*pi/lambda0;
    cases = [21 17;31 27;41 37;61 51];
    Tval = zeros(size(cases,1),1);
    ndof = zeros(size(cases,1),1);
    nelem = zeros(size(cases,1),1);

    for ic = 1:size(cases,1)
        [nodes, triangles, edgeInfo] = make_mesh(p,cases(ic,1),cases(ic,2));
        [Aglobal,rhs] = assemble_system(nodes,triangles,edgeInfo,p,k0);
        U = solve_periodic(Aglobal,rhs,nodes);
        [~,tcoef] = extract_spectrum_coefficients(U,nodes,p,k0);
        Tval(ic) = abs(tcoef)^2;
        ndof(ic) = size(nodes,1);
        nelem(ic) = size(triangles,1);
    end

    writetable(table(ndof,nelem,Tval, ...
        'VariableNames',{'Nodes','Triangles','Transmission'}), ...
        fullfile(outdir,'mesh_convergence.csv'));

    figure('Visible','off','Position',[100 100 900 600]);
    plot(ndof,Tval,'b-o','LineWidth',1.5,'MarkerSize',7); grid on;
    xlabel('Number of nodes');
    ylabel('Transmission at 550 nm');
    title('Mesh-convergence study');
    exportgraphics(gcf,fullfile(outdir,'mesh_convergence.png'),'Resolution',300);
    close(gcf);
end
