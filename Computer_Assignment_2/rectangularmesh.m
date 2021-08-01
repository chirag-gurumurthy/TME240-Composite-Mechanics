%% function [mesh, coord, Edof] = rectangularmesh(xmin, xmax, ymin, ymax, nelx, nely)
%----------------------------------------------------------------------------
% INPUT:
%          xmin,xmax: min/max x-coordinates of planar plate in xy-plane
%          ymin,ymax: min/max y-coordinates of planar plate in xy-plane
%
%          nelx: number of elements along the x-direction
%          nely: number of elements along the y-direction
%
% OUTPUT:
%
%          mesh: array where each column corresponds to an element and
%                where the rows indicates the nodes associated with that
%                element.
%
%          coord: Nnode x 2 array containing the x and y coordinates of
%                 each node (Nnode = total number of nodes)
%
%          Edof: Element topology matrix
%                Each column corresponds to the associated degrees-of-freedom
%                (5 dofs per node)
%
%----------------------------------------------------------------------------

function [mesh, coord, Edof] = rectangularmesh(xmin, xmax, ymin, ymax, nelx, nely)
j=1;
k=1;
for i=1:nelx*nely
    mesh(:,i)=[(j-1)+nelx*j-(nelx-k);(j-1)+nelx*j+1-(nelx-k);(nelx+1)*(j+1)-(nelx-k);(nelx+1)*(j+1)-1-(nelx-k)];
    k=k+1;
    if (k>nelx)
        k=1;
        j=j+1;
        if(j>nely)
            break
        end
    end
end

[X,Y]=meshgrid(linspace(xmin,xmax,nelx+1),linspace(ymin,ymax,nely+1));
X=X';
Y=Y';
coord=[X(:),Y(:)];
nodedof=5;
Edof=zeros(nelx*nely,nodedof*4+1);
Edof(:,1)=[[1:1:(nelx)*(nely)]'];
for i=1:nodedof
Edof(:,[i+1:nodedof:4*(nodedof)+1])=mesh'*nodedof-(nodedof-i);
end
Edof = Edof(:,2:end).';
% permute here instead of confusing students...
Edof = Edof([1, 2, 6, 7, 11, 12, 16, 17, 3, 8, 13, 18, 5, 4, 10, 9, 15, 14, 20, 19], :);
