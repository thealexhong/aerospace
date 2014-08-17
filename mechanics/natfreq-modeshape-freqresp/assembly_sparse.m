% Alexander Hong (997584706)
% November 20, 2013

% assembly.m
% This code assembles the global stiffness and mass matrices in sparse
... matrix format

function [K, M] = assembly_sparse (EA, EI, RHO, nElements, ...
    coords, connec, connec_dof, nDOF)
% Variables
% EA: axial rigidity
% EI: flexural rigidity
% RHO: mass per unit length
% nElements: number of elements
% coods: nodal coordinates
% connec: connectivity matrix
% connec_dof: relationship between local and global DOF numbering
% nDOF: number of DOF


%% Assemble local stiffness and mass matrices
K = zeros(nDOF, nDOF);
M = zeros(nDOF, nDOF);
nDOF_element = size(connec_dof,2);
K_element = {};
M_element = {};
for i = 1:nElements
    X1 = coords(connec(i,1),1);
    Y1 = coords(connec(i,1),2);
    X2 = coords(connec(i,2),1);
    Y2 = coords(connec(i,2),2);
    [K_element{i}] = ElementStiffMat(EA, EI, X1, Y1, X2, Y2);
    [M_element{i}] = ElementMassMat(RHO, X1, Y1, X2, Y2);
end

%% Assemble global stiffness and mass matrices in sparse matrix format
% row and column indices
row_col{1} = []; % for stiffness matrix
row_col{2} = []; % for mass matrix

% global stiffness and mass matrices in sparse matrix form
K = [];
M = [];
for i = 1:nElements
    for j = 1:6
        for k = 1:6
            % Assemble the global stiffness matrix
            if (K_element{i}(j,k) ~= 0)
                index1 = 0;
                sym_flag1 = 0; % notes if the current value is symmetric
                ip = connec_dof(i,j);
                jp = connec_dof(i,k);
                [row1 col1] = size(row_col{1});
                for l = 1:col1
                    if(ip ~= jp && row_col{1}(1,l) == jp...
                            && row_col{1}(2,l) == ip)                        
                            sym_flag1 = 1;
                            break;
                    end
                    if (row_col{1}(1,l) == ip...
                            && row_col{1}(2,l) == jp)
                        index1 = l;
                        break;
                    end
                end
                
                if (index1 == 0 &&  sym_flag1 == 0)                  
                    col1 = [ip; jp];
                    row_col{1} = [row_col{1} col1];                 
                    K = [K K_element{i}(j,k)];
                end
                if (index1 ~= 0 && sym_flag1 == 0)
                    K(index1) = K(index1) + K_element{i}(j,k);
                end
            end          
        end
    end
end
K = [K; row_col{1}];
% insertion into matrix form
K_tmp = zeros(nDOF, nDOF);
for m = 1:size(K,2)
    K_tmp (K(2,m),K(3,m)) = K(1,m); 
end
K = K_tmp;

for i = 1:nElements
    for j = 1:6
        for k = 1:6
            % Assemble the global mass matrix
            if (M_element{i}(j,k) ~= 0)
                index2 = 0;
                sym_flag2 = 0;
                ip = connec_dof(i,j);
                jp = connec_dof(i,k);
                [row2 col2] = size(row_col{2});
                for l = 1:col2
                    if(ip ~= jp && row_col{2}(1,l) == jp...
                            && row_col{2}(2,l) == ip)                        
                            sym_flag2 = 1;
                            break;
                    end
                    if (row_col{2}(1,l) == ip...
                            && row_col{2}(2,l) == jp)
                        index2 = l;
                        break;
                    end
                end
                if (index2 == 0 &&  sym_flag2 == 0)                  
                    col2 = [ip; jp];
                    row_col{2} = [row_col{2} col2];                 
                    M = [M M_element{i}(j,k)];
                end
                if (index2 ~= 0 && sym_flag2 == 0)
                    M(index2) = M(index2) + M_element{i}(j,k);
                end
            end
        end
    end
end
M = [M; row_col{2}];
% insertion into matrix form
M_tmp = zeros(nDOF, nDOF);
for m = 1:size(M,2)
    M_tmp (M(2,m),M(3,m)) = M(1,m); 
end
M = M_tmp;