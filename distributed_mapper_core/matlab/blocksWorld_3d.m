clear all; close all; clc; 
% figure_format; %Change the default figure formats

import gtsam.*

%% Configurations
sigma_x = 0.2; % X Covariance
sigma_y = 0.2; % Y Covariance
sigma_z = 0.2; % Z Covariance
sigma_yaw = 5*pi/180; % Yaw Covariance
sigma_pitch =5*pi/180; % Pitch Covariance
sigma_roll = 5*pi/180; % Roll Covariance
priorNoise = noiseModel.Diagonal.Sigmas([0.05; 0.05; 0.05; 5*pi/180; 5*pi/180; 5*pi/180]);
p0 = Pose3;
betweenNoise = noiseModel.Diagonal.Sigmas([1; 1; 1; 1; 1; 1; ]); %sigma_x; sigma_y; sigma_z; sigma_yaw; sigma_pitch; sigma_roll]);  
upperTriangularNoise = [sigma_x, 0, 0, 0, 0, 0, ...
                        0, sigma_y, 0, 0, 0, 0, ...
                        0, 0, sigma_z, 0, 0, 0, ...
                        0, 0, 0, sigma_yaw, 0, 0, ...
                        0, 0, 0, 0, sigma_pitch, 0, ...
                        0, 0, 0, 0, 0, sigma_roll]; 
width = 4; % Dimensions (nr of nodes) of block world  -- keep it even
height = 4; % Dimensions (nr of nodes) of block world  -- keep it even
depth = 2; % Dimensions (nr of nodes) of block world  -- keep it even
tile_height = 1; % Height of each tile -- num_edges -- keep it odd
tile_width = 1; % Width of each tile -- num_edges -- keep it odd
tile_depth = 1; % Depth of each tile -- num_edges -- keep it odd
num_subgraphs = (width)/(tile_width+1) * (height)/(tile_height+1) * (depth)/(tile_depth+1); 
graph_full = []; %For saving the full graph to disk
initial_full = []; %For saving the full initial to disk

% Checks
if mod(width,2) == 1
  error('width must be off')
end

if rem(width,tile_width+1) ~= 0
  error('width should be divisible by tile_width+1')
end

write_to_disk = 1;  % Write data to disk using writeG2oDataset2D
folder_name = 'width8_9robots'; %horzcat('width',num2str(width));
dir_name = horzcat('../data/blocks_world/',folder_name,'/'); % Directory in which files are written, files are named as 1.g2o, 2.g2o ...
mkdir(dir_name);

nrNodes = height*width*depth
graph_full = zeros(3*nrNodes, 44); %For saving the full graph to disk (rows = nr of factors)


%% create ground truth - 
% Graph (separators in paranthesis)
% 0 ---  1 ---  2 ---- ..... width-1
% |      |      |              |
% width+1  ....             2*width-1
% |
%           ...
%      
% |      |      |              |
%           ....         (height-1)*(width-1)

disp('setting ground truth and initial values')

ground_truth = gtsam.Values;
initial = gtsam.Values;
for k = 0:depth-1
    for i = 0:height-1
        for j = 0:width-1      
            node_id = height*width*k + width*i + j;
            ground_truth.insert(node_id, Pose3(Rot3, Point3(i*2,j*2,k*2)));

            if(j~=0)
                prev_node_id = height*width*k + width*i + j-1;
            elseif(i~=0)
                prev_node_id = height*width*k + width*(i-1) + j;
            elseif(k~=0)
                prev_node_id = height*width*(k-1) + width*(i) + j;
            end
            
            if(i== 0 && j == 0 && k == 0)
                initial.insert(node_id, ground_truth.atPose3(node_id).compose(Pose3(Rot3.Ypr(sigma_yaw*randn(1), sigma_pitch*randn(1), sigma_roll*randn(1)), Point3(sigma_x*randn(1), sigma_y*randn(1), sigma_z*randn(1)))));
            else
                delta = ground_truth.atPose3(prev_node_id).between( ground_truth.atPose3(node_id));
                noisy_delta = delta.compose(Pose3(Rot3.Ypr(sigma_yaw*randn(1), sigma_pitch*randn(1), sigma_roll*randn(1)), Point3(sigma_x*randn(1), sigma_y*randn(1), sigma_z*randn(1))));
                initial.insert(node_id, initial.atPose3(prev_node_id).compose(noisy_delta));
            end
        end
    end
end

disp('creating factor graph')
%% create relative pose measurements
fg = NonlinearFactorGraph;
fg.add(PriorFactorPose3(0, p0, priorNoise)); 
factorId = 1;
% Connect each edge with its neighbors in i-1, j-1, k-1
for k = 0:depth-1
    for i = 0:height-1
        for j = 0:width-1
            
            % delta between (i-1,j,k) and (i,j,k)
            if i~=0
                k1 = k*width*height + (i-1)*(width) +j;
                k2 = k*width*height + i*width + j;
                delta = ground_truth.atPose3(k1).between( ground_truth.atPose3(k2));
                noisy_delta = delta.compose(Pose3(Rot3.Ypr(sigma_yaw*randn(1), sigma_pitch*randn(1), sigma_roll*randn(1)), Point3(sigma_x*randn(1), sigma_y*randn(1), sigma_z*randn(1))));
                Dxyth = [noisy_delta.x  noisy_delta.y  noisy_delta.z noisy_delta.rotation.yaw noisy_delta.rotation.pitch noisy_delta.rotation.roll];
                fg.add(BetweenFactorPose3(k1, k2, noisy_delta, betweenNoise));
                graph_full(factorId,:) = [k1 k2 Dxyth, upperTriangularNoise];
                factorId = factorId+1;
            end
            
            % delta between (i,j-1,k) and (i,j,k)
            if j~=0
                k1 = k*width*height + i*width + j-1;
                k2 = k*width*height + i*width + j;
                delta = ground_truth.atPose3(k1).between( ground_truth.atPose3(k2));
                noisy_delta = delta.compose(Pose3(Rot3.Ypr(sigma_yaw*randn(1), sigma_pitch*randn(1), sigma_roll*randn(1)), Point3(sigma_x*randn(1), sigma_y*randn(1), sigma_z*randn(1))));
                Dxyth = [noisy_delta.x  noisy_delta.y  noisy_delta.z noisy_delta.rotation.yaw noisy_delta.rotation.pitch noisy_delta.rotation.roll];
                fg.add(BetweenFactorPose3(k1, k2, noisy_delta, betweenNoise));                
                graph_full(factorId,:) =  [k1 k2 Dxyth, upperTriangularNoise];
                factorId = factorId+1;                
            end
            
            
            % delta between (i,j,k-1) and (i,j,k)
            if k~=0
                k1 = (k-1)*width*height + i*width + j;
                k2 = k*width*height + i*width + j;
                delta = ground_truth.atPose3(k1).between(ground_truth.atPose3(k2));
                noisy_delta = delta.compose(Pose3(Rot3.Ypr(sigma_yaw*randn(1), sigma_pitch*randn(1), sigma_roll*randn(1)), Point3(sigma_x*randn(1), sigma_y*randn(1), sigma_z*randn(1))));
                Dxyth = [noisy_delta.x  noisy_delta.y  noisy_delta.z noisy_delta.rotation.yaw noisy_delta.rotation.pitch noisy_delta.rotation.roll];
                fg.add(BetweenFactorPose3(k1, k2, noisy_delta, betweenNoise));                
                graph_full(factorId,:) = [k1 k2 Dxyth, upperTriangularNoise];
                factorId = factorId+1;               
            end
        end
    end
end
       
nrBetween = factorId - 1;
graph_full = graph_full(1:nrBetween,:); % cut the structure to the actual number of factors

disp('writing to disk');
% write to disk
if(write_to_disk)
    file_path = [dir_name,'fullGraph.g2o'];
    writeG2o(fg, initial, file_path);
end
    
%% Centralized optimization
optimizer = GaussNewtonOptimizer(fg, initial);
result = optimizer.optimize;
marginals = Marginals(fg, result);

% % Plotting (does not plot out of order factors)
% figure(1); cla; hold on;
% plot3DSubGraph(fg, result, '-', 'r');
% plot3DSubGraph(fg, ground_truth, '-','g');
% plot3DSubGraph(fg, initial, '-', 'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISTRIBUTED OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_tiles_width = width/(tile_width+1);
num_tiles_height = height/(tile_height+1);
num_tiles_depth = depth/(tile_depth+1);
subgraph_symbols = ['a':'z','A':'Z'];
robots= [1:num_subgraphs];


%% Create subgraphs
for subgraph_id = 1:num_subgraphs
    
    % sub-initial
    subgraph_symbol = subgraph_symbols(subgraph_id)
   %robot_id = robots(subgraph_id);
    initial_subgraph_curr = gtsam.Values;    
    
    loc = (subgraph_id -1);
    subgraph_k = tile_depth*floor(loc/(num_tiles_width*num_tiles_height)); 
    subgraph_ij = mod(loc, (num_tiles_width*num_tiles_height));    
    subgraph_i = tile_height*floor(subgraph_ij/num_tiles_width);
    subgraph_j = tile_width*mod(subgraph_ij, num_tiles_width);
    
    % +1 to skip the boundary edge
    if(subgraph_i ~= 0)
        subgraph_i = 2*subgraph_i;
    end
    
    if(subgraph_j ~= 0)
        subgraph_j = 2*subgraph_j;
    end
    
    if(subgraph_k ~= 0)
        subgraph_k = 2*subgraph_k;
    end
           
    for k = subgraph_k:subgraph_k+tile_depth
        for i = subgraph_i:subgraph_i+tile_height
            for j = subgraph_j:subgraph_j+tile_width
                pose_id = k*width*height + i*width+j;
                pose_symbol = symbol(subgraph_symbol, pose_id);
                pose = initial.atPose3(pose_id);
                initial_subgraph_curr.insert(pose_symbol, pose);                
            end
        end
    end
    initials(subgraph_id) = initial_subgraph_curr;
    
    % sub-graph
    subgraph_curr = gtsam.NonlinearFactorGraph;
    
    
    if subgraph_id == 1
        subgraph_curr.add(PriorFactorPose3(0, p0, priorNoise));
    end
    
    for k = subgraph_k:subgraph_k+tile_depth
        for i = subgraph_i:subgraph_i+tile_height
            for j = subgraph_j:subgraph_j+tile_width
                if i~=subgraph_i
                    k1 = k*width*height + (i-1)*(width) +j;
                    k2 = k*width*height + (i)*(width) +j;
                    index  = find(graph_full(:,1) == k1 & graph_full(:,2) == k2);
                    symbol_1 = symbol(subgraph_symbol, graph_full(index,1));
                    symbol_2 = symbol(subgraph_symbol, graph_full(index,2));
                    noisy_delta = Pose3(Rot3.Ypr(graph_full(index,6), graph_full(index,7), graph_full(index,8)), Point3(graph_full(index,3), graph_full(index,4), graph_full(index,5)));
                    subgraph_curr.add(BetweenFactorPose3(symbol_1, symbol_2, noisy_delta, betweenNoise));                    
                end
                
                if j~=subgraph_j
                    k1 = k*width*height + (i)*(width) +j-1;
                    k2 = k*width*height + (i)*(width) +j;
                    index  = find(graph_full(:,1) == k1 & graph_full(:,2) == k2);
                    symbol_1 = symbol(subgraph_symbol, graph_full(index,1));
                    symbol_2 = symbol(subgraph_symbol, graph_full(index,2));
                    noisy_delta = Pose3(Rot3.Ypr(graph_full(index,6), graph_full(index,7), graph_full(index,8)), Point3(graph_full(index,3), graph_full(index,4), graph_full(index,5)));
                    subgraph_curr.add(BetweenFactorPose3(symbol_1, symbol_2, noisy_delta, betweenNoise));
                    
                end
                
                if k~=subgraph_k
                    k1 = (k-1)*width*height + (i)*(width) +j;
                    k2 = k*width*height + (i)*(width) +j;
                    index  = find(graph_full(:,1) == k1 & graph_full(:,2) == k2);                    
                    symbol_1 = symbol(subgraph_symbol, graph_full(index,1));
                    symbol_2 = symbol(subgraph_symbol, graph_full(index,2));
                    noisy_delta = Pose3(Rot3.Ypr(graph_full(index,6), graph_full(index,7), graph_full(index,8)), Point3(graph_full(index,3), graph_full(index,4), graph_full(index,5)));
                    subgraph_curr.add(BetweenFactorPose3(symbol_1, symbol_2, noisy_delta, betweenNoise));                    
                end
                
                % Separators
                
                % Previous k subgraph
                if k == subgraph_k && k~=0
                    k1 = (k-1)*width*height + (i)*(width) +j;
                    k2 = k*width*height + (i)*(width) +j;
                    index  = find(graph_full(:,1) == k1 & graph_full(:,2) == k2);                    
                    previous_subgraph = subgraph_id - num_tiles_width*num_tiles_height;
                    prev_subgraph_symbol = subgraph_symbols(previous_subgraph);
                    symbol_1 = symbol(prev_subgraph_symbol, graph_full(index,1));
                    symbol_2 = symbol(subgraph_symbol, graph_full(index,2));                    
                    noisy_delta = Pose3(Rot3.Ypr(graph_full(index,6), graph_full(index,7), graph_full(index,8)), Point3(graph_full(index,3), graph_full(index,4), graph_full(index,5)));
                    subgraph_curr.add(BetweenFactorPose3(symbol_1, symbol_2, noisy_delta, betweenNoise));                    
                end
                
                
                % Previous j subgraph
                if j == subgraph_j && j~=0
                    k1 = k*width*height + (i)*(width) +j-1;
                    k2 = k*width*height + (i)*(width) +j;                    
                    index  = find(graph_full(:,1) == k1 & graph_full(:,2) == k2);                    
                    previous_subgraph = subgraph_id -1;
                    prev_subgraph_symbol = subgraph_symbols(previous_subgraph);
                    symbol_1 = symbol(prev_subgraph_symbol, graph_full(index,1));
                    symbol_2 = symbol(subgraph_symbol, graph_full(index,2));                    
                    noisy_delta = Pose3(Rot3.Ypr(graph_full(index,6), graph_full(index,7), graph_full(index,8)), Point3(graph_full(index,3), graph_full(index,4), graph_full(index,5)));
                    subgraph_curr.add(BetweenFactorPose3(symbol_1, symbol_2, noisy_delta, betweenNoise));                    
                end
                
                % Previous i subgraph
                if i == subgraph_i && i~=0
                    k1 = k*width*height + (i-1)*(width) +j;
                    k2 = k*width*height + (i)*(width) +j;
                    index  = find(graph_full(:,1) == k1 & graph_full(:,2) == k2);
                    previous_subgraph = subgraph_id - num_tiles_width;
                    prev_subgraph_symbol = subgraph_symbols(previous_subgraph);
                    symbol_1 = symbol(prev_subgraph_symbol, graph_full(index,1));
                    symbol_2 = symbol(subgraph_symbol, graph_full(index,2));                    
                    noisy_delta = Pose3(Rot3.Ypr(graph_full(index,6), graph_full(index,7), graph_full(index,8)), Point3(graph_full(index,3), graph_full(index,4), graph_full(index,5)));
                    subgraph_curr.add(BetweenFactorPose3(symbol_1, symbol_2, noisy_delta, betweenNoise));                    
                end
                
                % Next k subgraph
                if k == subgraph_k + tile_depth && (k+1) < depth
                    k1 = (k)*width*height + (i)*(width) +j;
                    k2 = (k+1)*width*height + (i)*(width) +j;
                    index  = find(graph_full(:,1) == k1 & graph_full(:,2) == k2);                    
                    next_subgraph = subgraph_id + num_tiles_width*num_tiles_height;
                    next_subgraph_symbol = subgraph_symbols(next_subgraph);
                    symbol_1 = symbol(subgraph_symbol, graph_full(index,1));
                    symbol_2 = symbol(next_subgraph_symbol, graph_full(index,2));                    
                    noisy_delta = Pose3(Rot3.Ypr(graph_full(index,6), graph_full(index,7), graph_full(index,8)), Point3(graph_full(index,3), graph_full(index,4), graph_full(index,5)));
                    subgraph_curr.add(BetweenFactorPose3(symbol_1, symbol_2, noisy_delta, betweenNoise));                    
                end
                
                
                % Next j subgraph
                if j == subgraph_j + tile_width && (j+1) < width
                    k1 = k*width*height + (i)*(width) +j;
                    k2 = k*width*height + (i)*(width) +j+1;                    
                    index  = find(graph_full(:,1) == k1 & graph_full(:,2) == k2);                    
                    next_subgraph = subgraph_id +1;
                    next_subgraph_symbol = subgraph_symbols(next_subgraph);
                    symbol_1 = symbol(subgraph_symbol, graph_full(index,1));
                    symbol_2 = symbol(next_subgraph_symbol, graph_full(index,2));                    
                    noisy_delta = Pose3(Rot3.Ypr(graph_full(index,6), graph_full(index,7), graph_full(index,8)), Point3(graph_full(index,3), graph_full(index,4), graph_full(index,5)));
                    subgraph_curr.add(BetweenFactorPose3(symbol_1, symbol_2, noisy_delta, betweenNoise));                    
                end
                
                % Next i subgraph
                if i == subgraph_i + tile_height && (i+1) < height
                    k1 = k*width*height + (i)*(width) +j;
                    k2 = k*width*height + (i+1)*(width) +j;
                    index  = find(graph_full(:,1) == k1 & graph_full(:,2) == k2);
                    next_subgraph = subgraph_id + num_tiles_width;
                    next_subgraph_symbol = subgraph_symbols(next_subgraph);
                    symbol_1 = symbol(subgraph_symbol, graph_full(index,1));
                    symbol_2 = symbol(next_subgraph_symbol, graph_full(index,2));                    
                    noisy_delta = Pose3(Rot3.Ypr(graph_full(index,6), graph_full(index,7), graph_full(index,8)), Point3(graph_full(index,3), graph_full(index,4), graph_full(index,5)));
                    subgraph_curr.add(BetweenFactorPose3(symbol_1, symbol_2, noisy_delta, betweenNoise));                    
                end
                
            end
        end
    end
    
    
    subgraphs(subgraph_id) = subgraph_curr;        

    % write to disk
    if(write_to_disk)
         file_path = [dir_name,num2str(subgraph_id-1),'.g2o'];
        writeG2o(subgraph_curr, initial_subgraph_curr, file_path);
    end
    
end



%% Error centralized
fg_error = fg.error(result);
fprintf('Final error centralized fg: %g\n', fg_error)
    
%% Solve using Distributed Optimization
%[initial_x, y_k] = ADMM_parallel(subgraphs, initials, separators, 'linearADMM', fg, fg_error); % fg is only passed for debug
%errADMM = computeCentralizedErrorFromSubgraphs(fg, initial_x);
%fprintf('Final error centralized graph w.r.t ADMM result:%g\n', errADMM)
