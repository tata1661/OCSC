clear
dbstop if error
addpath('basic_tool'); 
addpath('OCSC');
addpath('mtimesx');%**
%% set para
K = [100];
psf_s=11;                                                                                                      
psf_radius = floor( psf_s/2 );
precS = 1;
use_gpu = 1;
data = 'city_10';
data = 'fruit_10';
%% load data
load (sprintf('datasets/%s/train/train_lcne.mat',data)) %%% 
padB = padarray(b, [psf_radius, psf_radius, 0], 0, 'both');
PARA= auto_para(K,psf_s,b,'no',1e-3,precS,use_gpu);
%% run
t1 = tic;
[ d,d_hat]  = alt_min_online(padB,PARA,[],b);
tt = toc(t1);
%% save
repo_name = 'result';
repo_path = sprintf('%s/%s',repo_name,data);
save_name = sprintf('K%d_psf%d',K,psf_s);
save_me = sprintf('%s/record_%s.mat',repo_path,save_name);
save(save_me,'d_hat','d','tt');
fprintf('Done sparse coding learning! --> Time %2.2f sec.\n\n', tt)