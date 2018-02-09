function [d_curr, d_hat] = alt_min_online(Mtb,para,init,b)
%% Initialize variables
if ~isempty(init)
    if isfield(init, 'd')     
        d_small = init.d;
        d = dsmall2d(d_small,para);
        d_hat = fft2(d); 
    elseif isfield(init, 'd_hat')
        d_hat=init.d_hat;
    end
    if isfield(init,'A_h')
        A_h=init.A_h;
    else
        A_h = [];
    end
    if isfield(init,'B_h')
        B_h=init.B_h;
    else
        B_h=[];
    end    
else
    A_h = [];
    B_h=[];
    d_small = init_dic(para);
    %d_small = randn(para.size_k);
    if (para.precS ==1)
        d_small = single(d_small);
    end
    if (para.gpu ==1)
        d_small = gpuArray(d_small);
    end   
    d = dsmall2d(d_small,para);
    d_hat = fft2(d); 
end
if (para.precS ==1)
    Mtb = single(Mtb);
end
if (para.gpu ==1)
    Mtb = gpuArray(Mtb);
end  
b_hat = fft2(Mtb);
%%
s=[];
y=[];
scale = 1;
for s_i=1:para.N
    temp_b = b(:,:,s_i);
    temp_b_hat = b_hat(:,:,s_i);
    %% 1.pre-process Z
    t_Z = tic;%~~~~!!!
    [stat_Z] = precompute_H_hat_Z(d_hat, para);    
    %% 2.update Z
    [z_si,z_hat_si] = updateZ_ocsc(temp_b_hat,para,d_hat,stat_Z);
    timeZ = toc(t_Z);
    objZ = objective_online(z_hat_si,d_hat, temp_b_hat,para );
   if strcmp( para.verbose, 'all')
       if (mod(s_i,scale)==0)
            [ps] = eval_psnr(d_hat, z_hat_si,temp_b,para); 
            fprintf('Z: no.img: %d, obj: %2.2f, psnr: %2.2f\n', s_i,objZ,ps)
        end 
    end
    clear stat_Z           
    %% 1.pre-process D
    if isempty(A_h)
        init_or_not = 1;
    else
        init_or_not = 2;
    end
    t_D =tic;
    if para.gpu ==1
        [A_h,B_h] = hist_ocsc_gpu(temp_b_hat,z_hat_si, para, A_h,B_h,init_or_not);
    else
        [A_h,B_h] = hist_ocsc_cpu(temp_b_hat,z_hat_si, para, A_h,B_h,init_or_not);
    end
    %% 2.update D
    [d,d_hat,s,y] = updateD_ocsc(para,A_h,B_h,s,y,d_hat);    
    timeD =toc(t_D);
    d_curr = d2dsmall(d,para);
    if strcmp( para.verbose, 'all')
        if para.gpu==1
            d_show = gather(d_curr);
            show_dic(d_show,para,0,0);
        else
            show_dic(d_curr,para,0,0); 
        end
    end
end
end