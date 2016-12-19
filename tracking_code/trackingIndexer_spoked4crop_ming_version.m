function [ndxr_struct] = trackingIndexer_spoked4crop_ming_version(src_leg,trk_hort_ref,pos_ops,theta_ops,spoke_count)
% trackingIndexer_spoked_LR Generates linear indexing variables for
% tracking flies in pez videos
%   Inputs required: radius of resample (src_leg) which will be reduced by
%   one-third, removing the inner-most third, and head versus tail
%   reference (trk_hort_ref).  Outputs include variables needed for tracking
%   and visualization, condensed in a structure class variable.  IMPORTANT:
%   THE TOP HALF OF THE REINDEXED OUTPUT IS A SYMMETRIC MIRROR OF THE BOTTOM.

%%
if nargin == 0 || isempty(mfilename)
    src_leg = 31;
    trk_hort_ref = 1;
    theta_ops = [fliplr(thetaLeg*(-1)) 0 thetaLeg];
    heading_factor = 3;
    spoke_count = 360/heading_factor;
    theta_ops = round([theta_ops spoke_count/2 spoke_count/4 -spoke_count/4]);
    pos_ops = [fliplr(posLeg*(-1)) 0 posLeg];
end
im_leg = max(abs(pos_ops))+src_leg;
im_dim = im_leg*2+1;
col_indexer = ((1:im_dim)-1).*im_dim;
% spoke_span_factor = 1/2;% fraction of coverage on either side of center
% spoke_span_factor = 2/3;% fraction of coverage on either side of center
spoke_span_factor = 1;% fraction of coverage on either side of center
pos_count = numel(pos_ops);
theta_count = numel(theta_ops);

spoke_leg = round(spoke_count*spoke_span_factor/2);
spoke_span = spoke_leg*2+1;

% spoke_span_ops = [spoke_count/2-spoke_leg+1:spoke_count/2+spoke_leg+1;
%     spoke_count-spoke_leg+1:spoke_count,1:spoke_leg+1];%points toward fly center
spoke_span_ops = [spoke_count-spoke_leg+1:spoke_count,1:spoke_leg+1;
    spoke_count/2-spoke_leg+1:spoke_count/2+spoke_leg+1];%points away from center

y_ops = repmat(pos_ops',[1 pos_count theta_count]);
x_ops = repmat(pos_ops,[pos_count 1 theta_count]);
theta_findr = repmat(theta_ops,pos_count^2,1);
mastr_findr = [x_ops(:) y_ops(:) theta_findr(:)];

rho_min = src_leg*(0);%0.1
rho_span = src_leg-rho_min;
spoke_length = round(rho_span);
% rho_adj = linspace(0,1,spoke_length).^(1);
% rho_adj = linspace(0,1,spoke_length).^(0.7);
% rho_adj = linspace(0,1,spoke_length).^(1.5);

resize_factor = 0.7;
indexer_bounds = [src_leg src_leg]*2-4;
coverage = 2;
pixel_last_of_small = indexer_bounds(1)/2;%mid-point to end of smallest flies
head_half_small = round(resize_factor*pixel_last_of_small);
half_large = indexer_bounds(2)/2;%mid-point to end of largest flies, defines image block dim!!!!
head_half_large = round(resize_factor*half_large);
pixel_last_of_large = head_half_large;%adjusted for indexing vector
pts_of_origin_count = spoke_length;
pts_including_boundaries_count = pts_of_origin_count*2+1;
factor = 1+(1/pts_including_boundaries_count)*coverage;
pts_max_last = round(log(pixel_last_of_large)/log(factor));
pts_distr = factor.^(1:pts_max_last);%redistribute points based on factor
pts_min_last = find(max(pts_distr(pts_distr < head_half_small))== pts_distr);
if pts_including_boundaries_count > pts_min_last-1
    pts_including_boundaries_count = pts_min_last-1;
end
pts_min_first = pts_min_last-pts_including_boundaries_count;
pix_avg_ind_mids = fliplr(pts_distr(pts_min_first:pts_max_last));
if max(round(pix_avg_ind_mids)) > src_leg
    midRef = find(round(pix_avg_ind_mids) > src_leg,1,'last')+1;
else
    midRef = 2;
end
ndx_rhos = round(fliplr(pix_avg_ind_mids(midRef:2:end-1)));
ndx_rhos = unique(ndx_rhos);
% ndx_rhos = round(rho_adj.*rho_span+rho_min)

spoke_length = numel(ndx_rhos);

deg2rad = 2*pi/spoke_count;
rhos_array = repmat(ndx_rhos(:),1,spoke_count);
thetas_array = repmat((1:spoke_count)-1,spoke_length,1).*deg2rad;
rot_ndx_x = round(rhos_array.*cos(thetas_array));
rot_ndx_y = -round(rhos_array.*sin(thetas_array));

sampl_ndx_r1 = repmat(pos_ops,1,pos_count);
sampl_ndx_c1 = repmat(pos_ops,pos_count,1);
sampl_ndx_r2 = repmat(sampl_ndx_r1(:),[1 spoke_length spoke_count]);
sampl_ndx_r3 = permute(sampl_ndx_r2,[2 3 1]);
sampl_ndx_c2 = repmat(sampl_ndx_c1(:),[1 spoke_length spoke_count]);
sampl_ndx_c3 = permute(sampl_ndx_c2,[2 3 1]);
sampl_ndx_r4 = repmat(rot_ndx_y,[1 1 pos_count^2])+sampl_ndx_r3+im_leg+1;
sampl_ndx_c4 = repmat(rot_ndx_x,[1 1 pos_count^2])+sampl_ndx_c3+im_leg+1;
sampl_ndx_m1 = sampl_ndx_r4+col_indexer(sampl_ndx_c4);
sampl_ndx_m2 = reshape(sampl_ndx_m1,spoke_length*spoke_count,pos_count^2);
sampl_ndx_m3 = zeros(spoke_length*spoke_count,pos_count^2,theta_count);
for j = 1:theta_count
    shift_val = theta_ops(j)*spoke_length;
    sampl_ndx_m3(:,:,j) = circshift(sampl_ndx_m2,-shift_val);
    if theta_ops(j) == 0
        templ_ndxr_mastr = sampl_ndx_m3(:,median(1:pos_count^2),j);
    end
end
sampl_ndxr_mastr = reshape(sampl_ndx_m3,spoke_length*spoke_count,pos_count^2*theta_count);
templ_rendxr = repmat((1:spoke_length*spoke_span)',1,pos_count^2*theta_count);

spoke_ndxr1 = spoke_span_ops(trk_hort_ref,:);
spoke_ndxr2a = repmat(spoke_ndxr1(:),1,spoke_count);
spoke_ndxr2b = repmat((1:spoke_count)-1,spoke_span,1);
spoke_ndxr3 = spoke_ndxr2a+spoke_ndxr2b;
spoke_adj = zeros(spoke_span,spoke_count);
spoke_adj(spoke_ndxr3 > spoke_count) = -spoke_count;
spoke_ndxr4 = spoke_ndxr3+spoke_adj;
spoke_ndxr5 = reshape((1:spoke_length*spoke_count),...
    spoke_length,spoke_count);
spoke_ndxr_mastr = zeros(spoke_span*spoke_length,spoke_count);
for  iter = 1:spoke_count
    spoke_ndxr6 = spoke_ndxr5(:,spoke_ndxr4(:,iter)');
    spoke_ndxr_mastr(:,iter) = spoke_ndxr6(:);
end

ndxr_struct = struct('mastr_findr',mastr_findr,'spoke_length',...
    spoke_length,'spoke_leg',spoke_leg,'spoke_count',spoke_count,'templ_mastr',templ_ndxr_mastr,...
    'sampl_mastr',sampl_ndxr_mastr,'im_leg',im_leg,...
    'spoke_mastr',spoke_ndxr_mastr,'templ_rendxr',templ_rendxr);
