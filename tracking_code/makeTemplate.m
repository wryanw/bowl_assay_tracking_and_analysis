function makeTemplate(templCell,templPath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%
tmplCt = numel(templCell);
fly_length = zeros(tmplCt,1);
for iterT = 1:tmplCt
    tmplStr = templCell{iterT};
    hPos = tmplStr.hPos;
    tPos = tmplStr.tPos;
    fly_length(iterT) = sqrt(sum((hPos-tPos).^2));
end
[fly_length,maxRef] = max(fly_length);
tmplStr = templCell{maxRef};
ndxr_struct = tmplStr.ndxr_struct;
spoke_length = ndxr_struct.spoke_length;
spoke_leg = ndxr_struct.spoke_leg;
spoke_count = ndxr_struct.spoke_count;
templ_ndxr_mastr = (ndxr_struct.templ_mastr);
spoke_ndxr_mastr = (ndxr_struct.spoke_mastr);

whole_XData = [1 fly_length+im_leg*2];
whole_YData = whole_XData;
whole_dests = [im_leg,...
    median(whole_YData);
    fly_length+im_leg,...
    median(whole_YData)];

frmCell = cell(tmplCt,1);
for iterT = 1:tmplCt
    
    tmplStr = templCell{iterT};
    ex_pts = [tmplStr.hPos;tmplStr.tPos];
    Q_tform = cp2tform(ex_pts,whole_dests,'nonreflective similarity');
    frm_bot = tmplStr.frm_bot;
    frmTrans = imtransform(frm_bot,Q_tform,'XData',whole_XData,...
        'YData',whole_YData,'FillValues',0.01);
    frmCell{iterT} = frmTrans;
    %         imshow(frmTrans)
    %         uiwait(gcf)
end
frm_bot = (mean(cat(3,frmCell{:}),3));
%     imshow(frm_bot)
%     size(frm_bot)
%     hold on
%     plot(whole_dests(:,1),whole_dests(:,2))
%     [sqrt(sum((whole_dests(1,:)-whole_dests(2,:)).^2)) fly_length]

spoke_span = spoke_leg*2+1;
radians_per_spoke = 2*pi/spoke_count;
fly_theta = pi;%template default
layer_ref = round(fly_theta/radians_per_spoke);
layer_ref(layer_ref > spoke_count) = layer_ref-spoke_count;
layer_ref(layer_ref <= 0 ) = layer_ref+spoke_count;

frm_pad = padarray(frm_bot,[im_leg im_leg]);
fly_pos = mean(whole_dests);
neg_dim = round(fly_pos);
pos_dim = round(fly_pos+im_leg*2);
templUniversal = frm_pad(neg_dim(2):pos_dim(2),neg_dim(1):pos_dim(1));

templ_ndxr_layer = templ_ndxr_mastr(spoke_ndxr_mastr(:,layer_ref),:);
templ_singl = templUniversal(templ_ndxr_layer);
templIm = reshape(templ_singl,spoke_length,spoke_span);
%     imshow(templIm)

%%

save(templPath,'frm_bot','whole_dests')

end

