function masks = review_masks(img, masks)
%REVIEW_MASKS  Click to delete inaccurate nucleus masks; keep the rest.
%
%   masks = review_masks(img, masks)
%
%   Left-click inside a mask to toggle it for deletion (turns red).
%   Press ENTER / close the figure to apply. Deleted masks are removed;
%   survivors are returned renumbered. Useful for pruning a few bad
%   Cellpose masks before topping up with draw_nuclei_manual.

masks = masks(:).';
if isempty(masks), return; end

disp_img = imadjust(mat2gray(double(img)));
keep     = true(1, numel(masks));

f  = figure('Name','Review masks — click bad ones, ENTER to apply', ...
            'NumberTitle','off', 'Position',[80 80 950 950]);
ax = axes('Parent', f);
imshow(disp_img, 'Parent', ax, 'InitialMagnification','fit'); hold(ax, 'on');

hB = gobjects(1, numel(masks));
hT = gobjects(1, numel(masks));
for k = 1:numel(masks), [hB(k), hT(k)] = draw_one(ax, masks{k}, k, keep(k)); end
title(ax, 'Click a mask to toggle delete (red = delete). ENTER when done.');

set(f, 'WindowKeyPressFcn', @(~,e) strcmp(e.Key,'return') && (delete(f)==0));
set(f, 'WindowButtonDownFcn', @onclick);
uiwait(f);
masks = masks(keep);
masks = masks(:).';

    function onclick(~, ~)
        cp = get(ax, 'CurrentPoint');
        x  = round(cp(1,1)); y = round(cp(1,2));
        [H, W] = size(disp_img);
        if x < 1 || y < 1 || x > W || y > H, return; end
        for k = 1:numel(masks)
            if masks{k}(y, x)
                keep(k) = ~keep(k);
                delete(hB(k)); if isgraphics(hT(k)), delete(hT(k)); end
                col = [0 1 0]; if ~keep(k), col = [1 0 0]; end
                [hB(k), hT(k)] = draw_one(ax, masks{k}, k, keep(k), col);
                break
            end
        end
    end
end

% -------------------------------------------------------------------------
function [hb, ht] = draw_one(ax, bw, idx, keepFlag, col)
if nargin < 5, col = [0 1 0]; if ~keepFlag, col = [1 0 0]; end, end
B  = bwboundaries(bw);
hb = plot(ax, B{1}(:,2), B{1}(:,1), '-', 'Color', col, 'LineWidth', 2);
for b = 2:numel(B)
    plot(ax, B{b}(:,2), B{b}(:,1), '-', 'Color', col, 'LineWidth', 2);
end
st = regionprops(bw, 'Centroid');
ht = gobjects(1);
if ~isempty(st)
    ht = text(ax, st(1).Centroid(1), st(1).Centroid(2), num2str(idx), ...
        'Color','w', 'FontSize',13, 'FontWeight','bold', ...
        'HorizontalAlignment','center');
end
end
