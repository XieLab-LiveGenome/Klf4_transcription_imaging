function masks = draw_nuclei_manual(img, masks, tool)
%DRAW_NUCLEI_MANUAL  Hand-trace nucleus ROIs on an image.
%
%   masks = draw_nuclei_manual(img)                  % start from scratch
%   masks = draw_nuclei_manual(img, masks)           % add to existing ROIs
%   masks = draw_nuclei_manual(img, masks, 'polygon')
%
%   img   : 2-D image (any numeric class)
%   masks : 1xN cell array of logical masks, each size(img)
%   tool  : 'freehand' (default) | 'polygon'
%
%   Draw ONE nucleus at a time.
%     freehand : click-drag-release
%     polygon  : click vertices, double-click to close
%   Press ESC (instead of drawing) to stop; you then get a
%   Done / Undo last / Keep drawing prompt.
%
%   Output is a cell row of logical masks, drop-in compatible with the
%   output of cellpose_seg{1}.
%
%   Requires R2018b+ (drawfreehand / drawpolygon).

if nargin < 2 || isempty(masks), masks = {}; end
if nargin < 3 || isempty(tool),  tool  = 'freehand'; end
if ~exist('drawfreehand', 'file')
    error('draw_nuclei_manual:oldMATLAB', ...
        'drawfreehand/drawpolygon require R2018b or newer.');
end

masks    = masks(:).';                            % force cell row
disp_img = imadjust(mat2gray(double(img)));
MIN_PX   = 100;                                   % discard accidental micro-ROIs

f   = figure('Name','Manual nucleus ROIs','NumberTitle','off', ...
             'Position',[80 80 950 950]);
ax  = axes('Parent', f);
him = imshow(disp_img, 'Parent', ax, 'InitialMagnification','fit');
hold(ax, 'on');
for k = 1:numel(masks), local_outline(ax, masks{k}, k, [0 1 0]); end

while true
    title(ax, sprintf('%d nucleus ROI(s).  Draw the next one (%s) — ESC when finished.', ...
        numel(masks), tool), 'FontSize', 12);
    drawnow;

    switch lower(tool)
        case 'polygon'
            h = drawpolygon(ax, 'Color','y', 'LineWidth',1.5);
        otherwise
            h = drawfreehand(ax, 'Color','y', 'LineWidth',1.5, 'Closed',true);
    end

    % ESC while drawing -> ROI returns with empty/degenerate Position
    if ~isvalid(h) || isempty(h.Position) || size(h.Position,1) < 3
        if isvalid(h), delete(h); end
        choice = questdlg(sprintf('Finish with %d ROI(s)?', numel(masks)), ...
                          'Manual ROIs', 'Done', 'Undo last', 'Keep drawing', 'Done');
        switch choice
            case 'Done'
                break
            case 'Undo last'
                if ~isempty(masks)
                    masks(end) = [];
                    cla(ax);
                    him = imshow(disp_img, 'Parent', ax, 'InitialMagnification','fit');
                    hold(ax, 'on');
                    for k = 1:numel(masks), local_outline(ax, masks{k}, k, [0 1 0]); end
                end
            otherwise   % 'Keep drawing' or dialog closed
        end
        continue
    end

    bw = createMask(h, him);
    delete(h);

    if nnz(bw) < MIN_PX
        warning('draw_nuclei_manual:tinyROI', ...
            'ROI had only %d px (< %d) — discarded.', nnz(bw), MIN_PX);
        continue
    end

    masks{end+1} = bw;                                     %#ok<AGROW>
    local_outline(ax, bw, numel(masks), [0 1 0]);
end

title(ax, sprintf('Final: %d nucleus ROI(s)', numel(masks)), 'FontSize', 12);
masks = masks(:).';
end

% -------------------------------------------------------------------------
function local_outline(ax, bw, idx, col)
B = bwboundaries(bw);
for b = 1:numel(B)
    plot(ax, B{b}(:,2), B{b}(:,1), '-', 'Color', col, 'LineWidth', 1.5);
end
st = regionprops(bw, 'Centroid');
if ~isempty(st)
    text(ax, st(1).Centroid(1), st(1).Centroid(2), num2str(idx), ...
        'Color','w', 'FontSize',14, 'FontWeight','bold', ...
        'HorizontalAlignment','center');
end
end
