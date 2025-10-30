%For E-fields

function CPplot(H,orthogonalmodes_string,x_axis,y_axis)

    % Generate random indeces to plot
    Nk = size(H,3);
    freq_idcs = round(linspace(1,Nk-5,6));
    
    % Spatial parameters
    res = length(x_axis);
    
    % Pattern details
    [maskNos, orthogonalmodes_frequencies] = DoA_getModeMaskFrequencies(orthogonalmodes_string(freq_idcs));
    
    figure
    for i = 1:6
        % Extract the pattern and normalize
        freq_idx = freq_idcs(i);
        pattern = H(:,:,freq_idx);
        pattern_n = pattern / max(abs(pattern), [], 'all');

        % Convert to dB scale
        pat_dB = 20*log10(abs(pattern_n));

        % Plot
        maskNo = maskNos(i);
        frequency = orthogonalmodes_frequencies(i);
        subplot(2,3,i)
        imagesc(x_axis,y_axis,pat_dB);
        %% Font size of plot
        fontSize = 20; % Whatever you want.
        caption = sprintf('Mask '+string(maskNo)+', f='+string(frequency/1e9)+' GHz');
        title(caption, 'FontSize', fontSize,'FontWeight','normal');
        ylabel('Y (m)', 'FontSize', fontSize)
        xlabel('Z (m)', 'FontSize', fontSize)
        zlabel('Pattern', 'FontSize', fontSize)
        %%
        xticks(linspace(x_axis(1),x_axis(end),5))
        yticks(linspace(y_axis(1),y_axis(end),5))
        colormap(jet); % Use the same colormap as in the paper
        cb = colorbar(FontSize=fontSize);
        cb.Title.String = 'dBi';
        zMax = max(pat_dB,[],'all');
        zMin = zMax-15;
        clim([zMin, zMax]) 
        cb.Ticks = [-15 -10 -5 0];
    end
    % sgtitle('|E-Fields| on Characteristic Plane, res='+string(res)) % Set subplot title
    set(gcf, 'Position', get(0, 'Screensize')); % Set figure to fullscreen
    %% Save image as pdf
    fig = gcf;
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    fileName = 'EField_CP_mag';
    print(fig,fileName,'-dpdf','-r0')
    %%
end