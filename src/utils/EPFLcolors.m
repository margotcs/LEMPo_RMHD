function EPFLcolors = EPFLcolors()
% EPFLcolors - returns a struct mapping EPFL color names to RGB triplets
% Usage: colors = EPFLcolors(); plot(x, y, 'Color', colors.canard);


colortab = {'groseille', 'canard', 'taupe', 'rose', 'montRose', 'zinzolin', ...
            'vertDeau', 'leman', 'rougeSuisse', 'chartreuse', 'ardoise', 'marron'};

sz = numel(colortab);
colormatHex = cell(sz, 1);
colormatRGB = zeros(sz, 3);

for ii = 1:sz
    switch colortab{ii}
        case 'groseille'
            colormatHex{ii} = '#B51F1F';
        case 'canard'
            colormatHex{ii} = '#007480';
        case 'taupe'
            colormatHex{ii} = '#413D3A';
        case 'rose'
            colormatHex{ii} = '#ED6E9C';
        case 'montRose'
            colormatHex{ii} = '#F39869';
        case 'zinzolin'
            colormatHex{ii} = '#5C2483';
        case 'vertDeau'
            colormatHex{ii} = '#C2DD80';
        case 'leman'
            colormatHex{ii} = '#00A797';
        case 'rougeSuisse'
            colormatHex{ii} = '#FF0000';
        case 'chartreuse'
            colormatHex{ii} = '#C8D300';
        case 'ardoise'
            colormatHex{ii} = '#453A4C';
        case 'marron'
            colormatHex{ii} = '#5B3428';
        otherwise
            colormatHex{ii} = '#000000'; % fallback: black
    end
end

for i = 1:sz
    str = colormatHex{i};
    RGB = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    colormatRGB(i, :) = RGB;
end

for i = 1:sz
    EPFLcolors.(colortab{i}) = colormatRGB(i,:);
end

end