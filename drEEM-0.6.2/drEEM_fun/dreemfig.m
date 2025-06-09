function [fighandleout] = dreemfig(fighandlein,nrow,ncol)
% Dreemfig.m: (c) Urban J. Wuensch (2020)
% This function is part of the drEEM toolbox.
% It contains no documntation since it is a function not intended for use
% by the end user directly, but is only used by other functions within
% drEEM. dreemfig ensures that all figures produced in the toolbox are
% formatted equally.

FontSize=8;
FontName='Arial';
if ~exist('fighandlein','var')
    fighandlein=[];
end

if isempty(fighandlein) % Case: Make new figure
    fighandleout=figure('InvertHardcopy','off','Color',[1 1 1]);
    try
        set(fighandleout,...
            'defaultUicontrolFontName',FontName,...
            'defaultUitableFontName',FontName,...
            'defaultAxesFontName',FontName,...
            'defaultTextFontName',FontName,...
            'defaultUipanelFontName',FontName,...
            'defaultUicontrolFontSize',FontSize,...
            'defaultUitableFontSize',FontSize,...
            'defaultAxesFontSize',FontSize,...
            'defaultTextFontSize',FontSize,...
            'defaultUipanelFontSize',FontSize);
    catch
    end
%     if exist('nrow','var')&&exist('ncol','var')
%         pos=sizefigure(nrow,ncol);
%         % Check if position arguments extend beyond screen boundary
%         if any(pos>1)
%             pos(pos>1)=1;
%         end
%         
%         set(fighandleout,'units','normalized','position',pos)
%     end
elseif isgraphics(fighandlein) % Case: Make sure formatting is how we want it to be.
    ax = (findobj(fighandlein, 'type', 'axes'));
    for n=numel(ax):-1:1
        set(ax(n),'TickDir','out');
        set(ax(n),'FontSize',FontSize,'FontName',FontName);
        set(ax(n),'LineWidth',0.5);
    end
    fighandleout=fighandlein;
end

end
% 
% function [positionArg] = sizefigure(nrow,ncol,wthr)
% 
% 
% 
% scrsz = get(0,'ScreenSize');
% scrzyorg=scrsz;
% scrsz(4) = scrsz(4)*0.7;
% wth=scrsz(3)./scrsz(4);
% 
% 
% ppc = scrsz(3)./ncol;
% ppr = scrsz(4)./nrow;
% 
% 
% 
% wrel = (scrsz(3)/(ppc/ppr))/scrzyorg(3);
% hrel = scrsz(4)/scrzyorg(4);
% lrel = (1-wrel)./2;
% brel = (1-hrel)./2;
% 
% positionArg = [lrel brel wrel hrel];
% 
% % figure('units','normalized','position',[lrel brel wrel hrel])
% % hf=figure;
% % set(hf,'Position',[scrsz(1) scrsz(2) scrsz(3)/(ppc/ppr) scrsz(4)]);
% %
% % for n=1:nrow*ncol
% %     subplot(nrow,ncol,n)
% % end
% 
% end
% 
