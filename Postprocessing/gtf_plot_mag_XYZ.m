function gtf_plot_mag_XYZ(set_fig_nbr,TF,F_f,total_plot_width,plot_box_height,plot_lin,fmin,fmax,ymin,ymax,xlabel_str,ylabel_str,text_loc,FN_A,FS_A,ffbc,ffbr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ri = 0.99; gi = 0.8; bi = 0.8; alpha = 0.5 ; % patch color
fhv  = 'on'                                ; % HandleVisibility 'on' | 'off'

TF_X =  TF(:,1);
TF_Y =  TF(:,2);
TF_Z =  TF(:,3);


lws              = total_plot_width/300 ;

space_l          = FS_A*4 ;                % for Y-axis labeling
space_r          = FS_A/2 ;                % right margin
space_b          = FS_A*3.4375 ;           % for X-axis labeling
space_i          = FS_A*1.25 ;             % space between the plot boxes

Xplot_ypos     = space_b+(plot_box_height+space_i)*2;
Yplot_ypos     = space_b+(plot_box_height+space_i)  ;
Zplot_ypos     = space_b;
plot_box_width = total_plot_width-space_l-space_r ;  % points

figure (set_fig_nbr)
clf
set(gcf, 'Units', 'points', 'position', [5,55,total_plot_width,space_b+(plot_box_height+space_i)*3+ 12*0]);


% X-Gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TF_     = TF_X;
freq_   = F_f;

gdir  ='X';

ha_X = axes() ;
set(ha_X, 'Units', 'points', 'Position', [space_l, Xplot_ypos, plot_box_width, plot_box_height])

for n=1:1:length(ffbc)
  patch([ffbc(n)-ffbr(n)/2, ffbc(n)-ffbr(n)/2, ffbc(n)+ffbr(n)/2, ffbc(n)+ffbr(n)/2],[ymax+2,ymin-2,ymin-2,ymax+2],[ri,gi,bi], ...
         'linestyle','none','facealpha',alpha,'HandleVisibility',fhv)
  hold on
end
%hold on
set(gca,'ColorOrderIndex',1)
if plot_lin
  plot(freq_,            (abs(TF_)),'LineWidth',lws);
else
  semilogx(freq_,20*log10(abs(TF_)) ,'LineWidth',lws);
  set(gca, 'XScale', 'log')
end
xlim( [fmin   fmax])
ylim( [ymin ymax])
grid on
box on

ht_X = text(text_loc(1),text_loc(2),[gdir '-Gradient']);
xticklabels({})

% Y-Gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TF_     = TF_Y;
freq_   = F_f;

gdir  ='Y';

ha_Y = axes() ;
set(ha_Y, 'Units', 'points', 'Position', [space_l, Yplot_ypos, plot_box_width, plot_box_height])

for n=1:1:length(ffbc)
  patch([ffbc(n)-ffbr(n)/2, ffbc(n)-ffbr(n)/2, ffbc(n)+ffbr(n)/2, ffbc(n)+ffbr(n)/2],[ymax+2,ymin-2,ymin-2,ymax+2],[ri,gi,bi], ...
         'linestyle','none','facealpha',alpha,'HandleVisibility',fhv)
  hold on
end
%hold on
set(gca,'ColorOrderIndex',1)
if plot_lin
  plot(freq_,            (abs(TF_)),'LineWidth',lws);
else
  semilogx(freq_,20*log10(abs(TF_)) ,'LineWidth',lws);
  set(gca, 'XScale', 'log')
end
xlim( [fmin   fmax])
ylim( [ymin ymax])
grid on
box on

ht_Y = text(text_loc(1),text_loc(2),[gdir '-Gradient']);
xticklabels({})
ylabel(ylabel_str,'Interpreter', 'tex' )



% Z-Gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TF_     = TF_Z;
freq_   = F_f;

gdir  = 'Z';

ha_Z = axes() ;
set(ha_Z, 'Units', 'points', 'Position', [space_l, Zplot_ypos, plot_box_width, plot_box_height])

for n=1:1:length(ffbc)
  patch([ffbc(n)-ffbr(n)/2, ffbc(n)-ffbr(n)/2, ffbc(n)+ffbr(n)/2, ffbc(n)+ffbr(n)/2],[ymax+2,ymin-2,ymin-2,ymax+2],[ri,gi,bi], ...
         'linestyle','none','facealpha',alpha,'HandleVisibility',fhv)
  hold on
end
%hold on
set(gca,'ColorOrderIndex',1)
if plot_lin
  plot(freq_,            (abs(TF_)),'LineWidth',lws);
else
  semilogx(freq_,20*log10(abs(TF_)) ,'LineWidth',lws);
  set(gca, 'XScale', 'log')
end
xlim( [fmin   fmax])
ylim( [ymin ymax])
grid on
box on

ht_Z = text(text_loc(1),text_loc(2),[gdir '-Gradient']);
xlabel(xlabel_str,'Interpreter', 'tex' )


set(ha_X, 'FontName',FN_A,'FontSize',FS_A);
set(ht_X, 'FontName',FN_A,'FontSize',FS_A);
set(ha_Y, 'FontName',FN_A,'FontSize',FS_A);
set(ht_Y, 'FontName',FN_A,'FontSize',FS_A);
set(ha_Z, 'FontName',FN_A,'FontSize',FS_A);
set(ht_Z, 'FontName',FN_A,'FontSize',FS_A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



