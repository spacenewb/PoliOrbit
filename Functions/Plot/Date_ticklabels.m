function [] = Date_ticklabels( figureNum, Date_xmin, Date_xmax, Date_ymin, Date_ymax )
% Adds ticklabels with dates to plot
%
%   PROTOTYPE:
%   date_ticklabels( figureNum, Date_xmin, Date_xmax, Date_ymin, Date_ymax )

figure(figureNum);
hold on;

mjd_xmin = date2mjd2000(Date_xmin);
mjd_xmax = date2mjd2000(Date_xmax);
mjd_ymin = date2mjd2000(Date_ymin);
mjd_ymax = date2mjd2000(Date_ymax);

N = 10; % Number of x-labels
x_ = linspace( mjd_xmin, mjd_xmax, N );   % Label position
y_ = linspace( mjd_ymin, mjd_ymax, N );   % Label position

tick_label_x = cell( 1, N );
for s = 1:N
    date_array_x = mjd20002date( x_(s) );
    tick_label_x{s} = datestr( date_array_x,'dd mmm yyyy' );
    date_array_x = datevec( tick_label_x{s},'dd mmm yyyy' );
    x_(s) = date2mjd2000( date_array_x );
end

tick_label_y = cell( 1, N );
for s = 1:N
    date_array_y = mjd20002date( y_(s) );
    tick_label_y{s} = datestr( date_array_y, 'dd mmm yyyy' );
    date_array_y = datevec( tick_label_y{s}, 'dd mmm yyyy' );
    y_(s) = date2mjd2000( date_array_y );
end

xticks( x_ );
xticklabels( tick_label_x )
xtickangle( 45 )
yticks( y_ );
yticklabels( tick_label_y )

end