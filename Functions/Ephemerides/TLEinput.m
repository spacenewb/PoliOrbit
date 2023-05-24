function [ a, e, i, OM, om, theta, otherData ] = TLEinput(fname, type)
% TLE2Kp Inputs TLE data from file into MatLab
% The various TLE's can either be separated by the s/c names, OR they can
% be all attached to eachother
% mu = 398600 km^3/s^2
%
% otherData is a struct:
%   otherData(i).field
%
% Example:
% ---------------------------------------------------------------------
% SOYUZ 19
% 1 49269U 21089A   21287.53081300  .00005594  00000-0  11084-3 0  9990
% 2 49269  51.6413 112.6192 0004180 111.0667 297.8505 15.48678783307091
% ISS
% 1 49269U 21089A   21287.53081300  .00005594  00000-0  11084-3 0  9990
% 2 49269  51.6413 112.6192 0004180 111.0667 297.8505 15.48678783307091
% CREW-DRAGON
% 1 49269U 21089A   21287.53081300  .00005594  00000-0  11084-3 0  9990
% 2 49269  51.6413 112.6192 0004180 111.0667 297.8505 15.48678783307091
%
% OR
%
% 1 49269U 21089A   21287.53081300  .00005594  00000-0  11084-3 0  9990
% 2 49269  51.6413 112.6192 0004180 111.0667 297.8505 15.48678783307091
% 1 49269U 21089A   21287.53081300  .00005594  00000-0  11084-3 0  9990
% 2 49269  51.6413 112.6192 0004180 111.0667 297.8505 15.48678783307091
% 1 49269U 21089A   21287.53081300  .00005594  00000-0  11084-3 0  9990
% 2 49269  51.6413 112.6192 0004180 111.0667 297.8505 15.48678783307091
% ---------------------------------------------------------------------
%
% PROTOTYPE:
%       [ OrbitalElements, otherData ] = TLEinput(fname, type)
%
% INPUTS:
%       fileName            String containing the name of the .txt file
%       type                String 'nameseparated' or 'unseparated'
%                               If nameseparated is specified, function
%                               assumes file is like the first example with
%                               TLE's separated only by one line with s/c
%                               names. If unseparated is specified,
%                               function assumes there is no separation
%                               between TLE's
%
% OUTPUT:
%       a           [km]    Semi-major axis
%       e           [-]     Eccentricity
%       i           [deg]   Inclination
%       OM          [deg]   Right ascension of the ascending node
%       om          [deg]   Argument of pericenter
%       theta       [deg]   True anomaly -> (NOTE: THIS IF ONLY FOUND FOR THE CASE OF ELLIPTIC ORBITS)
%       otherData           Complete data stored in TLE (NOTE: function is still missing some parameters)
%
%
%
% CREATOR:
%       Matteo D'Ambrosio

%% File

fid = fopen(fname, 'rb');

if fid == -1
    error('Failed to open file --> Check if name is correct and includes .txt')
end

% Get length of file
len = 0;
while ~feof(fid)
    fgetl(fid);
    len = len + 1;
end
fclose(fid);

%% Get elements

fid = fopen(fname, 'rb');

switch type

    case 'unseparated'

        % Check length
        TLE_num = len/2;
        if mod(TLE_num,1)
            error('Error with file type --> Check ''type'' input or TLE file')
        end

        % Create vectors

        a = zeros(length(TLE_num),1);
        e = zeros(length(TLE_num),1);
        i = zeros(length(TLE_num),1);
        OM = zeros(length(TLE_num),1);
        om = zeros(length(TLE_num),1);
        theta = zeros(length(TLE_num),1);

        % Get strings
        for k = 1:TLE_num

            l1 = fgetl(fid);
            sc(k).sat_number = str2double(l1(3:7));
            sc(k).class = l1(8);
            sc(k).launch_year = str2double(l1(10:11));
            sc(k).launch_number = str2double(l1(12:14));
            sc(k).launch_piece = l1(15:17);
            sc(k).epoch_year = str2double(l1(19:20));
            sc(k).epoch = str2double(l1(21:32)) * 24 * 3600;
            % Missing first derivative of mean motion
            % Missing second derivative of mean motion
            % Missing BSTAR drag term
            sc(k).ephemeris_type = str2double(l1(63));
            sc(k).element_number = str2double(l1(65:68));

            l2 = fgetl(fid);
            sc(k).i = str2double(l2(9:16));
            sc(k).OM = str2double(l2(18:25));
            sc(k).e = str2double(l2(27:33))/1e7;
            sc(k).om = str2double(l2(35:42));
            sc(k).M = str2double(l2(44:51));
            sc(k).n = str2double(l2(53:63));
            sc(k).rev_number_atEpoch = str2double(l2(64:68));

            % Find theta from mean anomaly
            M = deg2rad(sc(k).M);
            FUN = @(E) E - sc(k).e .* sin(E) - M;
            E_guess = M + sc(k).e .* sin( M ) ./ ( 1 - sin( M + sc(k).e ) + sin( M ) );
            E = fzero(FUN,E_guess);
            arg = sqrt(( 1 + sc(k).e )./( 1 - sc(k).e ) ) .* tan(E/2);


            mu = astroConstants(13);
            % Create orbital element vectors
            a(k) = (mu/(sc(k).n*2*pi/(24*3600)).^2).^(1/3); sc(k).a = a(k);
            e(k) = sc(k).e;
            i(k) = sc(k).i;
            OM(k) = sc(k).OM;
            om(k) = sc(k).om;
            theta(k) = rad2deg(2 * atan2( sin(arg), cos(arg) )); sc(k).theta = theta(k);

        end

        a=a';
        e=e';
        i=i';
        OM=OM';
        om=om';
        theta=theta';

    case 'nameseparated'

        % Check length
        TLE_num = (len - len/3)/2;
        if mod(TLE_num,1)
            error('Error with file type --> Check ''type'' input or TLE file')
        end

        % Create vectors

        a = zeros(length(TLE_num),1);
        e = zeros(length(TLE_num),1);
        i = zeros(length(TLE_num),1);
        OM = zeros(length(TLE_num),1);
        om = zeros(length(TLE_num),1);
        theta = zeros(length(TLE_num),1);

        % Get strings
        for k = 1:TLE_num

            scname = fgetl(fid);
            sc(k).name = replace(scname,' ','-');

            l1 = fgetl(fid);
            sc(k).sat_number = str2double(l1(3:7));
            sc(k).class = l1(8);
            sc(k).launch_year = str2double(l1(10:11));
            sc(k).launch_number = str2double(l1(12:14));
            sc(k).launch_piece = l1(15:17);
            sc(k).epoch_year = str2double(l1(19:20));
            sc(k).epoch = str2double(l1(21:32)) * 24 * 3600;
            % Missing first derivative of mean motion
            % Missing second derivative of mean motion
            % Missing BSTAR drag term
            sc(k).ephemeris_type = str2double(l1(63));
            sc(k).element_number = str2double(l1(65:68));

            l2 = fgetl(fid);
            sc(k).i = str2double(l2(9:16));
            sc(k).OM = str2double(l2(18:25));
            sc(k).e = str2double(l2(27:33))/1e7;
            sc(k).om = str2double(l2(35:42));
            sc(k).M = str2double(l2(44:51));
            sc(k).n = str2double(l2(53:63));
            sc(k).rev_number_atEpoch = str2double(l2(64:68));

            % Find theta from mean anomaly
            M = deg2rad(sc(k).M);
            FUN = @(E) E - sc(k).e .* sin(E) - M;
            E_guess = M + sc(k).e .* sin( M ) ./ ( 1 - sin( M + sc(k).e ) + sin( M ) );
            E = fzero(FUN,E_guess);

            arg = sqrt(( 1 + sc(k).e )./( 1 - sc(k).e )  ) .* tan(E/2);

            mu = astroConstants(13);
            % Create orbital element vectors

            a(k) = (mu/(sc(k).n*2*pi/(24*3600)).^2).^(1/3); sc(k).a = a(k);
            e(k) = sc(k).e';
            i(k) = sc(k).i';
            OM(k) = sc(k).OM';
            om(k) = sc(k).om';
            theta(k) = rad2deg(2 * atan2( sin(arg), cos(arg) ))'; sc(k).theta = theta(k);

        end

        a=a';
        e=e';
        i=i';
        OM=OM';
        om=om';
        theta=theta';

    otherwise

        error('Type must be a string containing ''nameseparated'' or ''unseparated''')

end

% Close file
fclose(fid);
otherData = sc;

end