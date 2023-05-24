function TLE = TLE_read(fileName)

% FieldNames = {'SatName',...
%            'Ln-1', 'CatNum', 'Class', 'Intl. Desig. LnchYr', 'Intl. Desig. LnchNum', 'Intl. Desig. LnchPiece', 'EpochYr', 'EpochDay', 'MeanMot1D', 'MeanMot2D', 'B*', 'EphType', 'ElmntSetNum', 'Chksm-1',...
%            'Ln-2', 'CatNum', 'InclDeg', 'RAANDeg', 'Ecc', 'ArgPeriDeg', 'MeanAnomDeg', 'MeanMotRevDay', 'RevNumEpoch', 'Chksm-2'};

% Example:
% ISS (ZARYA)
% 1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
% 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537

% Field   Columns     Content                                                                             Example
% 1	        01–24	    Satellite name	                                                                    ISS (ZARYA)

% Field   Columns     Content                                                                             Example
% 1	        01	        Line number	                                                                        1
% 2	        03–07	    Satellite catalog number	                                                        25544
% 3	        08	        Classification (U: unclassified, C: classified, S: secret) 	                        U
% 4	        10–11	    International Designator (last two digits of launch year)	                        98
% 5	        12–14	    International Designator (launch number of the year)	                            067
% 6	        15–17	    International Designator (piece of the launch)	                                    A
% 7	        19–20	    Epoch year (last two digits of year)	                                            08
% 8	        21–32	    Epoch (day of the year and fractional portion of the day)	                        264.51782528
% 9	        34–43	    First derivative of mean motion; the ballistic coefficient	                        -.00002182
% 10	    45–52	    Second derivative of mean motion (decimal point assumed)	                        00000-0
% 11	    54–61	    B*, the drag term, or radiation pressure coefficient (decimal point assumed)	    -11606-4
% 12	    63–63	    Ephemeris type (always zero; only used in undistributed TLE data)	                0
% 13	    65–68	    Element set number. Incremented when a new TLE is generated for this object.	    292
% 14	    69–69	    Checksum (modulo 10)	                                                            7

% Field   Columns     Content                                                                             Example
% 1	        01	        Line number	                                                                        2
% 2	03–07	Satellite Catalog number	                                                                    25544
% 3	09–16	Inclination (degrees)	                                                                        51.6416
% 4	18–25	Right ascension of the ascending node (degrees)	                                                247.4627
% 5	27–33	Eccentricity (decimal point assumed)	                                                        0006703
% 6	35–42	Argument of perigee (degrees)	                                                                130.5360
% 7	44–51	Mean anomaly (degrees)	                                                                        325.0288
% 8	53–63	Mean motion (revolutions per day)	                                                            15.72125391
% 9	64–68	Revolution number at epoch (revolutions)	                                                    56353
% 10	69	Checksum (modulo 10)	                                                                        7

% Where decimal points are assumed, they are leading decimal points.
% (-11606-4) translates to −0.11606E−4 (−0.11606×10−4).

% TLE data uses a very simple checksum algorithm - basically the sum of all 
% the digits in the line mod 10, with minus signs counting as 1, letters 
% and whitespace as 0

% Two-digit Epoch Years from 57-99 correspond to 1957-1999 and those 
% from 00-56 correspond to 2000-2056.

%% Read File
fileID = fopen(fileName);
Lines = repmat('',3,69);
for i = 1:3
    Currline = fgetl(fileID);
    Lines(i,1:length(Currline)) = Currline;
end
fclose('all');

%% Checksum
letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
Checksum = zeros([1,3]);

for k = 2:3
    line = Lines(k,:);
    cksum = 0;
    for i = 1:68
        c = line(i);
        if c == ' ' || c == '.' || c == '+' || ismember(c, letters)
            continue
        elseif c == '-'
            cksum = cksum + 1;
        else
            cksum = cksum + str2num(c);
        end
    end
    Checksum(k) = mod(cksum,10);
end

if Checksum(2) == str2num(Lines(2,end)) && Checksum(3) == str2num(Lines(3,end))
    disp('Checksum Verified');
else
    disp('Checksum Verification Error');
end

%% TLE as Cell
TLE = repmat({''}, 3, 14);
% Line 0
TLE(1,1) = {deblank(Lines(1,:))}; % Non-Numeric
% Line 1
TLE(2,1) = {Lines(2,1)};        TLE(2,2) = {Lines(2,3:7)};      TLE(2,3) = {Lines(2,8)}; % Non-Numeric
TLE(2,4) = {Lines(2,10:11)};    TLE(2,5) = {Lines(2,12:14)};    TLE(2,6) = {deblank(Lines(2,15:17))}; % Non-Numeric
TLE(2,7) = {Lines(2,19:20)};    TLE(2,8) = {Lines(2,21:32)};    TLE(2,9) = {Lines(2,34:43)};
TLE(2,10) = {Lines(2,45:52)}; % Assumed Decimal
TLE(2,11) = {Lines(2,54:61)}; % Assumed Decimal
TLE(2,12) = {Lines(2,63)};      TLE(2,13) = {Lines(2,65:68)};   TLE(2,14) = {Lines(2,69)};
% Line 2
TLE(3,1) = {Lines(3,1)};        TLE(3,2) = {Lines(3,3:7)};      TLE(3,3) = {Lines(3,9:16)};
TLE(3,4) = {Lines(3,18:25)};    TLE(3,5) = {Lines(3,27:33)};    TLE(3,6) = {Lines(3,35:42)};
TLE(3,7) = {Lines(3,44:51)};    TLE(3,8) = {Lines(3,53:63)};    TLE(3,9) = {Lines(3,64:68)};
TLE(3,10) = {Lines(3,69)};

%% Assumed Decimal
for i = 10:11
    Temp_char = char(TLE(2,i));
    idxNeg = strfind(Temp_char,'-');
    if length(idxNeg) == 1 && ~ismember(1,idxNeg)
        Temp_char = strrep(Temp_char,'-','E-');
    elseif length(idxNeg) > 1
        Temp_char(idxNeg(2)) = 'k';
        Temp_char = strrep(Temp_char,'k','E-');
    end
    TLE(2,i) = {strrep(Temp_char, '.', '0.')};
end

TLE(2,9) = {strrep(TLE(2,9), '.', '0.')};
TLE(3,5) = {strcat('0.',string(TLE(3,5)))};

%% Non-Numeric
TLE(2,3) = {cellfun(@double ,TLE(2,3))};
TLE(2,6) = {cellfun(@double ,TLE(2,6))};

TLE(2:3,:) = cellfun(@(x) str2num(string(x)), TLE(2:3,:),'UniformOutput', false);
TLE(1,2:end) = {''};

TLE(2,3) = {cellfun(@char ,TLE(2,3))};
TLE(2,6) = {cellfun(@char ,TLE(2,6))};

end