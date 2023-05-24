function [ planet_object ] = planet3D_ADD(planetNum, coords, figureNum, multiplier)
% planet3D Creates a 3D model of the planet in a pre-existing figure
%
% PROTOTYPE: 
%   [] = planet3D_ADD(planetNum, coords, figureNum, multiplier)
%
% INPUT:
%   planetNum[1x1]  Number from 1 to 11, to choose which body to display:
%             1  Mercury      
%             2  Venus       
%             3  Earth       
%             4  Mars      
%             5  Jupiter       
%             6  Saturn       
%             7  Uranus       
%             8  Neptune       
%             9  Pluto       
%             10 Moon  
%             11 Sun
%   coords[1x3]     If specified, centers the planet at coordinates [ x, y, z ] [km]
%   figureNum[1x1]  Figure handle of the figure to add the planet to
%
% OUTPUT: 
%   planet_object[1]        Object used to rotate the planet with rotate(planet_object, ...)
% 
% CALLED FUNCTIONS:
%   astroConstants();
%
% REFERENCES:
%   Mercury: https://stevealbers.net/albers/sos/saturn/mimas/mimas_rgb_cyl_www.jpg
%   Venus: https://spacescenting.com/wp-content/uploads/2019/05/planet-journey.jpg
%   Earth: https://www.solarsystemscope.com/textures/download/8k_earth_daymap.jpg
%   Mars: https://www.solarsystemscope.com/textures/download/2k_mars.jpg
%   Jupiter: https://d2pn8kiwq2w21t.cloudfront.net/original_images/jpegPIA19643.jpg
%   Saturn: https://www.solarsystemscope.com/textures/download/2k_saturn.jpg
%   Uranus: https://static.wikia.nocookie.net/planet-texture-maps/images/b/b3/Uranus-2.jpg
%   Neptune: https://static.wikia.nocookie.net/planet-texture-maps/images/c/c1/Planetarium_neptune.jpg
%   Pluto: https://scx2.b-cdn.net/gfx/news/hires/2018/5b44d20643ef6.jpg
%   Moon: https://p0.pikist.com/photos/36/520/map-moon-world-satellite.jpg
%   Sun: https://www.solarsystemscope.com/textures/download/2k_sun.jpg
%--------------------------------------------------------------------------

switch planetNum % Set radius of planet, and select picture
    
    case 1, % Mercury
    R = astroConstants(21); 
    image = 'https://www.solarsystemscope.com/textures/download/2k_mercury.jpg';
    
    case 2, % Venus
    R = astroConstants(22); 
    image = 'https://www.solarsystemscope.com/textures/download/2k_venus_atmosphere.jpg';
    
    case 3, % Earth
    R = astroConstants(23);
    w = 7.2921159e-5;
    image = 'https://www.solarsystemscope.com/textures/download/2k_earth_daymap.jpg';
    
    case 4, % Mars
    R = astroConstants(24); 
    image = 'https://www.solarsystemscope.com/textures/download/2k_mars.jpg';
    
    case 5, % Jupiter
    R = astroConstants(25); 
    image = 'https://www.solarsystemscope.com/textures/download/2k_jupiter.jpg';
    
    case 6, % Saturn
    R = astroConstants(26); 
    image = 'https://www.solarsystemscope.com/textures/download/2k_saturn.jpg';
    
    case 7,% Uranus
    R = astroConstants(27); 
    image = 'https://www.solarsystemscope.com/textures/download/2k_uranus.jpg';
    
    case 8, % Neptune
    R = astroConstants(28); 
    image = 'https://www.solarsystemscope.com/textures/download/2k_neptune.jpg';
    
    case 9, % Pluto
    R = astroConstants(29); 
    image = 'https://scx2.b-cdn.net/gfx/news/hires/2018/5b44d20643ef6.jpg';
    
    case 10, % Sun
    R = 696340;
    image = 'https://www.solarsystemscope.com/textures/download/2k_sun.jpg';
    
    case 11, % Moon
    R = astroConstants(30); 
    image = 'https://p0.pikist.com/photos/36/520/map-moon-world-satellite.jpg';
    
        
    otherwise
    error('Input ''planetNum'' does not correspond to an object')
    
end

if nargin == 3
    multiplier = 1;
end

R = R * multiplier;

%% Check inputs and create figure

if ishandle(figureNum)
    figure(figureNum);
else
    error('Specified figure with handle figureNum does not exist');
end

if ~isequal(size(coords), [1,3]) 
    close
    error('Coordinates vector must be [1x3]')
end

%% Figure

hold on
grid on
axis equal;

%Set axis labels
xlabel('x [km]'); 
ylabel('y [km]'); 
zlabel('z [km]'); 

view(120, 30); %Set initial view

%% Create 3D surface for mapping

[x, y, z] = ellipsoid( coords(1), coords(2), -coords(3), R, R, R); % Create meshgrid for sphere

planet_object = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none', 'HandleVisibility', 'off'); % Create sphere as a surface

%% Map texture on the surface

image = imread(image); % Load image from website

alpha = 1; % Transparency of globe: 1 = opaque, 0 = invisible

% 'Facecolor' set to 'texturemap' lets image be mapped on the sphere, with
% the image coming from the 'CData' property 

set(planet_object, 'FaceColor', 'texturemap', 'CData', image, 'FaceAlpha', alpha, 'EdgeColor', 'none')

drawnow;

end