function [longdeg,latdeg] = nmea0183_to_degrees(longitude, latitude)
%nmea_to_degrees 
%   Converts given longitude and latitude nmea values to
%   degrees
longdeg = floor(longitude/100) + (longitude - floor(longitude/100)*100)/60;
latdeg = floor(latitude/100) + (latitude - floor(latitude/100)*100)/60;
end

