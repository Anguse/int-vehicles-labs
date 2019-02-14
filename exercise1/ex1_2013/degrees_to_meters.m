function [flon,flat] = degrees_to_meters(longdeg,latdeg)
%DEGREES_TO_METERS Summary of this function goes here
%   Detailed explanation goes here
flon = longdeg * 62393;
flat = latdeg * 111342;
end

