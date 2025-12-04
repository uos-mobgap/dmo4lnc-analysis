function [timestamp_UNIX,StartDateTime] = timestamp_startdatetime_conversion(timestamp_mcr,TimeZone)
% This function converts McRoberts timestamps from "McR local UNIX" to UTC UNIX
% And creates the StartDateTime from the first timestamp (UTC UNIX)

% INPUT: 
% timestamp_mcr: timestamps from dynaport, that are in "McR local UNIX"
% TimeZone: local timezone retrived from escience or subject ID
% OUTPUT:
% timestamp_UNIX:  standard UTC UNIX timestamps
% StartDateTime: char vector with the start date-time of the acquisition in local time (format ISO 8601)

% First, we convert the first timestamp_mcr using UTC time zone(although we know it is not an actual UTC), and then we convert it in datenum
firstDayUTC = datenum(datetime(timestamp_mcr(1),'convertfrom','posixtime','format','yyyy-MM-dd','Timezone','utc'));

%Then we reconvert datenum in datetime assigning to it the local timezone 
%(because we know that timestamp_mcr(1) is a McR local UNIX and not a standard UNIX).
%In this way we are sure to obtain the correct start day of the acquisition
firstDayLocal = datetime(firstDayUTC,'convertfrom','datenum','format','yyyy-MM-dd','Timezone',TimeZone);

% from the start day we extract the offset from UTC with tzoffset() which also takes into account the daylight-saving time
offsetUTC = tzoffset(firstDayLocal);

% we convert the offset into seconds (UNIX format)
offsetUNIX = seconds(offsetUTC);

% we subtract the UNIX offset from timestamp_mcr (the full vector) and obtain a timestamp in the UTC UNIX standard
timestamp_UNIX = timestamp_mcr-offsetUNIX;

StartDateTime = char(datetime(timestamp_UNIX(1),'ConvertFrom','posixtime','Format','yyyy-MM-dd''T''HH:mm:ss.SSSXXX','TimeZone',TimeZone));
end

