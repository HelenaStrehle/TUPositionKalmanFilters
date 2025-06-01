function T = GeotrackerRef(filename)
    % Læs KML-filen
    xmlDoc = xmlread(filename);

    % Find gx:Track uden namespace
    trackElements = xmlDoc.getElementsByTagName('gx:Track');
    if trackElements.getLength() < 1
        error('Ingen gx:Track element fundet');
    end
    trackElem = trackElements.item(0);

    % Find <when> og <gx:coord>
    whenElements  = trackElem.getElementsByTagName('when');
    coordElements = trackElem.getElementsByTagName('gx:coord');

    N = whenElements.getLength();
    if coordElements.getLength() ~= N
        error('Mismatch i antal tider og koordinater');
    end

    % Init arrays
    times = cell(N,1);
    latitudes = zeros(N,1);
    longitudes = zeros(N,1);

    % Tidsstempler (ensret til .SSS)
    for k = 0:(N-1)
        t = strtrim(char(whenElements.item(k).getFirstChild.getData));
        if endsWith(t, 'Z') && ~contains(t, '.')
            t = strrep(t, 'Z', '.000Z');
        end
        times{k+1} = t;
    end

    % Konverter til datetime med millisekunder og tidszoner
    timeUTC_dt = datetime(times, ...
        'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''', ...
        'TimeZone', 'UTC');

    timeDK_dt = datetime(timeUTC_dt, 'TimeZone', 'Europe/Copenhagen');

    % Konverter til tekstformat UDEN 'T', men MED millisekunder
    timeUTC = cellstr(datestr(timeUTC_dt, 'yyyy-mm-dd HH:MM:SS.FFF'));
    timeDK  = cellstr(datestr(timeDK_dt,  'yyyy-mm-dd HH:MM:SS.FFF'));

    % Koordinater
    for k = 0:(N-1)
        coordStr = char(coordElements.item(k).getFirstChild.getData);
        parts = sscanf(coordStr, '%f');
        longitudes(k+1) = parts(1);
        latitudes(k+1)  = parts(2);
    end

    % Returner tabel
    T = table(timeUTC, timeDK, latitudes, longitudes, ...
        'VariableNames', {'TimeUTC','TimeDK','Latitude','Longitude'});
end
