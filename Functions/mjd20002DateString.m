function DateString = mjd20002DateString(mjd)

date = round(mjd20002date(mjd));

DateString = [num2str(date(3)), '\', num2str(date(2)), '\', num2str(date(1))]

end