
; aztrack
;
; Generate a track of exposures at constant elevation. The track moves
; in increasing azimuth, except for a random frequency of negative
; azimuth moves.
;
; INPUT
;  year,month,day: date of observation
;  hours: array of UT times in hours for each exposure
;     e.g., hours=[1.0,1.5,2.0] means the exposures start at 
;              01:00UT,01:30UT,02:00UT
;  elevation: the constant elevation for the track
; OUTPUT
;  ra,dec,az,utout: vectors containing output tracks
;                   utout is a (possibly) truncated version of 'hours',
;                   in that the sequence may reach a limit in ra/dec 
; PARAMETERS
;  jumpFrac: fraction of exposures which get a negative AZ offset
;            default is 0.1 (10%). The negative offsets are at 3*offsetScale
;  offsetScale: maximum AZ offset between each exposure, in degrees
;            default is 3
;  decMin: The minimum declination allowed (default is 30), in degrees
;  decMax:     maximum   "             "   (default is 75)
;  raMin,raMax: same for RA, defaults are 30,300
;  gbMin: The minimum galactic latitude allowed (default is 17)
;  seed: fix the random seed for the track
;
pro aztrack,year,month,day,hours,elevation,ra,dec,az,utout, $
      jumpFrac=jumpFrac,offsetScale=offsetScale, $
      decMin=decMin,decMax=decMax,raMin=raMin,raMax=raMax, $
      gbMin=gbMin,seed=seed
	if not keyword_set(jumpFrac) then jumpFrac = 0.1
	if not keyword_set(offsetScale) then offsetScale = 3.0
	if not keyword_set(decMin) then decMin = 30.
	if not keyword_set(decMax) then decMax = 75.
	if not keyword_set(raMin) then raMin = 30.
	if not keyword_set(raMax) then raMax = 300.
	if not keyword_set(gbMin) then gbMin = 17.
	;if not keyword_set(seed) then seed = !NULL
	alt = elevation
	; convert the exposure start times to julian dates
	jdcnv,year,month,day,hours,julian
	; select a random declination from within the allowed range, this will
	; be the starting dec
	x1 = sin(decMin*!DTOR)
	x2 = sin(decMax*!DTOR)
	y = randomu(seed)
	decStart = asin(x1+(x2-x1)*y) * !RADEG
	; find the allowed RA range at this declination, by looking for the
	; boundaries in galactic latitude
	dRA = 0.1
	npts = long((raMax-raMin)/dRA)+1
	rab = raMin + dRA*findgen(npts)
	decb = replicate(decStart,npts)
	glactc,rab,decb,2000.,gl,gb,1,/degree
	ii = where(abs(gb) gt gbMin)
	i1 = ii[0]
	ii = where(abs(gb[i1:npts-1]) lt gbMin)
	if ii[0] lt 0 then begin
		i2 = npts-1
	endif else begin
		i2 = ii[0]-1
	endelse
	ra1 = rab[i1]
	ra2 = rab[i2]
	; convert the ra range at the fixed declination to an alt/az range,
	; then find the azimuth where the altitude matches the desired one
	x = rab[i1:i2]
	y = decb[i1:i2]
	eq2hor,x,y,julian[0],tmpalt,tmpaz,obsname='kpno'
	delalt = min(abs(tmpalt-alt),i)
	foo = where(abs(tmpalt-alt) lt 0.1)
	if delalt gt 0.1 then begin
		; there is no match, this track failed
		ra = [-9999.99]
		return
	endif
	azStart = tmpaz[i]
	; starting at the first azimuth, form a track by making random offsets
	; in azimuth, each one up to offsetScale. A specified fraction of the
	; moves (jumpFrac) will be in the negative direction at 3*offsetScale
	az = 0*hours
	az[0] = azStart
	for i=1,n_elements(az)-1 do begin
		if (randomu(seed) lt jumpFrac) then begin
			az[i] = az[i-1] - 3*offsetScale*randomu(seed)
		endif else begin
			az[i] = az[i-1] + offsetScale*randomu(seed)
		endelse
	endfor
	; convert the azimuth track to ra/dec
	hor2eq,alt,az,julian,ra,dec,obsname='kpno'
	; sometimes it can come out just below the allowed cut
	if dec[0] lt decMin then begin
		ra = [-9999.99]
		return
	endif
	; cut off the track when it extends out of the dec range
	ii = where((dec lt decMin) or (dec gt decMax))
	if ii[0] lt 0 then begin
		last = n_elements(ra)-1
	endif else begin
		last = ii[0]-1
	endelse
	ra = ra[0:last]
	dec = dec[0:last]
	utout = hours[0:last]
	; or extends beyond the specified RA range
	ii = where(ra gt raMax)
	if ii[0] lt 0 then begin
		last = n_elements(ra)-1
	endif else begin
		last = ii[0]-1
	endelse
	ra = ra[0:last]
	dec = dec[0:last]
	utout = hours[0:last]
end

