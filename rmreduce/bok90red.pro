
function filelist, imprefix, startframe, endframe, $
                   exclude_frames=exclude_frames, gzip=gzip
	;
	if keyword_set(gzip) then sfx='.gz' else sfx=''
	return_list = strarr(endframe-startframe+1)
	for i = startframe,endframe,1 do begin
		if keyword_set(exclude_frames) then begin
			reject = 0
			for j = 0,n_elements(exclude_frames)-1,1 do begin
				if i eq exclude_frames[j] then begin
					reject = 1
					break
				endif
			endfor
			if reject then continue
		endif
		i1 = i - startframe
		return_list[i1] = imprefix + string(i,format='(i4.4)') + '.fits' + sfx
	endfor
	return_list = return_list[where(return_list ne '')]
	return, return_list
end

;
; return a list of all UT dates with Bok observations, based on scanning
; the log directory
;
pro get_all_utdates,utdates
	logdir = getenv('BOK90PRIMEDIR')+'/../../github/rapala/survey/logs/'
	logfiles = file_search(logdir+'bass.*.log')
	logfiles = file_basename(logfiles,'.log')
	utdates = strmid(logfiles,6)
end

;
; read a single Bok observing log, returns some meta-data about each image
;
function read_log, utdate, all=all
	logdir = getenv('BOK90PRIMEDIR')+'/../../github/rapala/survey/logs/'
	logfn = logdir+'bass.'+utdate+'.log'
	readcol,logfn, $
	        fn,ut,filter,exptime,imtype,objname,airmass, $
	        format='A,A,A,F,A,A,F'
	w = where(strmid(fn,0,1) ne '#',nents)
	log = replicate({logent, fn:'',ut:'',filter:'',exptime:0.0, $
	                         imtype:'',objname:'',airmass:0.0, $
	                         reject:0},nents)
	log.fn = fn[w]
	log.ut = ut[w]
	log.filter = filter[w]
	log.exptime = exptime[w]
	log.imtype = imtype[w]
	log.objname = objname[w]
	log.airmass = airmass[w]
	; XXX temporary - should be read from file
;	reject_list = replicate({rejent, utdate:'',fn:''},2)
;	reject_list[0].utdate = '20140114'
;	reject_list[0].fn = 'bokrm.20140114.0001'
;	reject_list[1].utdate = '20140319'
;	reject_list[1].fn = 'bokrm.20140319.0026'
;	; XXX
;	for i = 0,n_elements(reject_list)-1,1 do begin
;		j = where((utdate eq reject_list[i].utdate) and $
;		          (log.fn eq reject_list[i].fn))
;		if j[0] ne -1 then begin
;			log[j].reject = 1
;		endif
;	endfor
	return,log
end

;
; default aperture radii in pixels 
;
pro get_default_apers,aperrad,skyrad
	aperrad=findgen(6)*1.5+2. ; range is 0.9''..4.3'' in steps of 0.7''
	aperrad=[aperrad,15.] ; add 6.8'' aperture as largest
	skyrad=[22.,33.] ; 10'' -> 15''
end

pro reduce_night, utdate, filters=filters, $
                  name=name,filename=filename, $
                  skipifexists=skipifexists, gzip=gzip, cleanup=cleanup

	if not keyword_set(skipifexists) then skipifexists=1
	if not keyword_set(gzip) then gzip=1

	rawdir = getenv('BOK90PRIMERAWDIR')+'/'+utdate+'/'
	print, rawdir

	datadir = getenv('BOK90PRIMEDIR')+'/data/'

	reduxdir = getenv('BOK90PRIMEOUTDIR')
	outdir = reduxdir+'/'+utdate+'/'
	print, outdir
	if file_test(outdir,/dir) eq 0 then file_mkdir,outdir

	log = read_log(utdate)

	; default to all filters used
	if not keyword_set(filters) then begin
		filters = log[uniq(log.filter)].filter
		index=where(filters eq 'g' or filters eq 'bokr')
		filters=filters[index]
	endif

	bias_list = where((log.imtype eq 'zero') and $
	                  (log.reject eq 0))
	if bias_list[0] ne -1 then begin
		bias_list = log[bias_list].fn + '.fits'
		if keyword_set(gzip) then bias_list += '.gz'
		bok4_mkbias, bias_list, $
		             rawdir=rawdir, outdir=outdir, skipifexists=skipifexists
		biasfile = outdir+'bias.fits'
	endif else begin
		; no biases taken this night. go backward in UT until one is found
		date = long(utdate)
		while date gt 20150101 do begin
			utdate2 = string(date,format='(i0)')
			utdir2 = getenv('BOK90PRIMEOUTDIR')+'/ut'+utdate2+'/'
			biasfile = utdir2+'bias.fits'
			if file_test(biasfile,/read) then break
			date -= 1
		endwhile
	endelse

	domeflats = strarr(n_elements(filters))

	for i=0,n_elements(filters)-1,1 do begin
		flat_list = where((log.imtype eq 'flat') and $
		                  (log.filter eq filters[i]) and $
		                  (log.reject eq 0))
		if flat_list[0] eq -1 then begin
			; no dome flats in this filter taken this night. 
			; go backward in UT until one is found.
			date = long(utdate)
			print,'no dome flat for '+filters[i]+' for '+utdate,date
			while date gt 20150101L do begin
				utdate2 = string(date,format='(i0)')
				utdir2 = getenv('BOK90PRIMEOUTDIR')+'/'+utdate2+'/'
				flatfile = utdir2+'domeflat_'+filters[i]+'.fits'
				if file_test(flatfile,/read) then begin
					domeflats[i] = flatfile
					break
				endif
				date -= 1
			endwhile
			print,'set to dome flat '+utdate2
			spawn,'cp '+getenv('BOK90PRIMEOUTDIR')+'/'+utdate2+$
				'/domeflat_'+filters[i]+'.fits '+$
				getenv('BOK90PRIMEOUTDIR')+'/'+utdate+$
				'/domeflat_'+filters[i]+'.fits'
		endif else begin
			flat_list = log[flat_list].fn + '.fits'
			if keyword_set(gzip) then flat_list += '.gz'
			bok4_mkflat, flat_list, 'domeflat_'+filters[i]+'.fits', $
			             biasname=biasfile, $
			             rawdir=rawdir, outdir=outdir, $
			             skipifexists=skipifexists
			domeflats[i] = outdir+'domeflat_'+filters[i]+'.fits'
		endelse
	endfor

;	bok4_mkbpm, flatname=outdir+'domeflat_'+filters[i]+'.fits', $
;	            outname=outdir+'bpm.fits', $
;	            skipifexists=skipifexists
;	spawn,'cp '+outdir+'ut'+utdate+'/bpm.fits '+outdir+'bpm.fits'

	for i=0,n_elements(filters)-1,1 do begin
		bok4_mkbpm, flatname=outdir+'domeflat_'+filters[i]+'.fits', $
                    outname=outdir+'bpm.fits', $
                    skipifexists=skipifexists

		if keyword_set(name) then begin
			sci_ims = where((log.imtype eq 'object') and $
			                (log.filter eq filters[i]) and $
			                (log.objname eq name) and $
			                (log.reject eq 0))
		endif else if keyword_set(filename) then begin
			sci_ims = where((log.imtype eq 'object') and $
			                (log.filter eq filters[i]) and $
			                (log.fn eq filename) and $
			                (log.reject eq 0))
		endif else begin
			sci_ims = where((log.imtype eq 'object') and $
			                (log.filter eq filters[i]) and $
			                (log.reject eq 0))
		endelse
		if sci_ims[0] eq -1 then continue
		sci_ims = log[sci_ims].fn
		sci_ims_raw = sci_ims + '.fits'
		if keyword_set(gzip) then sci_ims_raw += '.gz'
		print,'objects: ',filters[i],'  ',sci_ims

		;Reduce bias, flat, sky & bleeding points
		bok4_ccdproc1, sci_ims_raw, $
		               biasname=biasfile, flatname=domeflats[i], $
		               tnorm=50., $
		               rawdir=rawdir, outdir=outdir, skipifexists=skipifexists

		flat_list = sci_ims[0:n_elements(sci_ims)-1:2]
		print,flat_list
;		bok4_mkflat, flat_list, 'skyflat_g.fits', $
;		             biasname='bias.fits', $
;		             rawdir=rawdir, outdir=outdir, skipifexists=skipifexists

		if filters[i] eq 'i' or filters[i] eq 'z' then begin
			frgname = outdir+'fringe_'+filters[i]+'.fits'
			nsci = n_elements(sci_ims)
			bok4_mkfrg, sci_ims[0:nsci-1:4], frgname, $
			            bpmname=outdir+'/bpm.fits', $
			            outdir=outdir+'ccdproc1', skipifexists=skipifexists

                	;fix fringe, mask, and sky
                	bok4_ccdproc2, sci_ims, $
                               bpmname=outdir+'/bpm.fits', frgname=frgname, $
                               outdir=outdir, skipifexists=skipifexists
 
		endif else begin
			frgname = 0
			;Fringe is superflat
			flatname=outdir+'superflt_'+filters[i]+'.fits'
			nsci=min([n_elements(sci_ims),30])
			bok4_mkfrg,sci_ims[0:nsci-1],flatname, $
				bpmname=outdir+'/bpm.fits',$
				outdir=outdir+'ccdproc1',skipifexists=skipifexists
			
			; fix superflat, mask and sky			
			bok4_ccdproc2_superflt,sci_ims,$
				bpmname=outdir+'/bpm.fits',flatname=flatname, $
				outdir=outdir, skipifexists=skipifexists	
		endelse
		

		bok4_combineccds, sci_ims, $
		               outdir=outdir, skipifexists=skipifexists

;		residual=datadir+['xmap.fits','ymap.fits']

;		bok_astrom,sci_ims,sdssfn=datadir+'sdss.fits',$
;			ucac4fn=datadir+'ucac4.fits', $
;		        scale=8000.,outdir=outdir,/sex,residual=residual,$
;			skipifexists=skipifexists

;		get_default_apers,aperrad,skyrad

;		bok_photom,sci_ims,sdssfn=datadir+'sdss.fits',aperrad=aperrad,$
;			skyrad=skyrad,outdir=outdir,scale=8000.,/ps,$
;			skipifexists=skipifexists,/setslope
		
;		targetfn=datadir+'target_fibermap.fits'
;		rmaper,sci_ims,targetfn=targetfn,outdir=outdir,$
;		    scale=8000.,/rm,/all,residual=residual,skipifexists=skipifexists

	endfor

	if keyword_set(cleanup) then begin
		spawn, 'rm -rf '+outdir+'ccdproc[12]'
	endif

end


pro all_pho,astrom=astrom,photom=photom,aperphot=aperphot, $
            utrange=utrange,redo=redo
;  utrange can be a single utdate string (e.g., '20140414')
;          or a month (e.g., '201404')
;          or any such substring to restrict which nights are reduced
if keyword_set(redo) then skipifexists=0 else skipifexists=1
get_all_utdates,utdate
for i=0,n_elements(utdate)-1 do begin
	if keyword_set(utrange) then begin
		if (strpos(utdate[i],utrange) eq -1) then continue
	endif

	datadir = getenv('BOK90PRIMEDIR')+'/data/'

        reduxdir = getenv('BOK90PRIMEOUTDIR')
        outdir = reduxdir+'/ut'+utdate[i]+'/'
	if not file_test(outdir+'ccdproc3',/dir) then continue
	print, outdir

        log = read_log(utdate[i])
	index=where((log.imtype eq 'object') and (log.reject eq 0))
	log=log[index]

	residual=datadir+['xmap.fits','ymap.fits']
	
	if keyword_set(astrom) then begin
	bok_astrom,log.fn,sdssfn=datadir+'sdss.fits', $
	           ucac4fn=datadir+'ucac4.fits', $
	           scale=8000.,outdir=outdir,/sex,residual=residual, $
	           skipifexists=skipifexists
	endif

	get_default_apers,aperrad,skyrad

	if keyword_set(photom) then begin
	bok_photom,log.fn,sdssfn=datadir+'sdss.fits',aperrad=aperrad,$
        	skyrad=skyrad,outdir=outdir,scale=8000.,/ps,/setslope, $
	        skipifexists=skipifexists
	endif
endfor
end

