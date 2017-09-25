pro dustmap_binary_writer, filename, particle_IDs, x, y, z, vx=vx, vy=vy, vz=vz, intensity=intensity, append=append, silent=silent

;dustmap_binary_writer.pro
;Writes discrete particle coordinates to a file in the
;format required by dustmap.pro

;MANDATORY INPUTS
;filename: a string indicating the name of the output file to create
;particle_IDs: an n-element vector of integers specifying the ID# of each particle
;x, y, z: n-element vectors of 3D coordinates for each particle

;OPTIONAL INPUTS
;vx, vy, vz: n-element vectors of 3D velocities for each particle
;intensity: n-element vector of the "intensity" of each particle


;---------------------
;Make local copies of the input vectors
if n_elements(particle_IDs) gt 0 then uparticle_IDs = particle_IDs
if n_elements(x) gt 0 then ux = x
if n_elements(y) gt 0 then uy = y
if n_elements(z) gt 0 then uz = z
if n_elements(vx) gt 0 then uvx = vx
if n_elements(vy) gt 0 then uvy = vy
if n_elements(vz) gt 0 then uvz = vz
if n_elements(intensity) gt 0 then uintensity = intensity
;---------------------



;---------------------
;Perform some checks
if n_elements(filename) ne 1 then stop, 'ERROR: inputfile must be a single string.'
if (size(filename,/type) ne 7) then stop, 'ERROR: filename must be a string.'
if n_elements(uparticle_IDs) lt 1 then stop, 'ERROR: particle_IDs must be specified.'
if n_elements(ux) lt 1 then stop, 'ERROR: x must be specified.'
if n_elements(uy) lt 1 then stop, 'ERROR: y must be specified.'
if n_elements(uz) lt 1 then stop, 'ERROR: z must be specified.'
nparticles = n_elements(uparticle_IDs)
if n_elements(ux) ne nparticles or n_elements(uy) ne nparticles or n_elements(uz) ne nparticles then $
   stop,'ERROR: particle_IDs, x, y, and z arrays must be equal in length.'
if n_elements(uvx) ne n_elements(uvy) or n_elements(uvx) ne n_elements(uvz) or n_elements(uvy) ne n_elements(uvz) then $
   stop, 'ERROR: vx, vy, and vz must be equal in length.'

if n_elements(uvx) gt 0 then do_v = 1 else do_v = 0
if n_elements(uintensity) gt 0 then do_i = 1 else do_i = 0

datatype = 1
if do_v and not do_i then datatype = 2
if not do_v and do_i then datatype = 3
if do_v and do_i then datatype = 4
if not keyword_set(silent) then print,'Creating a binary file with datatype = '+strcompress(string(datatype),/remove_all)
;---------------------


;---------------------
;Convert the input values to the correct type

;dustmap requires 4 byte int for particle ID#
if size(uparticle_IDs,/type) ne 3 and not keyword_set(silent) then print,'Converted particle_IDs vector to 4 byte integer.'
uparticle_IDs = long(uparticle_IDs)

;dustmap requires 4 byte floats for x,y,z coords
if size(ux,/type) ne 4 and not keyword_set(silent) then print,'Converted x vector to 4 byte floating point.'
if size(uy,/type) ne 4 and not keyword_set(silent) then print,'Converted y vector to 4 byte floating point.'
if size(uz,/type) ne 4 and not keyword_set(silent) then print,'Converted z vector to 4 byte floating point.'
ux = float(ux)
uy = float(uy)
uz = float(uz)

;dustmap requires 4 byte floats for vx,vy,vz
if do_v then begin
   if size(uvx,/type) ne 4 and not keyword_set(silent) then print,'Converted vx vector to 4 byte floating point.'
   if size(uvy,/type) ne 4 and not keyword_set(silent) then print,'Converted vy vector to 4 byte floating point.'
   if size(uvz,/type) ne 4 and not keyword_set(silent) then print,'Converted vz vector to 4 byte floating point.'
   uvx = float(uvx)
   uvy = float(uvy)
   uvz = float(uvz)
endif

;dustmap requires 4 byte int for particle ID#
if do_i then begin
   if size(uintensity,/type) ne 5 and not keyword_set(silent) then print,'Converted intensity vector to 8 byte floating point.'
   uintensity = double(uintensity)
endif
;---------------------


;---------------------
;Open the output file and write it
openw, lun, filename, /get_lun, append=append
i = ulong(0)
while i lt nparticles do begin
   writeu, lun, uparticle_IDs[i], ux[i], uy[i], uz[i]
   if datatype eq 2 or datatype eq 4 then writeu, lun, uvx[i], uvy[i], uvz[i]
   if datatype eq 3 or datatype eq 4 then writeu, lun, uintensity[i]
   i++
endwhile
free_lun, lun
;---------------------


end
