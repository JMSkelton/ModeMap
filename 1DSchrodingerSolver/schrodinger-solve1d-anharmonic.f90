PROGRAM schrodinger1d
!
! This program solves the 1D Schrodinger equation for a given potential.
! The equation can be solved using different approaches - whichever
! subroutine is chosen determines the solution method.
!
! The main algorithm was taken from my PhD work (2010)
!
! J. Buckeridge November 2014
!
!
  implicit none
!
!  integer, parameter :: nmax = 1000
  integer, parameter :: nmax = 2048
  real*8, parameter :: small = 1.d-9
  real*8, parameter :: vlarge = 1.d300
  real*8, parameter :: bhfac = 1.d3
  real*8, parameter :: bwfac = 1.d-1
  real*8, parameter :: convergez = 1.d-6
  real*8, parameter :: pi=dacos(-1.0d0)     ! pi
  real*8, parameter :: q=1.6021765314d-19   ! Electron charge (C)
  real*8, parameter :: hbar=1.05457148d-34  ! Dirac's constant (Js)
  real*8, parameter :: kboltz = 8.617343d-5 ! Boltzmann constant (eV/K)  
  integer i, j, k, ngrid, nval, ntot, stat, numsol, itestz
  integer check_v, numcoeff, xpos(nmax)
  real*8 xval(nmax), v(nmax), x_old(nmax), v_old(nmax)
  real*8 newstep, vslope, temp, coeff(nmax), xmin, xmax
  real*8 vmin, vmax, test
  real*8 eigenvalues(nmax), eigenvectors(nmax,nmax)
  real*8 mass, temperature
  real*8 testz, testz1, partz, omega, omegaconv
  real*8 avgx(nmax), normx, thermx, occup
  character*100 line
  character*1 routine
  logical check_ft
!
! Write header
!
  write(*,'(a)') "******************************************************"
  write(*,'(a)') 
  write(*,'(a)') "      PROGRAM TO SOLVE 1D SCHROEDINGER EQUATION" 
  write(*,'(a)') 
  write(*,'(a)') "Potential assumed to be eV, distance in (a.m.u.)^1/2 A" 
  write(*,'(a)') 
  write(*,'(a)') "------"
  write(*,'(a)') "j.buckeridge@ucl.ac.uk 2015"
  write(*,'(a)') "******************************************************"
  write(*,'(a)') 
!
! Open input file
!
  open (unit=11, FILE='input.dat', status='old', iostat=stat)
  if(stat /= 0) then
     write(*,'(a)') "ERROR: input.dat file not found!!"
     stop
  endif
!
! Read in number of grid points from input file
!
  write(*,'(a)') "Reading in desired number of grid points..."
  do
     read(11,*) line
     if(line(1:1) == "#") then
        cycle
     endif
     read(line,*,iostat=stat) ngrid
     if(stat /= 0) then
        write(*,'(a)') "ERROR: grid point number not found in input file!!"
        goto 1000
     elseif(ngrid > nmax) then
        write(*,'(a)') "ERROR: grid point number is greater than max (2048)!!"
        goto 1000
     else
        exit
     endif
  enddo
  write(*,'(a)') "...done"
!
! Read in mass from input file
!
  write(*,'(a)') "Reading in mass factor..."
  do
     read(11,*) line
     if(line(1:1) == "#") then
        cycle
     endif
     read(line,*,iostat=stat) mass
     if(stat /= 0) then
        write(*,'(a)') "ERROR: mass not found in input file!!"
        goto 1000
     else
        exit
     endif
  enddo
  write(*,'(a)') "...done"
!
! Read in which subroutine is used from input file
!
  check_ft = .false.
  write(*,'(a)') "Reading in subroutine choice..."
  do
     read(11,*) line
     if(line(1:1) == "#") then
        cycle
     endif
     read(line(1:1),'(a1)',iostat=stat) routine
     if(stat /= 0) then
        write(*,'(a)') "ERROR: subroutine choice not found in input file!!"
        goto 1000
     else
        if(routine == "F" .or. routine == "f") then
           check_ft = .true.
           write(*,'(a)') "...done"
           write(*,'(a)') "Will attempt to solve Schroedinger equation in Fourier space:"
        elseif(routine == "S" .or. routine == "s") then
           write(*,'(a)') "...done"
           write(*,'(a)') "Will attempt to solve Schroedinger equation using &
                &shooting method:"
           write(*,'(a)') "(WARNING: This method is only applicable to quantum &
                &well-like potentials -"
           write(*,'(a)') "          solutions above Vmax are ignored)"
        else
           write(*,'(a)') "ERROR: syntax error when reading in subroutine choice!!"
           goto 1000
        endif
        exit
     endif
  enddo
!
! Read in temperature from input file
!
  write(*,'(a)') "Reading in temperature (in Kelvin)..."
  do
     read(11,*) line
     if(line(1:1) == "#") then
        cycle
     endif
     read(line,*,iostat=stat) temperature
     if(stat /= 0) then
        write(*,'(a)') "ERROR: temperature not found in input file!!"
        goto 1000
     else
        exit
     endif
  enddo
  write(*,'(a)') "...done"
!
! Read in whether potential is polynomial fit or x,y data file
!
  check_v = 2
  do
     read(11,*) line
     if(line(1:1) == "#") then
        cycle
     endif
     read(line,*,iostat=stat) check_v
     if(stat /= 0) then
        write(*,'(a)') "ERROR: syntax error when checking potential format!!"
        goto 1000
     else
        exit
     endif
  enddo
!
! x,y data file format - will be linearly interpolated
!
  if(check_v == 0) then
!
! Open potential file
!
     write(*,'(a)') "Opening potential.dat file to read potential..."
     open (unit=12, FILE='potential.dat', status='old', iostat=stat)
     if(stat /= 0) then
        write(*,'(a)') "ERROR: potential file not found!!"
        goto 1000
     endif
!
! Read in potential and check file length
!
     i = 0
     do
        i = i+1
        if(i > ngrid - 4) then    ! 4 extra values to be added for barriers
           write(*,'(a)') "ERROR: potential file has too many entries!!"
           goto 1000
        endif
        read(12,*,end=101) x_old(i), v_old(i)
     enddo
101  ntot = i - 1
     write(*,'(a23,i4,x,a8)') "...read successfully (", ntot, "entries)"
!
! Linearly interpolate the potential so that it has ngrid entries
!    
! First find the step size and store old arrays
!
     write(*,'(a49,x,i5,x,a10)') "Now linearly interpolating data so that there are",&
          &ngrid,"entries..."
     newstep = ( x_old(ntot) - x_old(1) ) / real(ngrid-1)
     v = 0.d0
     xval = 0.d0
!
! Figure out positions in range 1 -> ngrid of x and v values
!
     xpos(1) = 1
     xval(1) = x_old(1)
     v(1) = v_old(1)
     do i=2,ntot
        test = ( x_old(i) - x_old(1) ) / newstep
        xpos(i) = int(test) + 1
        call round_int(test,xpos(i))
        xval(xpos(i)) = x_old(i)
        v(xpos(i)) = v_old(i)
     enddo
!
! Now determine the step sizes between each pair of values (and slope for v)
! and interpolate linearly between them
!
     do i=1,ntot-1
        newstep = ( x_old(i+1) - x_old(i) ) / real(xpos(i+1) - xpos(i))
        vslope = ( v_old(i+1) - v_old(i) ) / ( x_old(i+1) - x_old(i) )
        do j=xpos(i)+1,xpos(i+1)-1
           xval(j) = xval(xpos(i)) + real(j - xpos(i)) * newstep
           v(j) = v(xpos(i)) + real(j - xpos(i)) * newstep * vslope
        enddo
     enddo
!
     write(*,'(a)') "...done"
!
! Polynomial fitted data
!
  elseif(check_v == 1) then
!
! Read in coefficients 
!
     write(*,'(a)') "Polynomial fitted function found..." 
     do
        read(11,*) line
        if(line(1:1) == "#") then
           cycle
        endif
        read(line,*,iostat=stat) numcoeff
        if(stat /= 0) then
           write(*,'(a)') "ERROR: syntax error when reading number of coefficients!!"
           goto 1000
        else
           exit
        endif
     enddo
     write(*,'(a20,x,i4)') "Order of polynomial:", numcoeff - 1
     do i=1,numcoeff
        do
           read(11,*) line
           if(line(1:1) == "#") then
              cycle
           endif
           read(line,*,iostat=stat) coeff(i)
           if(stat /= 0) then
              write(*,'(a)') "ERROR: syntax error when reading coefficients!!"
              goto 1000
           else
              exit
           endif
        enddo
     enddo
     write(*,'(a)') "Coefficients read in successfully..."
!
! Read in min and max x values
!
     write(*,'(a)') "...reading in max and min x values..."
     do
        read(11,*) line
        if(line(1:1) == "#") then
           cycle
        endif
        read(line,*,iostat=stat) xmin
        if(stat /= 0) then
           write(*,'(a)') "ERROR: syntax error when reading xmin!!"
           goto 1000
        else
           exit
        endif
     enddo  
     do
        read(11,*) line
        if(line(1:1) == "#") then
           cycle
        endif
        read(line,*,iostat=stat) xmax
        if(stat /= 0) then
           write(*,'(a)') "ERROR: syntax error when reading xmax!!"
           goto 1000
        else
           exit
        endif
     enddo
!
! Check that they make sense
!
     if(abs(xmax - xmin) < small) then
        write(*,'(a)') "ERROR: xmin and xmax are the same!!"
        goto 1000
     endif
     if(xmin > xmax) then
        temp = xmax
        xmax = xmin
        xmin = temp
     endif
!
! Calculate function
!
     write(*,'(a)') "...done. Now determining function..."
     newstep = ( xmax - xmin ) / real(ngrid-1)
     do i=1,ngrid
        xval(i) = xmin + real(i-1) * newstep
        v(i) = 0.d0
        do j=1,numcoeff
           v(i) = v(i) + coeff(j) * xval(i)**(j-1)
        enddo
     enddo
     write(*,'(a)') "...done"
!
! Otherwise there has been an error in reading in potential data
!
  else
     write(*,'(a)') "ERROR: failure to read in correctly status of potential data!!"
     goto 1000
  endif
!
! Find min and max values of potential
!
  vmax = v(1)
  vmin = v(1)
  do i=1,ngrid
     if(v(i) > vmax) vmax = v(i)
     if(v(i) < vmin) vmin = v(i)
  enddo
  write(*,'(a27,x,e21.13,x,a2)') "Minimum value of potential:", vmin, "eV"
  write(*,'(a27,x,e21.13,x,a2)') "Maximum value of potential:", vmax, "eV"
!
! Write interpolated potential to file
!
  write(*,'(a)') "Writing potential to file potential-new.dat..."
  open (unit=13, FILE='potential-new.dat', status='replace', iostat=stat)
  do i=1,ngrid
     write(13,'(2(e21.13,2x))') xval(i), v(i)
  enddo
  write(*,'(a)') "...done"
!
! Call routine to solve schrodinger equation
!
  numsol = ngrid
  if(check_ft) then
     call schro_solve_ft(xval, v, mass, eigenvalues(1:numsol), &
          &eigenvectors(1:numsol,1:numsol), numsol)
  else
     call schro_solve_sm(xval, v, mass, eigenvalues(1:numsol), &
          &eigenvectors(1:numsol,1:numsol), numsol)
  endif
!
! Write output, first writing eigenvalues
!
  open (unit=14, file='eigenvalue.dat', status='replace', iostat=stat)
  do i=1,numsol
     write(14,'(e21.13)') eigenvalues(i)
  enddo
!
! Now eigenvectors, shifting each by its corresponding eigenvalue (for nice
! plots!)
!
  open (unit=15, file='eigenvector.dat', status='replace', iostat=stat)
  do i=1,numsol
     do j=1,ngrid
        write(15,'(2(e21.13,x))') xval(j), eigenvectors(j,i) + eigenvalues(i)
     enddo
     write(15,*)
  enddo
!
! Now determine the mode frequency. The procedure is:
!
! First determine partition function Z = sum_i exp(-E_i/kT),
! where E_i are the eigenvalues from Schroedinger solver.
!
! Then assume the oscillator partition function is the same as 
! this, i.e. Z = sum_n exp(-(n + 1/2) * hbar * omega / kT)
! and get omega from this assumption.
!
! i.e. omega = (2kT / hbar) * arcsinh[ ( 1/2Z ) ]
!
!
! First check how many eigenvalues need to be included to let the partition 
! function sum converge to 10^-6
!
  write(*,'(a)')
  write(*,'(a)') "Determining corresponding mode frequency..."
!
! Calculate the first term in the partition function. Convergence will be 
! tested against this value
!
  if(temperature == 0.d0) then
     testz1 = 0.d0
  else
     testz1 = exp(-(eigenvalues(1) - vmin) / (kboltz * temperature))
  endif
  i = 1
  do
     i = i+1
     if(temperature == 0.d0) then
        testz = 0.d0
     else
        testz = exp(-(eigenvalues(i) - vmin) / (kboltz * temperature)) / testz1
     endif
     if(i > numsol) then
        write(*,'(a)') "WARNING: Partition function not converged using all eigenvalues -"
        write(*,'(a30,x,e21.13)') "       : Convergence in sum is", testz
        itestz = numsol
        exit
     endif
     if(testz <= convergez) then
        itestz = i - 1
        write(*,'(a4,x,i4,x,a64)') "Used", itestz, "eigenvalues to get convergence of 10^-6 &
             &in partition function..."
        exit
     endif
  enddo
!
! Calculate Z
!
  partz = 0.d0
  if(temperature /= 0.d0) then
     if(itestz > 0) then
        do i=1,itestz
           partz = partz + exp(-(eigenvalues(i) - vmin) / (kboltz * temperature))
        enddo
     endif
  endif
!
! Write the eigenvectors associated with eigenvalues used in summation for
! Z to a separate file for plotting
!
  if(itestz > 0 .and. itestz < numsol) then
     open (unit=16, file='wavefns.dat', status='replace', iostat=stat)
     write(*,'(a)') "Writing eigenvectors used in determining partition function to file"
     write(*,'(a)') "wavefns.dat for plotting..."
     do i=1,itestz
        do j=1,ngrid
           write(16,'(2(e21.13,x))') xval(j), eigenvectors(j,i) + eigenvalues(i)
        enddo
        write(16,*)
     enddo
  endif
!
  if(partz < 1.d-15) then
     write(*,'(a)')
     write(*,'(a)') "WARNING: partition function is very close to zero!!"
     write(*,'(a)') "       : (Not enough eigenvalues < kT)"
     write(*,'(a)')
     !goto 1000
  endif
!
! Now calculate omega in rad/s, convert to THz, and print result
!
! For T = 0 K, omega = 2 * E_1 / hbar (see labbook 5/1/16 - JB)
!
  if(temperature == 0.d0) then
     omega = 2.d0 * (eigenvalues(1) - vmin) / (hbar/q)
  elseif( (1.d0 / (2.d0 * partz)) > vlarge) then
     omega = (2.d0 * kboltz * temperature / (hbar/q)) * asinh(vlarge)
  else
     omega = (2.d0 * kboltz * temperature / (hbar/q)) * asinh(1.d0 / (2.d0 * partz))
  endif
  omegaconv = omega / (2.d0 * pi * 1.d12)
  write(*,'(a)') 
  write(*,'(a7,x,e13.5,x,a5)') "Omega =", omega, "rad/s"
  write(*,'(a7,x,f13.5,x,a3)') "Omega =", omegaconv, "THz"
  omegaconv = (omega / (2.d0 * pi)) * 4.13558d-12
  write(*,'(a7,x,f13.5,x,a3)') "Omega =", omegaconv, "meV"
  omegaconv = (omega / (2.d0 * pi)) * 3.33565d-11
  write(*,'(a7,x,f13.5,x,a4)') "Omega =", omegaconv, "cm-1"
  write(*,'(a)') 
!
! Calculate average position for each eigenmode from which the thermally
! averaged position can be determined
!
!  x_i  = < psi_i | x | psi_i > / < psi_i | psi_i>
!
! < x > = sum_i x_i exp( -E_i/kT ) / Z
!
  write(*,'(a)') "Calculating thermally averaged x position..."
  open (unit=17, file='avgxpos.dat', status='replace', iostat=stat)
!
! First get the x_i. Watch out for spurious eigenvectors with zero norm
! (may be a problem when using Fourier routine)
!
  do i=1,itestz
     avgx(i) = 0.d0
     normx = 0.d0
     do j=1,ngrid
        avgx(i) = avgx(i) + xval(j) * eigenvectors(j,i)**2
        normx = normx + eigenvectors(j,i)**2
     enddo
     if(abs(normx) < small) then
        avgx(i) = 0.d0
     else
        avgx(i) = avgx(i) / normx
     endif
     write(17,'(e21.13)') avgx(i)
  enddo
!
! Now get < x >
!
  thermx = 0.d0
  if(temperature == 0.d0) then
     thermx = 0.d0
  else
     do i=1,itestz
        thermx = thermx + avgx(i) * exp(-(eigenvalues(i) - vmin) / &
             &(kboltz * temperature))
     enddo
     thermx = thermx / partz
  endif
!
! Print result
!
  write(*,'(a)')
  write(*,'(a19,x,e13.5,x,a9)') "Thermal average x =", thermx, "amu^1/2 A"
  write(*,'(a)') 
!
! Calculate thermal average mode occupation
!
! < s > = 1 / exp[hbar omega / kT] - 1
!
  write(*,'(a)') "Calculating thermal occupation of mode..."
!
  if(temperature == 0.d0) then
     occup = 0.d0
  else
     occup = 1.d0 / ( exp(( hbar * omega / q) / ( kboltz * temperature )) - 1.d0)
  endif
!
! Print result
!
  write(*,'(a)')
  write(*,'(a28,x,e13.5)') "Thermal occupation of mode =", occup
  write(*,'(a)') 
!
1000 continue
!
! Close files
!
  do i=11,17
     close(unit=i)
  enddo
!
end program schrodinger1d
! **************************************************************************
!
subroutine schro_solve_sm(x, v, m, eigval, eigvec, num)
!
! Subroutine to solve the 1D Schrodinger equation, given a 1D potential
!
implicit none
!
  integer, parameter :: nmax = 2048
!  integer, parameter :: nmax = 1000
  integer, parameter :: short = 5
  real*8, parameter :: small = 1.d-5
  real*8, parameter :: coarse_small = 1.d-5
  real*8, parameter :: vsmall = 1.d-12
  real*8, parameter :: blowup = 1.d24
  real*8, parameter :: stepsize = 1.d3
  real*8, parameter :: pi=dacos(-1.0D0)    ! pi
  real*8, parameter :: m_e=9.10938188D-31  ! Electron mass (kg)
  real*8, parameter :: m_p=1.66053892d-27  ! Atomic mass scale (kg)
  real*8, parameter :: k_B=8.617343d-5     ! Boltzmann constant (eV/K)
  real*8, parameter :: q=1.6021765314D-19  ! Electron charge (C)
  real*8, parameter :: hbar=1.05457148D-34 ! Dirac's constant (Js)
  real*8, parameter :: eps=8.8541872D-21   ! Permittivity of free space (F/nm)
!
  integer i, j, k, num, ival1_fe, imax(num), max_subbands
  real*8 x(num), v(num), eigval(num), eigvec(num,num), m
  real*8 eff_mass(num), energy_sum, vmin, vmax, dz
  real*8  energy_fe, diff1_fe, energy_old, energy_old1
  real*8 energy_new, temp1, scanlength
!
! Initialise effective mass tensor
!
  eff_mass = m * m_p
!
! Find max and min values of potential
!
  vmin = v(1)
  vmax = v(1)
  do i=1,num
     if( v(i) < vmin ) vmin = v(i)
     if( v(i) > vmax ) vmax = v(i)
  enddo
  scanlength = 1.d0 * (vmax - vmin)
!
! Set step size along z
!
  dz = ( x(num) - x(1) ) / real(num-1)
!
! Max solutions to search for. Initialise as num
!
  max_subbands = num
!
! First do a coarse step through the energy range, finding the first max_subbands
! solutions by determining where the divergence of the wavefn changes sign
!
  write(*,'(a)') "Scanning energy range to find solutions via shooting method..."
  i = 0
  j = 1
  ival1_fe = 0
  energy_sum = coarse_small * (vmax - vmin)
  do
     i = i+1
     energy_fe = vmin + energy_sum * real(i)
     if(energy_fe >= vmin + scanlength) then
        if(j == 0) then
           write(*,'(a)') "ERROR: No solutions found!!"
           write(*,'(a)') "Are you sure the shooting method is appropriate for this potential?"
           stop
        endif
        write(*,'(a5,x,i4,x,a37)') "Found", j-1, "solutions in coarse scan. Refining..."
        max_subbands = j-1
        exit
     endif
     eigvec(:,j) = 0.d0
     eigvec(1,j) = 0.d0
     eigvec(2,j) = 1.d0
!
! Use equation 3.53 of "Quantum wells, wires and dots" by Harrison (2005 Wiley) 
! to calculate the wave function at energy_fe. Effective mass are in kg, keep 
! track of units. If the wavefn blows up (>10^10) don't calculate further
!
     imax(j) = num
     do k = 3,num
        eigvec(k,j) = (eff_mass(k)+eff_mass(k-1))/2.d0 * &
             &( ( ( 2.d0 * dz**2 * q * 1.d-20 ) * ( v(k-1) - energy_fe)/hbar**2 + 2.d0 / &
             &(eff_mass(k)+eff_mass(k-1)) + 2.d0 / (eff_mass(k-1)+eff_mass(k-2)) )&
             &* eigvec(k-1,j) - eigvec(k-2,j) * 2.d0 / (eff_mass(k-1)+eff_mass(k-2)) )
        if(abs(eigvec(k,j)) >= blowup) then
           imax(j) = k
           exit
        endif
     enddo
!
! Determine how the wave functions diverge at infinity (far to the right). Where the 
! direction of divergence changes, take the midpoint energy as the first approximation 
! to the eigenvalue. If the wave functions are flat and close to zero take the energy 
! as a first approximation to the eigenvalue
!
     diff1_fe = eigvec(imax(j),j) - eigvec(imax(j)-short,j)
     if (abs(diff1_fe) <= small .and. abs(eigvec(imax(j),j) - eigvec(1,j)) <= small) then
        eigval(j) = energy_fe
        if(j == max_subbands) exit
        j = j+1
        ival1_fe = 0
     elseif(diff1_fe > 0.d0) then
        if(ival1_fe == 0) ival1_fe = 1
        if(ival1_fe == -1) then
           eigval(j) = (energy_fe + vmin + energy_sum * real(i-1)) / 2.d0
           if(j == max_subbands) exit
           j = j+1
           ival1_fe = 0
        endif
     elseif(diff1_fe < 0.d0) then
        if(ival1_fe == 0) ival1_fe = -1
        if(ival1_fe == 1) then
           eigval(j) = (energy_fe + vmin + energy_sum * real(i-1)) / 2.d0
           if(j == max_subbands) exit
           j = j+1
           ival1_fe = 0
        endif
     endif
  enddo
!
! Now refine solutions by stepping through energy values on an increasingly fine grid
!
  do j=1,max_subbands
     energy_sum = coarse_small * (vmax - vmin)
     do
        energy_old = eigval(j) - 0.5d0 * energy_sum
        energy_old1 = eigval(j) + 0.5d0 * energy_sum
        i = 0
        ival1_fe = 0
        do
           i = i+1
           energy_fe = energy_old + real(i-1) * (energy_old1 - energy_old) / stepsize
           if(energy_fe > energy_old1) then
              write(*,'(a23,x,i4)') "ERROR refining solution", j
              goto 108
           endif
           imax(j) = num
           eigvec(1,j) = 0.d0
           eigvec(2,j) = 1.d0
           do k = 3,num
              eigvec(k,j) = (eff_mass(k)+eff_mass(k-1))/2.d0 * ( ( &
                   &( 2.d0 * dz**2 * q * 1.d-20 ) * ( v(k-1) - energy_fe)/hbar**2 + 2.d0 / &
                   &(eff_mass(k)+eff_mass(k-1)) + 2.d0 / (eff_mass(k-1)+eff_mass(k-2)) )&
                   &* eigvec(k-1,j) - eigvec(k-2,j) * 2.d0 / (eff_mass(k-1)+eff_mass(k-2)) )
              if(abs(eigvec(k,j)) >= blowup) then
                 imax(j) = k
                 exit
              endif
           enddo
           diff1_fe = eigvec(imax(j),j) - eigvec(imax(j)-short,j)
           if (abs(diff1_fe) <= small .and. abs(eigvec(imax(j),j) - eigvec(1,j)) <= small) then
              eigval(j) = energy_fe
              ival1_fe = 0
              goto 108
           elseif(diff1_fe > 0.d0) then
              if(ival1_fe == 0) ival1_fe = 1
              if(ival1_fe == -1) then
                 eigval(j) = energy_fe - ((energy_old1 - energy_old) / stepsize) / 2.d0
                 goto 108
              endif
           elseif(diff1_fe < 0.d0) then
              if(ival1_fe == 0) ival1_fe = -1
              if(ival1_fe == 1) then
                 eigval(j) = energy_fe - ((energy_old1 - energy_old) / 1000.d0) / 2.d0
                 goto 108
              endif
           endif
        enddo
108     if(abs(eigvec(imax(j),j)) <= small .and. imax(j) == num) exit
        energy_sum = energy_sum / stepsize
        if(energy_sum < vsmall) exit
     enddo
  enddo
  write(*,'(a)') "...done"
  write(*,'(a)') "Now renormalising eigenvectors..."
!
! Now re-normalise eigenvectors. Watch out for any that 'blew up'
!
  do i=1,max_subbands
     energy_sum = 0.d0
     do j=1,imax(i)
        energy_sum = energy_sum + eigvec(j,i)**2
     enddo
     eigvec(:,i) = eigvec(:,i) / sqrt(energy_sum)
  enddo
!
! Wavefunctions are now normalized over the grid points, must re-normalize
! so that they are normalized over the length of the potential
! i.e. so that INT |psi(z)|^2 dz = 1 
! 
  eigvec = eigvec / sqrt(dz)
!
  write(*,'(a)') "...done"
  write(*,'(a)') "WARNING: highest energy eigenvector may not be accurate with this approach"
!
! Output number of solutions
!
  num = max_subbands
!
end subroutine schro_solve_sm
! **************************************************************************
!
subroutine schro_solve_ft(x, v, m, eigval, eigvec, num)
!
! Subroutine to solve the 1D Schrodinger equation, given a 1D potential
!
implicit none
!
  integer, parameter :: nmax = 2048
  integer, parameter :: nmaxf = 4096
!  integer, parameter :: nmax = 1000
!  integer, parameter :: nmaxf = 2048
  integer, parameter :: large = 1e8
  real*8, parameter :: small = 1.d-9
  real*8, parameter :: pi=dacos(-1.0D0)    ! pi
  real*8, parameter :: m_e=9.10938188D-31  ! Electron mass (kg)
  real*8, parameter :: m_p=1.66053892d-27  ! Atomic mass scale (kg)
  real*8, parameter :: k_B=8.617343d-5     ! Boltzmann constant (eV/K)
  real*8, parameter :: q=1.6021765314D-19  ! Electron charge (C)
  real*8, parameter :: hbar=1.05457148D-34 ! Dirac's constant (Js)
  real*8, parameter :: eps=8.8541872D-22   ! Permittivity of free space (F/A)
!
  integer i, j, k, num, numf, ierr
  real*8 x(num), v(num), eigval(num), eigvec(num,num), m
  real*8 vtmp(nmaxf), vf(2*nmaxf), eigvaltmp(nmaxf)
  real*8 eigr(nmaxf,nmaxf), eigi(nmaxf,nmaxf)
  real*8 hamr(nmaxf,nmaxf), hami(nmaxf,nmaxf)
  real*8 fv1(nmaxf), fv2(nmaxf), fm1(2,nmaxf)
  real*8 feigv(2*nmaxf,nmaxf)
  real*8 dz, kterm
!
! Get array sizes and number of grid steps for fourier routine
!
  i = 0
  do
     i = i+1
     numf = 2**i
     if(numf > int(real(num) * 1.2d0)) exit
     if(i > large .or. numf > nmaxf) then
        write(*,'(a)') "ERROR: setting fourier grid failed - too many points!!"
        stop
     endif
  enddo
  write(*,'(a16,x,i4,x,a33)') "Fourier grid has", numf, "points. &
       &Transforming potential..."
!
! Set effective mass
!
  m = m * m_p
!
! Set step size along z
!
  dz = ( x(num) - x(1) ) / real(num-1)
!
! Set up potential for fourier transform
!
  do i=1,num
     vtmp(i) = v(i)
  enddo
  do i=num+1,numf
     vtmp(i) = v(num)
  enddo
!
  vf = 0.d0
  do i=1,numf
     vf(2*i - 1) = vtmp(i)
  enddo
!
! Fourier transform vf using subroutine dfour1, and scale by 1/numf as 
! is required by this routine
!
  call dfour1 (vf(1:2*numf),numf,1)
!
  vf = vf * ( 1.d0 / numf )
  write(*,'(a)') "...done"
  write(*,'(a)') "Now solving Schroedinger equation in Fourier space..."
!
! Reinitialize variables   
!
  hamr = 0.d0
  hami = 0.d0
  eigr = 0.d0
  eigi = 0.d0
  eigvaltmp = 0.d0
  feigv = 0.d0      
!
! Construct the Hamiltonian for the system, which consists of a plane-wave
! part (numf * numf matrix) constructed from the fourier transform of the 
! potential vf, and a kinetic energy term which is added to the diagonal 
! elements of the plane wave part
!     
! First the plane-wave part. This block constructs an nxn matrix hamr that 
! consists of real parts of the fourier components of the potential in decreas-
! ing order along the rows and increasing order down the columns
!
  do i=1,numf
     do j=-i+numf+2,-i+2*numf+1
        k=j
        if (k > numf) then
           k = k-numf 
        endif
        hamr(j-numf+i-1,i) = vf(2*k-1)
     enddo
  enddo
!
! This block constructs an nxn matrix consisting of the imaginary parts of the 
! fourier components of the potential (same order as hamr)
!
  do i=1,numf
     do j=-i+numf+2,-i+2*numf+1 
        k=j
        if (k > numf) then
           k = k-numf
        endif
        hami(j-numf+i-1,i) = vf(2*k)
     enddo
  enddo
!
! Add on kinetic energy contribution (with k=2*pi*n/L, n=0,+/-1,..)
! by first calculating (hbar^2)/(2m*) in eV(A^2), multiplying by
! k^2 in 1/(A^2), and adding to each diagonal component of matrix
!     
! To get (hbar^2)/(2m*) in eV(A^2):
!    
  kterm=( (( (hbar / q)**2 ) * q )/ ( 2.d0 * m ) ) * 1.d20
!
! Adding to diagonal elements of matrix
!
  do i=1,numf
     hamr(i,i) = hamr(i,i) + kterm * &
          &( 2.d0 * pi * ( i - ( numf/2 + 1 ) ) / ( numf * dz ) )**2
  enddo
!
! Diagonalize Hamiltonian using eispack subroutine ch, to get the eigenve-
! ctors and eigenvalues
!
  call ch (numf,numf,hamr(1:numf,1:numf),hami(1:numf,1:numf),&
       eigvaltmp(1:numf),1,eigr(1:numf,1:numf),eigi(1:numf,1:numf),&
       fv1(1:numf),fv2(1:numf),fm1(1:2,1:numf),ierr)  
!
! Rearrange eigenvectors by placing real and imaginary
! parts of each component one after the other into array feigv, which is 
! an array of length 2*numf (so that routine dfour1 can inverse fourier 
! transform it)
!     
! Arranging real parts first (the elements are assigned from the middle of
! the k-space eigenvectors)
!
  do i=1,numf
     do j=1,numf/2                   
        feigv(2*j-1,i) = eigr(j+numf/2,i)
     enddo
  enddo
!
  do i=1,numf
     do j=numf/2+1,numf
        feigv(2*j-1,i) = eigr(j-numf/2,i)
     enddo
  enddo
!
! Now the imaginary parts
!       
  do i=1,numf
     do j=1,numf/2
        feigv(2*j,i) = eigi(j+numf/2,i)
     enddo
  enddo
!       
  do i=1,numf
     do j=numf/2+1,numf
        feigv(2*j,i) = eigi(j-numf/2,i)
     enddo
  enddo
  write(*,'(a)') "...done"
  write(*,'(a)') "Now inverse transforming and normalising eigenvectors..."
!
! Use subroutine dfour1 to calculate the inverse fourier transform of 
! the first eigenvector, and rescale as is required by the routine
!       
  do i=1,numf
     call dfour1(feigv(1:2*numf,i),numf,-1)
  enddo
!
  feigv(1:2*numf,:) = dsqrt( 1.d0 / numf ) * feigv(1:2*numf,:)
!
! Wavefunctions are now normalized over the grid points, must re-normalize
! so that they are normalized over the length of the potential
! i.e. so that INT |psi(z)|^2 dz = 1 
! 
  feigv = feigv / sqrt(dz)
!
! Return the first nmax solutions, returning just the real parts of 
! the eigenvectors
!
  do i=1,num
     eigval(i) = eigvaltmp(i)
  enddo
  do i=1,num
     do j=1,num
        eigvec(j,i) = feigv(2*j-1,i)
     enddo
  enddo
  write(*,'(a)') "...done"
!
end subroutine schro_solve_ft
! **************************************************************************
!
!
subroutine ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
!
  integer i,j,n,nm,ierr,matz
  double precision ar(nm,n),ai(nm,n),w(n),zr(nm,n),zi(nm,n),&
       fv1(n),fv2(n),fm1(2,n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a complex hermitian matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a=(ar,ai).
!
!        ar  and  ai  contain the real and imaginary parts,
!        respectively, of the complex hermitian matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        zr  and  zi  contain the real and imaginary parts,
!        respectively, of the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1, fv2, and  fm1  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
  if (n .le. nm) go to 10
  ierr = 10 * n
  go to 50
!
10 call  htridi(nm,n,ar,ai,w,fv1,fv2,fm1)
  if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
  call  tqlrat(n,w,fv2,ierr)
  go to 50
!     .......... find both eigenvalues and eigenvectors ..........
20 do i = 1, n
!
     do j = 1, n
        zr(j,i) = 0.0d0
     enddo
!
     zr(i,i) = 1.0d0
  enddo
!
  call  tql2(nm,n,w,fv1,zr,ierr)
  if (ierr .ne. 0) go to 50
  call  htribk(nm,n,ar,ai,fm1,n,zr,zi)
50 return
end subroutine ch
SUBROUTINE dfour1(data,nn,isign)
  INTEGER isign,nn
  DOUBLE PRECISION data(2*nn)
  INTEGER i,istep,j,m,mmax,n
  DOUBLE PRECISION tempi,tempr
  DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
  n=2*nn
  j=1
  do i=1,n,2
     if(j.gt.i)then
        tempr=data(j)
        tempi=data(j+1)
        data(j)=data(i)
        data(j+1)=data(i+1)
        data(i)=tempr
        data(i+1)=tempi
     endif
     m=n/2
1    if ((m.ge.2).and.(j.gt.m)) then
        j=j-m
        m=m/2
        goto 1
     endif
     j=j+m
  enddo
  mmax=2
2 if (n.gt.mmax) then
     istep=2*mmax
     theta=6.28318530717959d0/(isign*mmax)
     wpr=-2.d0*sin(0.5d0*theta)**2
     wpi=sin(theta)
     wr=1.d0
     wi=0.d0
     do m=1,mmax,2
        do i=m,n,istep
           j=i+mmax
           tempr=wr*data(j)-wi*data(j+1)
           tempi=wr*data(j+1)+wi*data(j)
           data(j)=data(i)-tempr
           data(j+1)=data(i+1)-tempi
           data(i)=data(i)+tempr
           data(i+1)=data(i+1)+tempi
        enddo
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
     enddo
     mmax=istep
     goto 2
  endif
  return
end subroutine dfour1
double precision function epslon (x)
  double precision x
!
!     estimate unit roundoff in quantities of size x.
!
  double precision a,b,c,eps
!
!     this program should function properly on all systems
!     satisfying the following two assumptions,
!        1.  the base used in representing floating point
!            numbers is not a power of three.
!        2.  the quantity  a  in statement 10 is represented to 
!            the accuracy used in floating point variables
!            that are stored in memory.
!     the statement number 10 and the go to 10 are intended to
!     force optimizing compilers to generate code satisfying 
!     assumption 2.
!     under these assumptions, it should be true that,
!            a  is not exactly equal to four-thirds,
!            b  has a zero for its last bit or digit,
!            c  is not exactly equal to one,
!            eps  measures the separation of 1.0 from
!                 the next larger floating point number.
!     the developers of eispack would appreciate being informed
!     about any systems where these assumptions do not hold.
!
!     this version dated 4/6/83.
!
  a = 4.0d0/3.0d0
10 b = a - 1.0d0
  c = b + b + b
  eps = dabs(c-1.0d0)
  if (eps .eq. 0.0d0) go to 10
  epslon = eps*dabs(x)
  return
end function epslon
subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)
!
  integer i,j,k,l,m,n,nm
  double precision ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
  double precision h,s,si
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure trbak1, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine forms the eigenvectors of a complex hermitian
!     matrix by back transforming those of the corresponding
!     real symmetric tridiagonal matrix determined by  htridi.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain information about the unitary trans-
!          formations used in the reduction by  htridi  in their
!          full lower triangles except for the diagonal of ar.
!
!        tau contains further information about the transformations.
!
!        m is the number of eigenvectors to be back transformed.
!
!        zr contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        zr and zi contain the real and imaginary parts,
!          respectively, of the transformed eigenvectors
!          in their first m columns.
!
!     note that the last component of each returned vector
!     is real and that vector euclidean norms are preserved.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
  if (m .eq. 0) go to 200
!     .......... transform the eigenvectors of the real symmetric
!                tridiagonal matrix to those of the hermitian
!                tridiagonal matrix. ..........
  do k = 1, n
!
     do j = 1, m
        zi(k,j) = -zr(k,j) * tau(2,k)
        zr(k,j) = zr(k,j) * tau(1,k)
     enddo
  enddo
!
  if (n .eq. 1) go to 200
!     .......... recover and apply the householder matrices ..........
  do i = 2, n
     l = i - 1
     h = ai(i,i)
     if (h .eq. 0.0d0) cycle
!
     do j = 1, m
        s = 0.0d0
        si = 0.0d0
!
        do k = 1, l
           s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
           si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
        enddo
!     .......... double divisions avoid possible underflow ..........
        s = (s / h) / h
        si = (si / h) / h
!
        do k = 1, l
           zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
           zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
        enddo
!
     enddo
!
  enddo
!
200 return
end subroutine htribk
subroutine htridi(nm,n,ar,ai,d,e,e2,tau)
!
  integer i,j,k,l,n,ii,nm,jp1
  double precision ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
  double precision f,g,h,fi,gi,hh,si,scale,pythag
!
!     this subroutine is a translation of a complex analogue of
!     the algol procedure tred1, num. math. 11, 181-195(1968)
!     by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a complex hermitian matrix
!     to a real symmetric tridiagonal matrix using
!     unitary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        ar and ai contain the real and imaginary parts,
!          respectively, of the complex hermitian input matrix.
!          only the lower triangle of the matrix need be supplied.
!
!     on output
!
!        ar and ai contain information about the unitary trans-
!          formations used in the reduction in their full lower
!          triangles.  their strict upper triangles and the
!          diagonal of ar are unaltered.
!
!        d contains the diagonal elements of the the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!        tau contains further information about the transformations.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
  tau(1,n) = 1.0d0
  tau(2,n) = 0.0d0
!
  do i = 1, n
     d(i) = ar(i,i)
  enddo
!     .......... for i=n step -1 until 1 do -- ..........
  do ii = 1, n
     i = n + 1 - ii
     l = i - 1
     h = 0.0d0
     scale = 0.0d0
     if (l .lt. 1) go to 130
!     .......... scale row (algol tol then not needed) ..........
     do k = 1, l
        scale = scale + dabs(ar(i,k)) + dabs(ai(i,k))
     enddo
!
     if (scale .ne. 0.0d0) go to 140
     tau(1,l) = 1.0d0
     tau(2,l) = 0.0d0
130  e(i) = 0.0d0
     e2(i) = 0.0d0
     go to 290
!
140  do k = 1, l
        ar(i,k) = ar(i,k) / scale
        ai(i,k) = ai(i,k) / scale
        h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
     enddo
!
     e2(i) = scale * scale * h
     g = dsqrt(h)
     e(i) = scale * g
     f = pythag(ar(i,l),ai(i,l))
!     .......... form next diagonal element of matrix t ..........
     if (f .eq. 0.0d0) go to 160
     tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
     si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
     h = h + f * g
     g = 1.0d0 + g / f
     ar(i,l) = g * ar(i,l)
     ai(i,l) = g * ai(i,l)
     if (l .eq. 1) go to 270
     go to 170
160  tau(1,l) = -tau(1,i)
     si = tau(2,i)
     ar(i,l) = g
170  f = 0.0d0
!
     do j = 1, l
        g = 0.0d0
        gi = 0.0d0
!     .......... form element of a*u ..........
        do k = 1, j
           g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
           gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
        enddo
!
        jp1 = j + 1
        if (l .lt. jp1) go to 220
!
        do k = jp1, l
           g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
           gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
        enddo
!     .......... form element of p ..........
220     e(j) = g / h
        tau(2,j) = gi / h
        f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
     enddo
!
     hh = f / (h + h)
!     .......... form reduced a ..........
     do j = 1, l
        f = ar(i,j)
        g = e(j) - hh * f
        e(j) = g
        fi = -ai(i,j)
        gi = tau(2,j) - hh * fi
        tau(2,j) = -gi
!
        do k = 1, j
           ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k)&
                + fi * tau(2,k) + gi * ai(i,k)
           ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k)&
                - fi * e(k) - gi * ar(i,k)
        enddo
     enddo
!
270  do k = 1, l
        ar(i,k) = scale * ar(i,k)
        ai(i,k) = scale * ai(i,k)
     enddo
!
     tau(2,l) = -si
290  hh = d(i)
     d(i) = ar(i,i)
     ar(i,i) = hh
     ai(i,i) = scale * dsqrt(h)
  enddo
!
  return
end subroutine htridi
subroutine tqlrat(n,d,e2,ierr)
!
  integer i,j,l,m,n,ii,l1,mml,ierr
  double precision d(n),e2(n)
  double precision b,c,f,g,h,p,r,s,t,epslon,pythag
!
!     this subroutine is a translation of the algol procedure tqlrat,
!     algorithm 464, comm. acm 16, 689(1973) by reinsch.
!
!     this subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the rational ql method.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e2 contains the squares of the subdiagonal elements of the
!          input matrix in its last n-1 positions.  e2(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        e2 has been destroyed.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
  ierr = 0
  if (n .eq. 1) go to 1001
!                                                             
  do i = 2, n
     e2(i-1) = e2(i)
  enddo
!
  f = 0.0d0
  t = 0.0d0
  e2(n) = 0.0d0
!                                                                
  do l = 1, n
     j = 0
     h = dabs(d(l)) + dsqrt(e2(l))
     if (t .gt. h) go to 105
     t = h
     b = epslon(t)
     c = b * b
!     .......... look for small squared sub-diagonal element ..........
105  do m = l, n
        if (e2(m) .le. c) go to 120
!     .......... e2(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
     enddo
!
120  if (m .eq. l) go to 210
130  if (j .eq. 30) go to 1000
     j = j + 1
!     .......... form shift ..........
     l1 = l + 1
     s = dsqrt(e2(l))
     g = d(l)
     p = (d(l1) - g) / (2.0d0 * s)
     r = pythag(p,1.0d0)
     d(l) = s / (p + dsign(r,p))
     h = g - d(l)
!
     do i = l1, n
        d(i) = d(i) - h
     enddo
!
     f = f + h
!     .......... rational ql transformation ..........
     g = d(m)
     if (g .eq. 0.0d0) g = b
     h = g
     s = 0.0d0
     mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
     do ii = 1, mml
        i = m - ii
        p = g * h
        r = p + e2(i)
        e2(i+1) = s * r
        s = e2(i) / r
        d(i+1) = h + s * (h + d(i))
        g = d(i) - e2(i) / g
        if (g .eq. 0.0d0) g = b
        h = g * p / r
     enddo
!
     e2(l) = s * g
     d(l) = h
!     .......... guard against underflow in convergence test ..........
     if (h .eq. 0.0d0) go to 210
     if (dabs(e2(l)) .le. dabs(c/h)) go to 210
     e2(l) = h * e2(l)
     if (e2(l) .ne. 0.0d0) go to 130
210  p = d(l) + f
!     .......... order eigenvalues ..........
     if (l .eq. 1) go to 250
!     .......... for i=l step -1 until 2 do -- ..........
     do ii = 2, l
        i = l + 2 - ii
        if (p .ge. d(i-1)) go to 270
        d(i) = d(i-1)
     enddo
!
250  i = 1
270  d(i) = p
  enddo
!
  go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 return
end subroutine tqlrat
double precision function pythag(a,b)
  double precision a,b
!
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!
  double precision p,r,s,t,u
  p = dmax1(dabs(a),dabs(b))
  if (p .eq. 0.0d0) go to 20
  r = (dmin1(dabs(a),dabs(b))/p)**2
10 continue
  t = 4.0d0 + r
  if (t .eq. 4.0d0) go to 20
  s = r/t
  u = 1.0d0 + 2.0d0*s
  p = u*p
  r = (s/u)**2 * r
  go to 10
20 pythag = p
  return
end function pythag
subroutine tql2(nm,n,d,e,z,ierr)
!
  integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
  double precision d(n),e(n),z(nm,n)
  double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
  ierr = 0
  if (n .eq. 1) go to 1001
!                                                                              
  do i = 2, n
     e(i-1) = e(i)
  enddo
!
     f = 0.0d0
     tst1 = 0.0d0
     e(n) = 0.0d0
!
     do l = 1, n
        j = 0
        h = dabs(d(l)) + dabs(e(l))
        if (tst1 .lt. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
        do m = l, n
           tst2 = tst1 + dabs(e(m))
           if (tst2 .eq. tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
        enddo
!
120     if (m .eq. l) go to 220
130     if (j .eq. 30) go to 1000
        j = j + 1
!     .......... form shift ..........
        l1 = l + 1
        l2 = l1 + 1
        g = d(l)
        p = (d(l1) - g) / (2.0d0 * e(l))
        r = pythag(p,1.0d0)
        d(l) = e(l) / (p + dsign(r,p))
        d(l1) = e(l) * (p + dsign(r,p))
        dl1 = d(l1)
        h = g - d(l)
        if (l2 .gt. n) go to 145
!
        do i = l2, n
           d(i) = d(i) - h
        enddo
!
145     f = f + h
!     .......... ql transformation ..........
        p = d(m)
        c = 1.0d0
        c2 = c
        el1 = e(l1)
        s = 0.0d0
        mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
        do ii = 1, mml
           c3 = c2
           c2 = c
           s2 = s
           i = m - ii
           g = c * e(i)
           h = c * p
           r = pythag(p,e(i))
           e(i+1) = s * r
           s = e(i) / r
           c = p / r
           p = c * d(i) - s * g
           d(i+1) = h + s * (c * g + s * d(i))
!     .......... form vector ..........
           do k = 1, n
              h = z(k,i+1)
              z(k,i+1) = s * z(k,i) + c * h
              z(k,i) = c * z(k,i) - s * h
           enddo
!
        enddo
!
        p = -s * s2 * c3 * el1 * e(l) / dl1
        e(l) = s * p
        d(l) = c * p
        tst2 = tst1 + dabs(e(l))
        if (tst2 .gt. tst1) go to 130
220     d(l) = d(l) + f
     enddo
!     .......... order eigenvalues and eigenvectors ..........
     do ii = 2, n
        i = ii - 1
        k = i
        p = d(i)
!
        do j = ii, n
           if (d(j) .ge. p) cycle
           k = j
           p = d(j)
        enddo
!
        if (k .eq. i) cycle
        d(k) = d(i)
        d(i) = p
!
        do j = 1, n
           p = z(j,i)
           z(j,i) = z(j,k)
           z(j,k) = p
        enddo
!                                                                                                         
     enddo
!
     go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 return
end subroutine tql2
! **************************************************************************
!
subroutine round_int(rval,ival)
!
! This subroutine rounds ival (where int(rval) = ival) up or down
! accordingly, to avoid integer contraction. It works by checking
! the sign of rval and either adding or subtracting 1 sequentially
! until zero is passed. The rounding up or down of ival is decided
! according to the abolute value of rval after all the additions 
! or subtractions of 1.
!
! It is assumed that ival is the result of int(rval). The max
! absolute value of rval is 10^12
!
! J. Buckeridge November 2014
!
  implicit none
!
  integer, parameter :: large = 1e8
  integer i, j, k, ival
  real*8 rval
  logical signcheck
!
  if(rval > 0.d0) then
     signcheck = .true.
  elseif(rval < 0.d0) then
     signcheck = .false.
  else
     ival = 0
     return
  endif
!
  i = 0
  do
     i = i+1
     if(signcheck) then
        rval = rval - 1.d0
        if(rval < 0.d0) then
           if(abs(rval) > 0.5d0) then
              exit
           else
              ival = ival + 1
              exit
           endif
        endif
     else
        rval = rval + 1.d0
        if(rval > 0.d0) then
           if(abs(rval) > 0.5d0) then
              exit
           else
              ival = ival - 1
              exit
           endif
        endif
     endif
     if(i > large) then
        write(*,'(a)') "ERROR: integer checking routine failed!!"
        stop
     endif
  enddo
  return
!
end subroutine round_int
