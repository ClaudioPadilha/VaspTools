module mod_dos
  implicit none
  
  public  :: col_count, row_count, atomic_species, smooth

  contains
    subroutine atomic_species (poscar, ntype, satoms, natoms)
      character (len=50), intent (in)                               :: poscar
      integer, intent (inout)                                       :: ntype
      character (len=2), dimension (:), allocatable, intent (inout) :: satoms
      integer, dimension (:), allocatable, intent (inout)           :: natoms
      
      integer  :: i = 1, j, u = 97
      character (len=50) :: temp

      open (u, FILE = poscar, STATUS = 'OLD')

      ! reads and ignores first 5 lines of POSCAR
      do i = 1, 5
        read (u, '(A)')
      end do
      
      ! reads the line that has the atomic symbols
      read (u, '(A)') temp

      ! counts the number of collumns (and atomic species)
      ntype = col_count (temp)

      ! put the symbols found in a vector
      allocate (satoms (ntype)) 
      read (temp, *) (satoms (i), i = 1, ntype)

      ! reads the line that has the number of each atomic species
      read (u, '(A)') temp

      ! puts the numbers in a vector
      allocate (natoms (ntype))
      read (temp, *) (natoms (i), i = 1, ntype)

      close (u)
      write (*, '(A,I3)') " Number of atomic species: ", ntype
      do i = 1, ntype
        write (*, '(A,I3,A,I3)') " species [",i, "] : " // trim(satoms (i)) // ', number = ', natoms (i)
      end do

    end subroutine atomic_species

    ! returns the number of registers in a single string
    function col_count (line)
      character (len=*), intent (in)  :: line
      integer                         :: n, i, j, ini, fim, col_count
      
      col_count = 1

      do i = 1, len (line)
        if (line (i:i) /= ' ') then 
          ini = i
          exit
        end if
        !if (c .and. line (i:i) /= ' ' .and. line (i + 1 : i + 1) == ' ') then 
        !  j = j + 1
        !  print *, "space found: [", line(i:i+4), "]"
        !end if
      end do
      fim = ini + len_trim(adjustl(line))
      !print *, ini, fim, "[", line(ini:fim), "]"

      do i = ini, fim - 1
        if (line(i:i) == ' ' .and. line (i+1:i+1) /= ' ') col_count = col_count + 1
      end do
      !col_count = j

    end function col_count

    function row_count (u, n)
      integer   :: n, u, row_count, i, stat = 0, j = 0
      real      :: temp

      do
        read (u, *, IOSTAT = stat) (temp, i = 1, n)
        if (stat /= 0) exit
        j = j + 1
      end do

      rewind (u)

      row_count = j

    end function row_count

    subroutine smooth (sigma, dos, ndiv, dos_s)
      integer, intent (in)                                :: ndiv
      real, intent (in)                                   :: sigma
      real, dimension (:, :), intent (in)                 :: dos
      real, dimension (:, :), allocatable, intent (out)   :: dos_s

      integer                                             :: i, j, k, n, m, p, q
      double precision                                    :: dE, Ei, Ef, dif = 1.0d0 
      integer, parameter                                  :: pi = 3.14159265

      n = size (dos, 1)
      m = size (dos, 2)
      p = ndiv * (m - 1)

      allocate (dos_s (n, p)) 
      
      dos_s = 0.0d0
      Ei = dos (1, 1)
      Ef = dos (1, m)

      dE = (Ef - Ei) / p
     
      do i = 1, p + 1
        dos_s (1, i) = Ei + (i - 1) * dE 
        if (abs (dos_s (1, i) - dos_s (1, 1) - 3 * sigma) .lt. dif) then
           q = i
           dif = abs (dos_s (1, i) - dos_s (1, 1) - 3 * sigma)
        end if
      end do

      do i = 2, n
        do j = 1, m
          if ((j - 1) * ndiv - q .lt. 1) then
            do k = 1, (j - 1) * ndiv + q
              dos_s (i, k) = dos_s (i, k) + dos (i, j) * exp (- 0.5 * ((dos_s (1, k) - dos (1, j)) / sigma) ** 2)
            end do
          else if ((j - 1) * ndiv + q .gt. p) then
            do k = (j - 1) * ndiv - q, p 
              dos_s (i, k) = dos_s (i, k) + dos (i, j) * exp (- 0.5 * ((dos_s (1, k) - dos (1, j)) / sigma) ** 2)
            end do
          else
            do k = (j - 1) * ndiv - q, (j - 1) * ndiv + q 
              dos_s (i, k) = dos_s (i, k) + dos (i, j) * exp (- 0.5 * ((dos_s (1, k) - dos (1, j)) / sigma) ** 2)
            end do
          end if 
        end do   
      end do
      
    end subroutine smooth

end module mod_dos

program ProjectedDensityOfStates
  use mod_dos
  implicit none

  integer                                        :: natom, n, nenergy, i, j, k, l, m, u = 97, v = 98, w = 99, lixoi, tempi = 0
  integer                                        :: cols, ispin, col_dost, col_dosp, stat = 0, lorbit
  real, dimension (:, :, :), allocatable         :: dosp
  real, dimension (:, :), allocatable            :: dost
  integer, dimension (:), allocatable            :: na
  real                                           :: fermi, lixor, tempr, temps
  character (len=50)                             :: doscar, poscar, saidat, saidag
  character (len=400)                            :: temp4
  character (len=50), dimension (:), allocatable :: saidap
  character (len=2), dimension (:), allocatable  :: tipo
  character (len=20)                             :: lixoc
  character (len=2)                              :: temp2
  character (len=3)                              :: temp, temp3 
  character (len=1)                              :: op
  logical                                        :: orboutput, atomoutput, forbital
  character (len=1), dimension (3)               :: orb0 = (/'s', 'p', 'd'/)
  character (len=1), dimension (4)               :: orb1 = (/'s', 'p', 'd', 'f'/)
  character (len=6), dimension (9)               :: orb2 = (/'s     ', 'py    ', 'pz    ', 'px    ', 'dxy   ', 'dyz   ', &
	&  'dz2   ', 'dxz   ', 'dx2-y2'/) 
  character (len=9), dimension (16)              :: orb3 = (/'s              ', 'py             ', 'pz             ', & 
	& 'px             ',  'dxy            ', 'dyz            ', 'dz2            ', 'dxz            ', 'dx2-y2         ', &
	& 'fy(3x2-y2)     ', 'fxyz           ', 'fy(4z2-x2-y2)  ',  'fz(2z2-3x2-3y2)', 'fx(4z2-x2-y2)  ', 'fz(x2-y2)      ', &
	& 'fx(x2-3y2)     '/)

  call get_command_argument (1, doscar)
  call get_command_argument (2, poscar)
  call get_command_argument (3, op)

  if (len_trim (doscar) .eq. 0) then
    print *, ' '
    print *, 'Error: '
    print *, ' '
    print *, 'Sintax: dos_vasp.x [input file DOSCAR] [POSCAR file] (output option)'
    print *, ' '
    print *, ' [ ... ] mandatory'
    print *, ' ( ... ) optional, A = pDOS per atom, O = pDOS per orbital, X = both'
    print *, ' '
    stop
  end if

  if (op == 'O' .or. op == 'o') then
    orboutput = .true.
    atomoutput = .false.
    print *, "pDOS per orbital selected (s, px, py, pz, dxy, dyz, dz^2, dxz, and dx^2)"
  else if (op == 'A' .or. op == 'a') then
    atomoutput = .true.
    orboutput = .false.
    print *, "pDOS per atom selected (eg. DOS-Ti1, DOS-Ti2, etc.)"
  else if (op == 'X' .or. op == 'x') then
    atomoutput = .true.
    orboutput = .true.
    print *, "pDOS per orbital AND atom selected"
    print *, "  (s, px, py, pz, dxy, dyz, dz^2, dxz, and dx^2)"
    print *, "  (eg. DOS-Ti1, DOS-Ti2, etc.)"
  else
    orboutput = .false.
    atomoutput = .false.
    print *, "simple output selected (but still you get s, p and d projections)"
  end if

  ! read atomic species info from POSCAR
  call atomic_species (poscar, n, tipo, na)
  
! Open the DOSCAR file
  open (u, FILE = doscar, STATUS = 'OLD')

  ! Ler o numero de atomos e depois lixos
  read (u, *) natom, lixoi, lixoi, lixoi 

  do i = 1, n
    tempi = tempi + na(i)
  enddo

  if (natom /= tempi) then
    print *, 'Error: the ammount of atoms in ', trim(doscar), ' and ', trim(poscar), ' are not the same!'
    stop 
  end if

! ler lixos  
  read (u, *)
  read (u, *)
  read (u, *) 
  read (u, *) 
  
! Ler numero de bins de energia, energia de Fermi e lixos  
  read (u, *) lixor, lixor, nenergy, fermi, lixor

  print *, '# of atoms:', natom, ' # of bins:', nenergy, ' Fermi Energy:', fermi

  !here we need to know if ISPIN=1 or 2
  ! and allocate the total dos (dost) array accordingly
  read (u, '(A)') temp4
  cols = col_count (temp4)
  backspace (u)
  if (cols == 3) then
    col_dost = 2
    print *, "Calculation without spin polarization (ISPIN=1)"
    ispin = 1
  else if (cols == 5) then
    col_dost = 3
    print *, "Calculation with spin polarization (ISPIN=2)"
    ispin = 2
  else
    print *, "Error: Check ISPIN flag. Not recognized!"
    stop
  end if

  ! alocate the vector with total dos
  allocate (dost (col_dost, nenergy))

  ! reads the DOSCAR and stores the info for total dos
  if (ispin == 1) then
    do i = 1, nenergy
      read (u, *) (dost (k, i), k = 1, col_dost), lixor
    end do
  else
    do i = 1, nenergy
      read (u, *) (dost (k, i), k = 1, col_dost), lixor, lixor
    end do
  end if

  saidag = 'plot-dos-' // trim(doscar) // '.gnu'

  ! gnuplot script is written
  open (w, FILE = saidag, STATUS = 'replace')

  ! header of script .gnu
  write (w, *) "set terminal postscript eps color enhanced font 22"
  write (w, *) "set xlabel 'E-E_F (eV)'"
  write (w, *) "set ylabel 'DOS'" 
  write (w, *) "set border lw 4"
  write (w, *) "Ef = ", fermi 
  write (w, *) "set output '", trim(doscar), ".eps'"
  write (w, *) "set size ratio 0.3"
  write (w, *) '# Feel free to edit the following lines!'
  write (w, *) '#xmin = -10.0'
  write (w, *) '#xmax =  10.0'
  write (w, *) '#set xrange [xmin:xmax]'
  write (w, *) '#set xtics xmin xmax 2.0 format "%.1f"'
  write (w, *) '#unset ytics'
  write (w, *) '#set key left top'
  write (w, *) 'set style line 1 lt 2 lw 2 lc rgb "red"'
  write (w, *) 'set style arrow 1 nohead ls 1'
  write (w, *) '#ymin = -70'
  write (w, *) '#ymax =  70'
  write (w, *) '#set yrange [ymin:ymax]'
  write (w, *) '#set arrow from 0.0,ymin to 0.0,ymax as 1'

  ! escrever a dos total
  saidat = trim(doscar) // '-dos-total.dat'

  open (v, FILE = saidat, STATUS = 'replace')

  if (ispin == 2) then
    do i = 1, nenergy
      write (v, *) dost (1, i), ' ', dost (2, i), ' ', -1 * dost (3, i)
    end do
  else if (ispin == 1) then
    do i = 1, nenergy
      write (v, *) dost (1, i), ' ', dost (2, i)
    end do
  end if
  
  write (w, *) "plot '" // trim(saidat) // "' u ($1-Ef):2 w l lw 3 lt 1 lc 0 title 'Total', \"
  if (ispin == 2) write (w, *) "     '" // trim(saidat) // "' u ($1-Ef):3 w l lw 3 lt 1 lc 0 notitle, \"

  close (v)

  ! here we decide if LORBIT = 10, 11, or none
  read (u, *, IOSTAT = stat)

  if (stat /= 0) then
    print *, "No LORBIT option!"
    stop
  else
    read (u, '(A)') temp4
    cols = col_count (temp4)
    backspace (u) 
    backspace (u)
    !print *, cols, "[", trim(temp4), "]"
  end if

  if (cols == 4) then
    ! ISPIN = 1, LORBIT = 10, no f orbitals in pseudo
    ! E s p d
    lorbit = 10
    forbital = .false.
  else if (cols == 5) then
    ! ISPIN = 1, LORBIT = 10, and f orbitals in pseud
    ! E s p d f
    lorbit = 10
    forbital = .true.
  else if (cols == 10) then
    ! ISPIN = 1, LORBIT = 11, no f orbitals in pseudo
    ! E s py pz px dxy dyz dz2 dxz dx2
    lorbit = 11
    forbital = .false.
  else if (cols == 17) then
    ! ISPIN = 1, LORBIT = 11, f orbitals present in pseudo
    ! E s py pz px dxy dyz dz2 dxz dx2 fy(3x2-y2) fxy fy(4z2-x2-y2) fz(2z2-3x2-3y2) fx(4z2-x2-y2) fz(x2-y2) fx(x2-3y2)
    lorbit = 11
    forbital = .true.
  else if (cols == 7) then
    ! ISPIN = 2, LORBIT = 10, no f orbitals
    ! E s(up) s(down) p(up) p(down) d(up) d(down)
    lorbit = 10
    forbital = .false.
  else if (cols == 9) then
    ! ISPIN = 2, LORBIT = 10, f orbitals present
    ! E s(up) s(down) p(up) p(down) d(up) d(down) f(up) f(down)
    lorbit = 10
    forbital = .true.
  else if (cols == 19) then
    ! ISPIN = 2, LORBIT = 11, no f orbitals
    ! E s(u) s(d) py(u) py(d) pz(u) pz(d) px(u) px(d) dxy(u) dxy(d) dyz(u) dyz(d) dz2(u) dz2(d) dxz(u) dxz(d) dx2(u) dx2(d)
    lorbit = 11
    forbital = .false.
  else if (cols == 33) then
    ! ISPIN = 2, LORBIT = 11, f orbitals
    ! E s(u) s(d) py(u) py(d) pz(u) pz(d) px(u) px(d) dxy(u) dxy(d) dyz(u) dyz(d) dz2(u) dz2(d) dxz(u) dxz(d) dx2(u) dx2(d)
    !     fy(3x2-y2)(u) fy(3x2-y2)(d) fxy(u) fxyz(d) fy(4z2-x2-y2)(u) fy(4z2-x2-y2)(d) fz(2z2-3x2-3y2)(u) fz(2z2-3x2-3y2)(d)
    !     fx(4z2-x2-y2)(u) fx(4z2-x2-y2)(d) fz(x2-y2)(u) fz(x2-y2)(d) fx(x2-3y2)(u) fx(x2-3y2)(d)
    lorbit = 11
    forbital = .true.
  else
    print *, 'Error: LORBIT flag not recognized!'
    stop
  end if
  write (*, '(A,I2)') " LORBIT=", lorbit
  if (forbital) then
    print *, "f orbitals present"
  else
    print *, "f orbitals not present"
  end if
  ! allocate pdos considering ISPIN and LORBIT  
  allocate (dosp (natom, cols, nenergy))
  allocate (saidap (natom))

  ! le a dos projetada e escreve a dos projetada por atomo
  tempi = 0
  do l = 1, n
    do j = 1, na(l)
      read (u, *)

      write (temp, '(I3)') j

      if (atomoutput .eqv. .true.) then
        saidap(j) = trim(doscar) // '-dos-proj-' // trim(tipo(l)) // '-' // trim(adjustl(temp)) &
            & // '.dat'
        open (v, FILE = saidap(j), STATUS = 'replace') 
      end if

      do i = 1, nenergy
        read (u, *) (dosp (j + tempi, k, i), k = 1, cols)
        if (atomoutput .eqv. .true.) then
          if (ispin == 1) then
            write (v, *) dosp (j + tempi, 1, i), (' ', dosp (j + tempi, k, i), k = 2, cols)
          else
            write (v, *) dosp (j + tempi, 1, i), (' ', dosp (j + tempi, k, i), ' ', -1 * dosp (j + tempi, k + 1, i), k = 2, & 
		& cols - 1, 2)
          end if
        end if
      end do

!     parte do script do gnuplot. Esses arquivos devem ser utilizados escolhendo-se a coluna a ser plotada, na ordem do
!     array orb.
      if (atomoutput .eqv. .true.) then
        write (w, *) "     '" // trim(saidap(j)) // "' u ($1-Ef):2 w l lw 1 lc rgb 'color'" // &
                   & " title '" // trim(tipo(l)) // "_{" // trim(temp) // "}', \"
        if (ispin == 2) then
          write (w, *) "     '" // trim(saidap(j)) // "' u ($1-Ef):3 w l lw 1 lc rgb 'color'" // &
                   & " notitle, \" 
        end if
        close (v)
      end if

    end do
    tempi = tempi + na(l)
  end do

! dos projetada por orbital de cada especie atomica
  if (orboutput .eqv. .false.) then
    m = 1 ! do it only for s
  else 
    if (lorbit == 10) then
      if (forbital) then
        m = 4 ! do it also for p, d, and f
      else
        m = 3 
      end if
    else 
      if (forbital) then
        m = 16
      else
        m = 9 ! do it also for px, py, pz, dxy, dxz, etc...
      end if
    end if
  end if
  
  do k = 1, m
    tempi = 0
    write (temp, '(I2)') k + 1
    do l = 1, n
      if (lorbit == 10) then
        if (forbital) then
          saidat = trim(doscar) // '-dos-proj-' // trim (tipo(l)) // '-' // trim (orb1(k)) // '.dat'
        else
          saidat = trim(doscar) // '-dos-proj-' // trim (tipo(l)) // '-' // trim (orb0(k)) // '.dat'
        end if
      else 
        if (forbital) then
          saidat = trim(doscar) // '-dos-proj-' // trim (tipo(l)) // '-' // trim (orb3(k)) // '.dat'
        else
          saidat = trim(doscar) // '-dos-proj-' // trim (tipo(l)) // '-' // trim (orb2(k)) // '.dat'
        end if
      end if
      open (v, FILE = saidat, STATUS = 'replace')
      
      do i = 1, nenergy
        tempr = 0.d0
        temps = 0.d0

        if (ispin == 2) then
          do j = tempi + 1, na(l) + tempi
            tempr = tempr + dosp (j, 2 * k, i)
            temps = temps + dosp (j, 2 * k + 1, i)
          end do

          write (v, *) dosp (1, 1, i), ' ', tempr, ' ', -1 * temps
        else
          do j = tempi + 1, na(l) + tempi
            tempr = tempr + dosp (j, k + 1, i)
          end do

          write (v, *) dosp (1, 1, i), ' ', tempr
        end if
      end do

      write (temp2, '(I2)') l

      if (lorbit == 10) then
        write (w, *) "     '" // trim(saidat) // "' u ($1-Ef):2 w l lw 2 lt 2 lc rgb 'color'" // & 
              &  " title '" // trim(tipo(l)) // "_{" // trim(orb1(k)) // "}', \"
      else 
        write (w, *) "     '" // trim(saidat) // "' u ($1-Ef):2 w l lw 2 lt 2 lc rgb 'color'" // & 
              &  " title '" // trim(tipo(l)) // "_{" // trim(orb2(k)) // "}', \"
      end if
      if (ispin == 2) then
        write (w, *) "     '" // trim(saidat) // "' u ($1-Ef):3 w l lw 2 lt 2 lc rgb 'color'" // &
              &  " notitle, \"
      end if

      tempi = tempi + na(l)
      close (v)
    enddo
  enddo

! orbital p total
  tempi = 0
  do l = 1, n
    saidat = trim(doscar) // '-dos-proj-' // trim (tipo(l)) // '-' // trim (orb1(2)) // '.dat'

    open (v, FILE = saidat, STATUS = 'replace')

    if (ispin == 2) then
      if (lorbit == 11) then
        do i = 1, nenergy
          tempr = 0.d0
          temps = 0.d0

          do j = tempi + 1, na(l) + tempi
            tempr = tempr + dosp (j, 4, i) + dosp (j, 6, i) + dosp (j, 8, i)
            temps = temps + dosp (j, 5, i) + dosp (j, 7, i) + dosp (j, 9, i)
          end do

          write (v, *) dosp (1, 1, i), ' ', tempr, ' ', -1 * temps
        end do
      else
        do i = 1, nenergy
          tempr = 0.d0
          temps = 0.d0

          do j = tempi + 1, na(l) + tempi
            tempr = tempr + dosp (j, 4, i)
            temps = temps + dosp (j, 5, i)
          end do

          write (v, *) dosp (1, 1, i), ' ', tempr, ' ', -1 * temps
        end do
      end if
    else
      if (lorbit == 11) then
        do i = 1, nenergy
          tempr = 0.d0

          do j = tempi + 1, na(l) + tempi
            tempr = tempr + dosp (j, 3, i) + dosp (j, 4, i) + dosp (j, 5, i)
          end do

          write (v, *) dosp (1, 1, i), ' ', tempr
        end do
      else
        do i = 1, nenergy
          tempr = 0.d0

          do j = tempi + 1, na(l) + tempi
            tempr = tempr + dosp (j, 3, i)
          end do

          write (v, *) dosp (1, 1, i), ' ', tempr
        end do
      end if
    end if

    write (temp2, '(I2)') l

    write (w, *) "     '" // trim(saidat) // "' u ($1-Ef):2 w l lw 2 lt 2 lc " &
      // "rgb 'color' title '" // trim(tipo(l)) // "_p', \"
    if (ispin == 2) then
      write (w, *) "     '" // trim(saidat) // "' u ($1-Ef):3 w l lw 2 lt 2 lc " &
        // "rgb 'color' notitle, \"
    end if

    tempi = tempi + na(l)
    close (v)
  end do


! orbital d total
  tempi = 0
  do l = 1, n
    saidat = trim(doscar) // '-dos-proj-' // trim (tipo(l)) // '-d.dat'

    open (v, FILE = saidat, STATUS = 'replace')

    if (ispin == 2) then
      if (lorbit == 11) then
        do i = 1, nenergy
          tempr = 0.d0
          temps = 0.d0

          do j = tempi + 1, na(l) + tempi
            tempr = tempr + dosp (j, 10, i) + dosp (j, 12, i) + dosp (j, 14, i) + dosp (j, 16, i) + dosp (j, 18, i)
            temps = temps + dosp (j, 11, i) + dosp (j, 13, i) + dosp (j, 15, i) + dosp (j, 17, i) + dosp (j, 19, i) 
          end do

          write (v, *) dosp (1, 1, i), ' ', tempr, ' ', -1 * temps
        end do
      else
        do i = 1, nenergy
          tempr = 0.d0
          temps = 0.d0

          do j = tempi + 1, na(l) + tempi
            tempr = tempr + dosp (j, 6, i)
            temps = temps + dosp (j, 7, i)
          end do

          write (v, *) dosp (1, 1, i), ' ', tempr, ' ', -1 * temps
        end do
      end if
    else
      if (lorbit == 11) then
        do i = 1, nenergy
          tempr = 0.d0

          do j = tempi + 1, na(l) + tempi
            tempr = tempr + dosp (j, 6, i) + dosp (j, 7, i) + dosp (j, 8, i) + dosp (j, 9, i) + dosp (j, 10, i)
          end do

          write (v, *) dosp (1, 1, i), ' ', tempr
        end do
      else
        do i = 1, nenergy
          tempr = 0.d0

          do j = tempi + 1, na(l) + tempi
            tempr = tempr + dosp (j, 4, i)
          end do

          write (v, *) dosp (1, 1, i), ' ', tempr
        end do
      end if
    end if

    write (temp2, '(I2)') l

    write (w, *) "     '" // trim(saidat) // "' u ($1-Ef):2 w l lw 2 lt 2 lc " &
      // "rgb 'color' title '"  // trim(tipo(l)) // "_d', \"
    if (ispin == 2) then
      write (w, *) "     '" // trim(saidat) // "' u ($1-Ef):3 w l lw 2 lt 2 lc " &
      // "rgb 'color' notitle, \"
    end if

    tempi = tempi + na(l)
    close (v)
  enddo
  
  if (forbital) then
    ! orbital f total
    tempi = 0
    do l = 1, n
      saidat = trim(doscar) // '-dos-proj-' // trim (tipo(l)) // '-f.dat'

      open (v, FILE = saidat, STATUS = 'replace')

      if (ispin == 2) then
        if (lorbit == 11) then
          do i = 1, nenergy
            tempr = 0.d0
            temps = 0.d0

            do j = tempi + 1, na(l) + tempi
              tempr = tempr + dosp (j, 20, i) + dosp (j, 22, i) + dosp (j, 24, i) + dosp (j, 26, i) + &
              & dosp (j, 28, i) + dosp (j, 30, i) + dosp (j, 32, i)
              temps = temps + dosp (j, 21, i) + dosp (j, 23, i) + dosp (j, 25, i) + dosp (j, 27, i) + & 
              & dosp (j, 29, i) + dosp (j, 31, i) + dosp (j, 33, i)  
            end do

            write (v, *) dosp (1, 1, i), ' ', tempr, ' ', -1 * temps
          end do
        else
          do i = 1, nenergy
            tempr = 0.d0
            temps = 0.d0

            do j = tempi + 1, na(l) + tempi
              tempr = tempr + dosp (j, 8, i)
              temps = temps + dosp (j, 9, i)
            end do

            write (v, *) dosp (1, 1, i), ' ', tempr, ' ', -1 * temps
          end do
        end if
      else
        if (lorbit == 11) then
          do i = 1, nenergy
            tempr = 0.d0

            do j = tempi + 1, na(l) + tempi
              tempr = tempr + dosp (j, 11, i) + dosp (j, 12, i) + dosp (j, 13, i) + dosp (j, 14, i) + &
              & dosp (j, 15, i) + dosp (j, 16, i) + dosp (j, 17, i)
            end do

            write (v, *) dosp (1, 1, i), ' ', tempr
          end do
        else
          do i = 1, nenergy
            tempr = 0.d0

            do j = tempi + 1, na(l) + tempi
              tempr = tempr + dosp (j, 5, i)
            end do

            write (v, *) dosp (1, 1, i), ' ', tempr
          end do
        end if
      end if

      write (temp2, '(I2)') l

      write (w, *) "     '" // trim(saidat) // "' u ($1-Ef):2 w l lw 2 lt 2 lc " &
        // "rgb 'color' title '"  // trim(tipo(l)) // "_f', \"
      if (ispin == 2) then
        write (w, *) "     '" // trim(saidat) // "' u ($1-Ef):3 w l lw 2 lt 2 lc " &
        // "rgb 'color' notitle, \"
      end if

      tempi = tempi + na(l)
      close (v)
    enddo
  end if

  close (w)

end program ProjectedDensityOfStates
