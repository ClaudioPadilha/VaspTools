program MeanPotential 
  implicit none

  integer                                       :: nspecies = 1, natoms = 0, nx, ny, nz, ngrid, i, j, k, l, n, u = 97, v = 98
  integer, dimension (:), allocatable           :: natom
  real, dimension (:, :, :), allocatable        :: pot
  real, dimension (3, 3)                        :: cell
  real, dimension (3)                           :: lcell, dl
  real, dimension (:, :), allocatable           :: temp_pot_xy, temp_pot_xz, temp_pot_yz
  real, dimension (:), allocatable              :: potx, poty, potz
  real                                          :: escala, potm
  character (len=50)                            :: entrada, saida, progname
  character (len=90)                            :: temps
  character (len=2), dimension (:), allocatable :: atom_symbol

  call get_command_argument (0, progname)
  call get_command_argument (1, entrada)

  if (len_trim (entrada) .eq. 0) then
    print *, ' '
    print *, 'Erro: '
    print *, ' '
    print *, 'Sintaxe: ', len_trim (adjustl (progname)), ' [arquivo de input]'
    print *, ' '
    stop
  end if

! Open the data file
  open (u, FILE = entrada, STATUS = 'OLD')

! Ler a primeira linha do LOCPOT
  read (u, '(A)')  

! Ler fator de escala da rede
  read (u, *) escala

! Ler os vetors da rede
  read (u, *) cell

  cell = cell * escala

  lcell = 0

  do i = 1, 3
    do j = 1, 3
      lcell (i) = lcell (i) + cell (i, j) ** 2
    enddo
    lcell (i) = sqrt (lcell (i))
  enddo

! Este trecho le o numero de especies atomicas e seus simbolos
! e os aloca nos vetores natom e atom_symbol, respectivamente
  read (u, '(A)') temps 

  n = len_trim (temps)
  k = n - len_trim (adjustl (temps)) + 2

  do j = k, n
    if ((temps (j - 1:j - 1) /= ' ') .and. (temps (j:j) == ' ')) nspecies = nspecies + 1 
  enddo  

  allocate (atom_symbol (nspecies))
  allocate (natom (nspecies)) 

  i = 1
  j = 0
  k = 0
  l = 1

  if (temps (i:i) /= ' ') j = i

  do while (i <= len_trim (temps))
    if ((temps (i:i) == ' ') .and. (temps (i + 1: i + 1) /= ' ')) j = i + 1 
    if ((temps (i:i) /= ' ') .and. ((temps (i + 1: i + 1) == ' ') .or. (i + 1 == len_trim (temps)))) then
      atom_symbol (l) = temps (j:i)
      l = l + 1
    endif
    i = i + 1
  enddo 

  read (u, *) natom

!  do i = 1, nspecies
!    print *, 'Especie: ', atom_symbol(i), ', Quantidade: ', natom (i)
!  enddo
 
! le a proxima linha e joga fora
  read (u, *)

! le e joga fora as posicoes de todos os atomos
  do i = 1, nspecies
    natoms = natoms + natom (i)
  enddo

  do i = 1, natoms
    read (u, *)
  enddo
  
! le as divisoes do grid em x, y e z
  read (u, *) nx, ny, nz

! calcula o numero total de pontos no grid
  ngrid = nx * ny * nz
 
! aloca o array do potencial e le os valores do mesmo 
  allocate (pot (nx, ny, nz))

  allocate (potz (nz))
  allocate (poty (ny))
  allocate (potx (nx))

  allocate (temp_pot_xy (nx, ny))
  allocate (temp_pot_xz (nx, nz))
  allocate (temp_pot_yz (ny, nz))

  dl (1) = lcell (1) / nx
  dl (2) = lcell (2) / ny
  dl (3) = lcell (3) / nz
   
  read (u, *) pot

  close (u)

! integra em x e y
  temp_pot_xz = 0
  potz = 0

  do i = 1, nx
    do j = 1, nz
      do k = 1, ny
        temp_pot_xz (i, j) = temp_pot_xz (i, j) + pot (i, k, j)
      enddo
    enddo
  enddo

  do i = 1, nz
    do j = 1, nx
      potz (i) = potz (i) + temp_pot_xz (j, i)
    enddo
  enddo

  potz = potz / (nx * ny) 

  open (v, FILE = 'potz.dat', STATUS = 'NEW')
  
  do i = 1, nz
    write (v, *) i, i * dl (3),  potz (i)
  enddo

  close (v)

! integra em x e z
  temp_pot_xy = 0
  poty = 0

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        temp_pot_xy (i, j) = temp_pot_xy (i, j) + pot (i, j, k)
      enddo
    enddo
  enddo

  do i = 1, ny
    do j = 1, nx
      poty (i) = poty (i) + temp_pot_xy (j, i)
    enddo
  enddo

  poty = poty / (nx * nz) 

  open (v, FILE = 'poty.dat', STATUS = 'NEW')
  
  do i = 1, ny
    write (v, *) i, i * dl (2), poty (i)
  enddo

  close (v)

! integra em y e z
  temp_pot_xy = 0
  potx = 0

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        temp_pot_xy (i, j) = temp_pot_xy (i, j) + pot (i, j, k)
      enddo
    enddo
  enddo

  do i = 1, nx
    do j = 1, ny
      potx (i) = potx (i) + temp_pot_xy (i, j)
    enddo
  enddo

  potx = potx / (ny * nz) 

  open (v, FILE = 'potx.dat', STATUS = 'NEW')
  
  do i = 1, nx
    write (v, *) i, i * dl (1), potx (i)
  enddo

  close (v)

! integra em x, y e z
  potm = 0

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        potm = potm + pot (i, j, k)
      enddo
    enddo
  enddo

  potm = potm / ngrid 

  open (v, FILE = 'potm.dat', STATUS = 'NEW')
  
  write (v, *) potm

  close (v)

end program MeanPotential 

