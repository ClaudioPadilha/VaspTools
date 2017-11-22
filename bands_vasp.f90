program VaspBands
  implicit none

  integer                            :: nband, nkpt, nkpath, npath, l, j, i, u = 99, v = 98, t = 97, s = 96
  real                               :: k = 0.d0, dk, ksimb = 0.d0, deltak, lixo, ef, dk0
  real, dimension (3)                :: ki = (/0.d0, 0.d0, 0.d0/), kf = (/0.d0, 0.d0, 0.d0/)
  real, dimension (3)                :: kpathi = (/0.d0, 0.d0, 0.d0/), kpathf = (/0.d0, 0.d0, 0.d0/), kpathtemp
  real, dimension (:,:), allocatable :: eigup, eigdown
  real, dimension (:), allocatable   :: simbolos
  character (20)                     :: strlixo, aux, simb, ksimbi = "IIIIII", ksimbf = "FFFFFF", ksimbtemp = "TTTTTT"
  character (50)                     :: entrada_eig, entrada_kpt, saida_dados, saida_script
  character                          :: sp

  call get_command_argument (1, entrada_eig)
  call get_command_argument (2, entrada_kpt)

  if ((len_trim (entrada_eig) .eq. 0) .or. (len_trim (entrada_kpt) .eq. 0)) then
    print *, ' '
    print *, 'Error: '
    print *, ' '
    print *, 'Sintax: bandas_vasp.x [input file EIGENVAL] [input file KPOINTS]'
    print *, ' '
    stop
  end if

!  print *, ' '
!  print *, 'Digite o nome do arquivo dos autovalores: (formato do EIGENVAL)'
!  read *, entrada_eig

!  print *, ' '
!  print *, 'Digite o nome do arquivo com o caminho dos pontos k: (formato do KPOINTS)'
!  read *, entrada_kpt
  
  ! Open the data file
  open (unit = u, FILE = entrada_eig, STATUS = 'OLD')
  open (unit = s, FILE = entrada_kpt, STATUS = 'OLD')

  ! Arquivo de output dos dados desse programa
  saida_dados = 'bands-'//trim(entrada_eig)//'.dat'
  open (unit = v, file = saida_dados, status ='replace') 

  ! Arquivo de script do gnuplot
  saida_script = 'plot-'//trim(entrada_eig)//'.gnu'
  open (unit = t, file = saida_script, status = 'replace') 

  ! primeira linha, ultimo registro do EIGENVAL tem o valor de NSPIN
  read (u, *) lixo, lixo, lixo, sp

  ! ler as 4 linhas seguintes do EIGENVAL e jogar fora
  do i = 1, 4
    read (u, *)
  enddo

  ! numero de bandas e numero de kpoints do EIGENVAL
  read (u, *) lixo, nkpt, nband  

  ! ler a primeira linha do KPOINTS e jogar fora 
  read (s, *)

  ! ler o numero de divisoes nos caminhos do KPOINTS
  read (s, *) nkpath

  ! calcula o numero de trechos no KPOINTS
  npath = nkpt / nkpath

  ! aloca os vetores que guardam os autovalores
  allocate (eigup (nkpath, nband + 1))

  if (sp .eq. '2') then
    allocate (eigdown (nkpath, nband + 1))
  endif

  ! aloca vetor com os valores dos k's dos pontos de alta simetria
  allocate (simbolos (npath - 1))

  ! ler as duas proximas linhas do KPOINTS e jogar fora
  read (s, *)
  read (s, *)

  print *, ' '
  print *, 'SP = ', sp
  print *, ' '
  write (*,*) 'Ef = ?'
  read (*,*) ef

  ! gnuplot script starts here
  write (t, *) "set terminal postscript eps color enhanced font 22"
  write (t, *) "set ylabel 'E-E_F (eV)'" 
  write (t, *) "set pointsize 0.5"
  write (t, *) "set border lw 4"
  write (t, *) "unset key"
  write (t, *) "Ef = ", ef
  write (t, *) "set output 'bandas-", trim(entrada_eig), ".eps'"
  write (t, *) '# Edite as linhas abaixo ao seu bel-prazer!'
  write (t, *) 'ymin = -10.0'
  write (t, *) 'ymax =  10.0'
  write (t, *) 'set yrange [ymin:ymax]'
  write (t, *) 'set ytics ymin ymax 1.0'
  write (t, *) 'set style line 1 lt 2 lw 1 lc 0'
  write (t, *) 'set style arrow 1 nohead ls 1'
  !write (t, *) 'set xtics ( \' !, trim(kpath), ')'

  ! Esse loop le as informacoes do KPOINTS e gera a string com os
  ! simbolos dos pontos de alta simetria para ser utilizada como
  ! label no eixo x do grafico, assim como as posicoes em x deles
  ! para que sejam geradas as linhas horizontais que deixam o aspecto
  ! da coisa um tanto quanto mais profissional! =p
  ! Detalhe importante! As linhas do KPOINTS devem estar no formato:
  !       0.000   0.000   0.000   !\Gamma 
  !       0.000   0.000   0.500   !Y 
  ! e nunca no formato:
  !   0.000   0.000   0.000   ! \Gamma 
  !   0.000   0.000   0.500   ! Y
  ! pois este espaco entre o ! e o simbolo ferra tudo! 
  
  do i = 1, npath
    read (s, *) kpathi, ksimbi
    read (s, *) kpathf, ksimbf 

!    print *, ksimbi, ksimbf

    deltak = sqrt ((kpathf(1)-kpathi(1))**2+(kpathf(2)-kpathi(2))**2+(kpathf(3)-kpathi(3))**2)

    dk0 = deltak / nkpath

    write (aux, '(F10.7)') ksimb
 
    if (i .eq. 1) then
      write (t, *) 'set xtics ("'//trim(ksimbi(2:len_trim(ksimbi)))//'" '//trim(aux)//', \'
    elseif (trim(ksimbtemp) .eq. trim(ksimbi)) then 
      write (t, *) '           "'//trim(ksimbtemp(2:len_trim(ksimbtemp)))//'" '//trim(aux)//', \'
      simbolos(i) = ksimb
      !print *, simbolos(i)
    else
      write (t, *) '           "'//trim(ksimbtemp(2:len_trim(ksimbtemp)))//'|' &
        & //trim(ksimbi(2:len_trim(ksimbi)))//'" '//trim(aux)//', \'
      simbolos(i) = ksimb
      !print *, simbolos(i)
    endif

    if (i .eq. npath) then
      write (aux, '(F10.7)') ksimb + deltak
      write (t, *) '           "'//trim(ksimbf(2:len_trim(ksimbf)))//'" '//trim(aux) //')'
      simbolos(i + 1) = ksimb + deltak
    endif

    ksimbtemp = ksimbf

    ksimb = ksimb + deltak 

  enddo

  write (aux, '(F4.1)') ksimb + 0.5 * dk0
  write (t, *) 'set xrange [0.0:' // trim (aux) // ']' 

  close(s)

  do i = 2, npath
    write (t, *) 'set arrow from ', simbolos (i), ', ymin to ', simbolos (i), ', ymax as 1'
  enddo

  if (sp .eq. '2') then
    write (t, *) 'plot "', trim(saida_dados), '" u 1:($2-Ef) pt 1 lt 0 lc 2, \'
    write (t, *) '     "', trim(saida_dados), '" u 1:($3-Ef) pt 2 lt 0 lc 3, \'
    write (t, *) '      0 w l lt 2 lw 4'
  else
    write (t, *) 'plot "', trim(saida_dados), '" u 1:($2-Ef) pt 1 lt 0 lc 2, \'
    write (t, *) '      0 w l lt 2 lw 4'
  endif
  
  close(t)

  ! Este loop le os autovalores por ponto k do arquivo EIGENVAL e os armazena em um
  ! vetor eigup (caso nao haja polarizacao de spin, sp = 1) ou nos vetores eigup e 
  ! eigdown (com polarizacao de spin, sp = 2) e coloca esses dados em um arquivo de 2 ou 
  ! colunas de dados com os valores de |k| no eixo x. Esse eh o arquivo de dados do 
  ! script do gnuplot.
  do i = 1, nkpt
    read (u, *) (kf(j), j = 1, 3), lixo

    dk = sqrt ((kf(1) - ki(1))**2 + (kf(2) - ki(2))**2 + (kf(3) - ki(3))**2)

    if (dk .gt. 10 * dk0) then
      print *, 'Pontos k descontinuos! dk = ', dk
      print *, 'dk esperado da ordem de: ', dk0
      dk = 0.d0
    endif

    ! guarda o valor do k-point na primeira posição dos vetores de autovalores
    eigup (i, 1) = k
    if (sp .eq. '2') then
      eigdown (i, 1) = k
    endif

    ! print *, 'k-point ', i, ': ', (kf(l), l = 1, 3), ' ... ok, dk = ', dk
    do j = 2, nband + 1
      if (sp .eq. '2') then
        read (u, *) lixo, eigup (i,j), eigdown (i,j)
        write (v, *) k, '  ', eigup (i,j), '  ', eigdown (i,j)
      else
        read (u, *) lixo, eigup (i,j)
        write (v, *) k, '  ', eigup (i,j)
      endif
      !print *, j, eigup (i, j), eigdown (i,j)
    enddo
    
    k = k + dk

    ki = kf

  enddo

  close(u)
  close(v)

  !print *, ' '
  !print *, trim(kpath)
  !print *, (simbolos(i), i = 2, npath)

  print *, ' '
  print *, 'Data wrote in ', trim(saida_dados)
  !print *, ' '
  print *, 'Use the script ', trim(saida_script), ' to plot the graph'

end program VaspBands
