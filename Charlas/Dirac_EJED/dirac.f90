! ==============================================================================
! RESOLUTOR DE LA ECUACIÓN DE DIRAC 1D - MALLA ESCALONADA (STAGGERED GRID)
! ==============================================================================
module params
    implicit none
    ! Precisión y constantes numéricas
    integer, parameter :: dp = selected_real_kind(15, 307)
    complex(dp), parameter :: img = cmplx(0.0_dp, 1.0_dp, kind=dp)
    
    ! Parámetros Físicos
    real(dp), parameter :: mass = 1.0_dp
    real(dp), parameter :: c_light = 1.0_dp
    
    ! Parámetros Topológicos del Dominio
    real(dp), parameter :: L = 50.0_dp
    integer, parameter  :: N = 2000
    real(dp), parameter :: dx = (2.0_dp * L) / real(N - 1, dp)
    
    ! Parámetros Temporales
    real(dp), parameter :: dt = 0.005_dp
    integer, parameter  :: Nt = 2000
    real(dp), parameter :: tau = dt / 2.0_dp
    
    ! Parámetros de la Condición Inicial (Gaussian Wave Packet)
    real(dp), parameter :: x0 = 0.0_dp
    real(dp), parameter :: sigma = 2.0_dp
    real(dp), parameter :: p0 = 5.0_dp
    
    ! Parámetros del Solver Numérico
    real(dp), parameter :: tol_gs = 1.0e-10_dp
    integer, parameter  :: max_iter_gs = 50000
end module params

! ==============================================================================
module physics_core
    use params
    implicit none
    
    ! Vectores Espaciales
    real(dp), allocatable :: x1(:)    ! Nodos enteros (para psi1)
    real(dp), allocatable :: x2(:)    ! Nodos semi-enteros (para psi2)
    
    ! Estado del Espinor Cuántico
    complex(dp), allocatable :: psi1(:)   ! Tamaño N
    complex(dp), allocatable :: psi2(:)   ! Tamaño N-1
    
contains

    subroutine allocate_and_initialize()
        integer :: j
        real(dp) :: E_kin, norm_factor, sum_prob
        complex(dp) :: C1, C2, phase_factor
        real(dp) :: envelope
        
        allocate(x1(N), x2(N-1))
        allocate(psi1(N), psi2(N-1))
        
        ! 1. Inicialización de la Topología de Mallas (Staggered Grid)
        do j = 1, N
            x1(j) = -L + real(j - 1, dp) * dx
        end do
        do j = 1, N - 1
            x2(j) = -L + (real(j, dp) - 0.5_dp) * dx
        end do
        
        ! 2. Álgebra Espinorial Relativista
        E_kin = sqrt(p0**2 + mass**2)
        ! Relación asintótica extraída de la resolución de onda plana
        C1 = cmplx(1.0_dp, 0.0_dp, kind=dp)
        C2 = cmplx(c_light * p0 / (E_kin + mass * c_light**2), 0.0_dp, kind=dp)
        
        ! 3. Construcción del Paquete Gaussiano
        sum_prob = 0.0_dp
        do j = 1, N
            envelope = exp(-((x1(j) - x0)**2) / (4.0_dp * sigma**2)) / ((2.0_dp * acos(-1.0_dp) * sigma**2)**0.25_dp)
            phase_factor = exp(img * p0 * x1(j))
            psi1(j) = C1 * envelope * phase_factor
            sum_prob = sum_prob + abs(psi1(j))**2 * dx
        end do
        
        do j = 1, N - 1
            envelope = exp(-((x2(j) - x0)**2) / (4.0_dp * sigma**2)) / ((2.0_dp * acos(-1.0_dp) * sigma**2)**0.25_dp)
            phase_factor = exp(img * p0 * x2(j))
            psi2(j) = C2 * envelope * phase_factor
            sum_prob = sum_prob + abs(psi2(j))**2 * dx
        end do
        
        ! 4. Normalización Implacable Inicial
        norm_factor = 1.0_dp / sqrt(sum_prob)
        psi1 = psi1 * norm_factor
        psi2 = psi2 * norm_factor
        
        ! Condiciones de Contorno de Dirichlet
        psi1(1) = (0.0_dp, 0.0_dp)
        psi1(N) = (0.0_dp, 0.0_dp)
    end subroutine allocate_and_initialize
    
    subroutine extract_observables(t_current, file_unit)
        real(dp), intent(in) :: t_current
        integer, intent(in) :: file_unit
        
        real(dp) :: norm, pos_x, current_J, mom_p, ener_E
        complex(dp) :: term1, term2
        integer :: j, j_c
        
        ! A. Norma (Conservación Prístina)
        norm = sum(abs(psi1(2:N-1))**2) * dx + sum(abs(psi2(1:N-1))**2) * dx
        
        ! B. Posición Esperada (Baricentro)
        pos_x = sum(x1(2:N-1) * abs(psi1(2:N-1))**2) * dx + sum(x2(1:N-1) * abs(psi2(1:N-1))**2) * dx
        
        ! C. Corriente de Probabilidad J (con interpolación espacial en los sub-nodos)
        current_J = 0.0_dp
        do j = 2, N - 1
            term1 = (psi2(j) + psi2(j-1)) * 0.5_dp
            current_J = current_J + real(conjg(psi1(j)) * term1 + conjg(term1) * psi1(j), dp) * dx
        end do
        
        ! D. Momento Esperado (Diferencias centrales adaptativas)
        mom_p = 0.0_dp
        do j = 2, N - 1
            term1 = -img * (psi1(j+1) - psi1(j-1)) / (2.0_dp * dx)
            mom_p = mom_p + real(conjg(psi1(j)) * term1, dp) * dx
        end do
        do j = 2, N - 2
            term2 = -img * (psi2(j+1) - psi2(j-1)) / (2.0_dp * dx)
            mom_p = mom_p + real(conjg(psi2(j)) * term2, dp) * dx
        end do
        
        ! E. Energía Esperada (Estructura Hermitiana Invariante)
        ener_E = 0.0_dp
        do j = 2, N - 1
            term1 = mass * psi1(j) - img * (psi2(j) - psi2(j-1)) / dx
            ener_E = ener_E + real(conjg(psi1(j)) * term1, dp) * dx
        end do
        do j = 1, N - 1
            term2 = -img * (psi1(j+1) - psi1(j)) / dx - mass * psi2(j)
            ener_E = ener_E + real(conjg(psi2(j)) * term2, dp) * dx
        end do
        
        ! F. Extracción Local de Nodos Centrales (x ~ 0)
        j_c = N / 2
        
        ! Exportación Secuencial al Documento .dat
        write(file_unit, '(F10.4, 11(ES16.7))') &
            t_current, norm, pos_x, current_J, mom_p, ener_E, &
            real(psi1(j_c), dp), aimag(psi1(j_c)), abs(psi1(j_c))**2, &
            real(psi2(j_c), dp), aimag(psi2(j_c)), abs(psi2(j_c))**2
            
    end subroutine extract_observables

end module physics_core

! ==============================================================================
module solvers
    use params
    use physics_core
    implicit none

contains

    subroutine crank_nicolson_step()
        complex(dp), allocatable :: b1(:), b2(:)
        complex(dp), allocatable :: psi1_old(:), psi2_old(:)
        complex(dp) :: A_diag1, A_diag2
        real(dp) :: diff_norm, tau_dx
        integer :: iter
        
        allocate(b1(N), b2(N-1))
        allocate(psi1_old(N), psi2_old(N-1))
        
        tau_dx = tau / dx
        A_diag1 = 1.0_dp + img * tau * mass
        A_diag2 = 1.0_dp - img * tau * mass
        
        ! 1. Evaluación del Vector Derecho (RHS - Explicit Part)
        b1 = (0.0_dp, 0.0_dp)
        b2 = (0.0_dp, 0.0_dp)
        
        ! Bucle implícito (Optimizado para SIMD Cache-lines)
        b1(2:N-1) = (1.0_dp - img * tau * mass) * psi1(2:N-1) - tau_dx * (psi2(2:N-1) - psi2(1:N-2))
        b2(1:N-1) = (1.0_dp + img * tau * mass) * psi2(1:N-1) - tau_dx * (psi1(2:N) - psi1(1:N-1))
        
        ! 2. Inversión Iterativa: Block Gauss-Seidel 
        ! (Vectorizable matemáticamente debido a la gráfica bipartita de la malla escalonada)
        do iter = 1, max_iter_gs
            psi1_old = psi1
            psi2_old = psi2
            
            ! Bloque A: Actualizar submalla psi1 usando old psi2 (sin dependencia recursiva inter-variable)
            psi1(2:N-1) = (b1(2:N-1) - tau_dx * (psi2(2:N-1) - psi2(1:N-2))) / A_diag1
            
            ! Bloque B: Actualizar submalla psi2 usando NEW psi1
            psi2(1:N-1) = (b2(1:N-1) - tau_dx * (psi1(2:N) - psi1(1:N-1))) / A_diag2
            
            ! Condiciones de Contorno Estrictas (Dirichlet)
            psi1(1) = (0.0_dp, 0.0_dp)
            psi1(N) = (0.0_dp, 0.0_dp)
            
            ! 3. Barrera de Tolerancia Espectral (Norma del residual dinámico)
            diff_norm = sqrt(sum(abs(psi1 - psi1_old)**2) + sum(abs(psi2 - psi2_old)**2))
            
            if (diff_norm <= tol_gs) exit
        end do
        
        if (iter >= max_iter_gs) then
            print *, "CRÍTICO: Gauss-Seidel no convergió en el paso temporal. diff_norm =", diff_norm
            stop 1
        end if
        
        deallocate(b1, b2, psi1_old, psi2_old)
    end subroutine crank_nicolson_step

end module solvers

! ==============================================================================
program dirac_1d_simulation
    use params
    use physics_core
    use solvers
    implicit none
    
    integer :: step, file_id
    real(dp) :: t_current
    
    ! Apertura y formato del documento sistemático de observación
    file_id = 10
    open(unit=file_id, file='dirac_observables.dat', status='replace', action='write')
    
    ! Cabecera Estructural
    write(file_id, '(A)') '# 1:t 2:Norma 3:PosX_Esp 4:CorrienteJ 5:Momento_Esp 6:Energia_Esp 7:Re(psi1_cnt) 8:Im(psi1_cnt) 9:Mod2(psi1_cnt) 10:Re(psi2_cnt) 11:Im(psi2_cnt) 12:Mod2(psi2_cnt)'
    
    call allocate_and_initialize()
    
    ! Estado Primitivo (t = 0)
    t_current = 0.0_dp
    call extract_observables(t_current, file_id)
    
    ! Integración Global Temporal
    print *, "Iniciando Integración de Crank-Nicolson..."
    do step = 1, Nt
        call crank_nicolson_step()
        
        t_current = real(step, dp) * dt
        call extract_observables(t_current, file_id)
        
        if (mod(step, 100) == 0) then
            print *, "Paso completado:", step, "/", Nt, " - Tiempo:", t_current
        end if
    end do
    
    close(file_id)
    print *, "Simulación Finalizada. Datos exportados exitosamente en dirac_observables.dat"
    
end program dirac_1d_simulation

