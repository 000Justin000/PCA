      program ppca
c
c     Message Passing Layer: MPI
c
c     Simple program to illustrate the idea of reverse communication
c     in regular mode for a standard symmetric eigenvalue problem.
c
c     We implement example one of ex-sym.doc in SRC directory
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is derived from the central difference discretization
c         of the 2-dimensional Laplacian on the unit square with
c         zero Dirichlet boundary condition.
c     ... OP = A  and  B = I.
c     ... Assume "call av (n,x,y)" computes y = A*x.
c     ... Use mode 1 of DSAUPD.
c
c\BeginLib
c
c\Routines called:
c     pdsaupd  Parallel ARPACK reverse communication interface routine.
c     pdseupd  Parallel ARPACK routine that returns Ritz values and (optionally)
c              Ritz vectors.
c     pdnorm2  Parallel version of Level 1 BLAS that computes the norm of a vector.
c     daxpy    Level 1 BLAS that computes y <- alpha*x+y.
c     av       Matrix vector multiplication routine that computes A*x.
c
c\Author
c     Richard Lehoucq
c     Danny Sorensen
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c     Junteng Jia
c
c\Parallel Modifications
c     Kristi Maschhoff
c
c\Revision history:
c     Starting Point: Serial Code FILE: sdrv1.F   SID: 2.2
c
c FILE: sdrv1.F   SID: 1.4   DATE OF SID: 3/19/97   RELEASE: 1
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      include 'mpif.h'
      include 'debug.h'
      include 'stat.h'
 
c     %---------------%
c     | MPI INTERFACE |
c     %---------------%
 
      integer           comm, myid, nprocs, rc, nloc
c
c     %-------------------------------%
c     | Define leading dimensions     |
c     | for all arrays.               |
c     | MAXNLOC: Maximum dimension of |
c     |          the local A allowed. |
c     | MAXNEV:  Maximum NEV allowed  |
c     | MAXNCV:  Maximum NCV allowed  |
c     %-------------------------------%
c
      integer          maxnloc, maxnev, maxncv, ldv
      parameter       (maxnloc=10000, maxnev=1000, maxncv=1100, 
     &                 ldv=maxnloc )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      double precision
     &                 v(ldv,maxncv), workl(maxncv*(maxncv+8)),
     &                 workd(3*maxnloc), d(maxncv,2), resid(maxnloc),
     &                 ax(maxnloc)
      logical          select(maxncv)
      integer          iparam(11), ipntr(11)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info, ierr, j, 
     &                 nx, nconv, maxitr, mode, ishfts
      logical          rvec
      double precision      
     &                 tol, sigma
c
c     %----------------------------------------------%
c     | Local Buffers needed for MPI communication   |
c     %----------------------------------------------%
c
      double precision
     &                  mpi_buf(300000)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                 zero
      parameter        (zero = 0.0)
c  
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      double precision           
     &                 pdnorm2 
      external         pdnorm2, daxpy
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic         abs

c
c     %-------------------%
c     | Covariance Matrix |
c     %-------------------%
c
      real(kind=8), allocatable ::  covmat_loc(:,:)
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      call MPI_INIT( ierr )
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
c
      ndigit = -3
      logfil = 6
      msaupd = 1
c
c     %-------------------------------------------------------%
c     | The number N(=NX*NX) is the dimension of the matrix.  |
c     | A standard eigenvalue problem is solved (BMAT = 'I'). | 
c     | NEV is the number of eigenvalues to be approximated.  |
c     | The user can modify NEV, NCV, WHICH to solve problems |
c     | of different sizes, and to get different parts of the |
c     | spectrum.  However, The following conditions must be  |
c     | satisfied:                                            |
c     |                   N <= MAXN,                          | 
c     |                 NEV <= MAXNEV,                        |
c     |             NEV + 2 <= NCV <= MAXNCV                  | 
c     %-------------------------------------------------------% 
c
      n   =   625     ! 107207
      nev =   100     ! 1000
      ncv =   200     ! 1100

c
c     %--------------------------------------%
c     | Set up distribution of data to nodes |
c     %--------------------------------------%
c
      nloc = n / nprocs
      if ( mod(n, nprocs) .ne. 0 ) then
         nloc = nloc + 1
        
         if ( myid .eq. nprocs - 1 ) then
            nloc = mod(n, nloc)
         end if
      end if
c
      if ( nloc .gt. maxnloc ) then
         print *, ' ERROR with _SDRV1: NLOC is greater than MAXNLOC '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'I'
      which = 'LM'
c
      allocate( covmat_loc(n, nloc) )
      call read_cov( comm, n, nloc, covmat_loc )
c
c     %--------------------------------------------------%
c     | The work array WORKL is used in PSSAUPD as       |
c     | workspace.  Its dimension LWORKL is set as       |
c     | illustrated below.  The parameter TOL determines |
c     | the stopping criterion.  If TOL<=0, machine      |
c     | precision is used.  The variable IDO is used for |
c     | reverse communication and is initially set to 0. |
c     | Setting INFO=0 indicates that a random vector is |
c     | generated in PSSAUPD to start the Arnoldi        |
c     | iteration.                                       |
c     %--------------------------------------------------%
c
      lworkl = ncv*(ncv+8)
      tol = zero 
      info = 0
      ido = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of PSSAUPD is used    |
c     | (IPARAM(7) = 1).  All these options may be        |
c     | changed by the user. For details, see the         |
c     | documentation in PSSAUPD.                         |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 1000
      mode   = 1
c      
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 
c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) |
c     %-------------------------------------------%
c
c

 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine PSSAUPD and take| 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call pdsaupd ( comm, ido, bmat, nloc, which, nev, tol, resid, 
     &                  ncv, v, ldv, iparam, ipntr, workd, workl,
     &                  lworkl, info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %--------------------------------------%
c           | Perform matrix vector multiplication |
c           |              y <--- OP*x             |
c           | The user should supply his/her own   |
c           | matrix vector multiplication routine |
c           | here that takes workd(ipntr(1)) as   |
c           | the input, and return the result to  |
c           | workd(ipntr(2)).                     |
c           %--------------------------------------%
c
            call av ( comm, n, nloc, covmat_loc, mpi_buf, 
     &                workd(ipntr(1)), workd(ipntr(2)) )
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call PSSAUPD again.|
c           %-----------------------------------------%
c
            go to 10
c
         end if 
c
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in PSSAUPD.|
c        %--------------------------%
c
         if ( myid .eq. 0 ) then
            print *, ' '
            print *, ' Error with _saupd, info = ', info
            print *, ' Check documentation in _saupd '
            print *, ' '
         endif
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using PSSEUPD.               |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |  
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        %-------------------------------------------%
c           
         rvec = .true.
c
         call pdseupd ( comm, rvec, 'All', select, d, v, 
     &        ldv, sigma, bmat, nloc, which, nev, tol, 
     &        resid, ncv, v, ldv, iparam, ipntr, workd, 
     &        workl, lworkl, ierr )
c        %----------------------------------------------%
c        | Eigenvalues are returned in the first column |
c        | of the two dimensional array D and the       |
c        | corresponding eigenvectors are returned in   |
c        | the first NEV columns of the two dimensional |
c        | array V if requested.  Otherwise, an         |
c        | orthogonal basis for the invariant subspace  |
c        | corresponding to the eigenvalues in D is     |
c        | returned in V.                               |
c        %----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of PSSEUPD.|
c            %------------------------------------%
c
c
            if ( myid .eq. 0 ) then
                print *, ' '
                print *, ' Error with _seupd, info = ', ierr
                print *, ' Check the documentation of _seupd. '
                print *, ' '
            endif
c
         else
c
             nconv =  iparam(5)
c
             do 20 j=1, nconv
c
c               %---------------------------%
c               | Compute the residual norm |
c               |                           |
c               |   ||  A*x - lambda*x ||   |
c               |                           |
c               | for the NCONV accurately  |
c               | computed eigenvalues and  |
c               | eigenvectors.  (iparam(5) |
c               | indicates how many are    |
c               | accurate to the requested |
c               | tolerance)                |
c               %---------------------------%
c
                call av(comm, n, nloc, covmat_loc, mpi_buf, v(1, j), ax)
                call daxpy(nloc, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = pdnorm2( comm, nloc, ax, 1 )
                if (myid .eq. 0) 
     &              write (*,*) "eigenvalue No.", j, "=", d(j,1)
c
 20          continue
c
c            %-------------------------------%
c            | Display computed residuals    |
c            %-------------------------------%
c
             call pdmout(comm, 6, nconv, 2, d, maxncv, -6,
     &            'Ritz values and direct residuals')
         end if
c
c        %------------------------------------------%
c        | Print additional convergence information |
c        %------------------------------------------%
c
         if (myid .eq. 0)then
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' ' 
            print *, ' No shifts could be applied during implicit
     &                 Arnoldi update, try increasing NCV.'
            print *, ' '
         end if      
c
         print *, ' '
         print *, '_SDRV1 '
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of processes is ', nprocs
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', 
     &              nconv 
         print *, ' The number of Implicit Arnoldi update',
     &            ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
         endif
c
      end if
c
c     %---------------------------%
c     |  Done with program ppca.  |
c     %---------------------------%
c
 9000 continue
c
      call MPI_FINALIZE(rc)
c
      end
c 
c=========================================================================
c ------------------------------------------------------------------
c     parallel matrix vector subroutine
c
c     The matrix used is the 2 dimensional discrete Laplacian on unit
c     square with zero Dirichlet boundary condition.
c
c     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block
c     tridiagonal matrix
c
c-------------------------------------------------------------------
c
      subroutine av (comm, n, nloc, covmat_loc, mpi_buf, v, w)
c
c     .. MPI Declarations ...
      include           'mpif.h'
      integer           comm, nprocs, myid, ierr
      integer           nloc, itr, gather_size
      double precision
     &                  v(nloc), w(nloc), mpi_buf(300000),
     &                  covmat_loc(n, nloc), sumv
c
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
c
      gather_size = n / nprocs
      if ( mod(n, nprocs) .ne. 0 ) then
         gather_size = gather_size + 1
      end if
c
      call MPI_ALLGATHER( v, gather_size, MPI_DOUBLE, mpi_buf, 
     &                    gather_size, MPI_DOUBLE, comm, ierr)
c
      DO itr = 1, nloc
        w(itr) = dot_product(covmat_loc(:,itr),mpi_buf(1:n)) 
      END DO
    
      return
      end
c=========================================================================

c=========================================================================
c ------------------------------------------------------------------
c     parallel read covariance matrix
c-------------------------------------------------------------------
c
      subroutine read_cov (comm, n, nloc, covmat_loc)
c
c     .. MPI Declarations ...
c
      include           'mpif.h'
      integer           comm, myid, nprocs, ierr
      integer           n, nloc
      integer           itr
      double precision
     &                  covmat_loc(n, nloc)
      integer           my_unit
      character         my_file*20

c
      call MPI_COMM_RANK( comm, myid, ierr )
c
      WRITE (my_file, "(A, I1)") "covmat", myid
c      WRITE (my_file, "(A, I1)") "result", myid

      my_unit = myid
      
      OPEN (UNIT=my_unit, FILE=my_file)
      
      READ (my_unit, *) covmat_loc

      CLOSE (my_unit)

      return
      end
c=========================================================================
