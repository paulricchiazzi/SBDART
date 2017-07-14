c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: /home/paul/rt/sbrt/RCS/disutil.f,v 1.7 2002/08/14 22:56:35 paul Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c ---------------------------------------------------------------------
c  Fortran-90 versions of machine-constant routines R1MACH, D1MACH, I1MACH
c
c  {R,D,I}1MACH revisited: no more uncommenting DATA statements
c  
c  Presented at the IFIP WG 2.5 International Workshop on 
c  "Current Directions in Numerical Software and High Performance 
c  Computing", 19 - 20 October 1995, Kyoto, Japan. 
c  
c  The widely-used original routines were modified to use Fortran-90 
c  intrinsic functions.  This was not completely possible with I1MACH, 
c  which returns some parameters (logical unit numbers of standard
c  input, standard output, and standard error) that may require
c  user customization. 
c  
c  David Gay (dmg@bell-labs.com)
c  Eric Grosse (ehg@bell-labs.com)
c  Bell Laboratories
c  700 Mountain Avenue
c  Murray Hill, New Jersey 07974-0636
c  USA 
c  
c  References:
c  
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996.
c
c http://www.nsc.liu.se/~boein/ifip/kyoto/workshop-info/proceedings/einarsson
c    /d1mach.html  (THIS WEB SITE WORKED AS OF APR 2000)
c -------------------------------------------------------------------------


      FUNCTION R1MACH (I)
c
c   R1MACH can be used to obtain machine-dependent parameters for
c   single precision numbers.  The results for various values of I are:
c
c   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c   R1MACH(3) = B**(-T), the smallest relative spacing.
c   R1MACH(4) = B**(1-T), the largest relative spacing.
c   R1MACH(5) = LOG10(B)
c
c   Assume single precision numbers are represented in the T-digit,
c   base-B form
c
c              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c   where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c   The values of B, T, EMIN and EMAX are provided in I1MACH as follows:
c   I1MACH(10) = B, the base.
c   I1MACH(11) = T, the number of base-B digits.
c   I1MACH(12) = EMIN, the smallest exponent E.
c   I1MACH(13) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c     Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c     August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   790101  DATE WRITTEN
c   960329  Modified for Fortran 90 (BE after suggestions by Eric Grosse)      
c --------------------------------------------------------------------

      use params, only: kr
      IMPLICIT NONE
      INTEGER :: I
      REAL(KR) :: B, X = 1.0, r1mach

      B = RADIX(X)

      SELECT CASE (I)
        CASE (1) 
          R1MACH = TINY(X)            ! smallest positive magnitude.
        CASE (2)
          R1MACH = HUGE(X)            ! largest magnitude.
        CASE (3)
          R1MACH = B**(-DIGITS(X))    ! smallest relative spacing.
        CASE (4)
          R1MACH = B**(1-DIGITS(X))   ! largest relative spacing.
        CASE (5)
          R1MACH = LOG10(B)
        CASE DEFAULT
          STOP 'R1MACH -- input argument out of bounds'
      END SELECT

      RETURN
      END FUNCTION R1MACH


      DOUBLE PRECISION FUNCTION D1MACH (I)
c
c   D1MACH can be used to obtain machine-dependent parameters for
c   double precision numbers.  The results for various values of I are:
c
c   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c   D1MACH(3) = B**(-T), the smallest relative spacing.
c   D1MACH(4) = B**(1-T), the largest relative spacing.
c   D1MACH(5) = LOG10(B)
c
c   Assume double precision numbers are represented in the T-digit,
c   base-B form
c
c        sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c   where 0 <= X(I) < B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c   The values of B, T, EMIN and EMAX are provided in I1MACH as follows:
c   I1MACH(10) = B, the base.
c   I1MACH(11) = T, the number of base-B digits.
c   I1MACH(12) = EMIN, the smallest exponent E.
c   I1MACH(13) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   790101  DATE WRITTEN
c   960329  Modified for Fortran 90 (BE after suggestions by Eric Grosse)      
c --------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER :: I
      DOUBLE PRECISION :: B, X = 1.D0

      B = RADIX(X)

      SELECT CASE (I)
        CASE (1)
          D1MACH = TINY(X)            ! smallest positive magnitude.
        CASE (2)
          D1MACH = HUGE(X)            ! largest magnitude.
        CASE (3)
          D1MACH = B**(-DIGITS(X))    ! smallest relative spacing.
        CASE (4)
          D1MACH = B**(1-DIGITS(X))   ! largest relative spacing.
        CASE (5)
          D1MACH = LOG10(B)
        CASE DEFAULT
          STOP 'D1MACH -- input arg out of bounds'
      END SELECT

      RETURN
      END FUNCTION D1MACH


      INTEGER FUNCTION I1MACH (I)
c
c   I1MACH can be used to obtain machine-dependent parameters for the
c   local machine environment.  The results for various values of I are:
c
c   I/O unit numbers (**MAY REQUIRE USER CUSTOMIZATION**):
c     I1MACH( 1) = the standard input unit.
c     I1MACH( 2) = the standard output unit.
c     I1MACH( 3) = the standard punch unit (obsolete, will cause error)
c     I1MACH( 4) = the standard error message unit.
c                  (the error message unit is usually 0 in UNIX systems)
c
c   Words:
c     I1MACH( 5) = the number of bits per integer storage unit.
c     I1MACH( 6) = the number of characters per integer storage unit.
c                  (obsolete, will cause an error)
c
c   Integers:
c     assume integers are represented in the S-digit, base-A form
c
c          sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
c
c     where 0 <= X(I) < A for I=0,...,S-1.
c
c     I1MACH( 7) = A, the base.
c     I1MACH( 8) = S, the number of base-A digits.
c     I1MACH( 9) = A**S - 1, the largest magnitude.
c
c   Floating-Point Numbers:
c     Assume floating-point numbers are represented in the T-digit,
c     base-B form
c                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
c
c     where 0 <= X(I) .LT. B for I=1,...,T; 0 < X(1); and EMIN <= E <= EMAX.
c
c     I1MACH(10) = B, the base.
c
c   Single-Precision:
c     I1MACH(11) = T, the number of base-B digits.
c     I1MACH(12) = EMIN, the smallest exponent E.
c     I1MACH(13) = EMAX, the largest exponent E.
c
c   Double-Precision:
c     I1MACH(14) = T, the number of base-B digits.
c     I1MACH(15) = EMIN, the smallest exponent E.
c     I1MACH(16) = EMAX, the largest exponent E.
c
c***REFERENCES  
c
c  P. Fox, A. Hall and N. Schryer, Framework for a portable library,
c     ACM Transactions on Mathematical Software 4, 177-188 (1978).
c
c  David Gay and Eric Grosse, Comment on Algorithm 528, Bell Labs, Murray 
c    Hill, NJ. submitted to ACM Transactions on Mathematical Software,
c    August 1996. 
c
c***REVISION HISTORY  (YYMMDD)
c   750101  DATE WRITTEN
c   960411  Modified for Fortran 90 (BE after suggestions by Eric Grosse)    
c --------------------------------------------------------------------

      use params, only: kr
      IMPLICIT NONE
      INTEGER :: I
      REAL(KR) :: X_single  = 1.0
      DOUBLE PRECISION :: X_double = 1.D0

      SELECT CASE (I)
        CASE (1)
          I1MACH = 5 ! Input unit
        CASE (2)
          I1MACH = 6 ! Output unit
        CASE (3)
          STOP 'I1MACH: input arg = 3 is obsolete'
        CASE (4)
          I1MACH = 0 ! Error message unit
        CASE (5)
          I1MACH = BIT_SIZE(I)
        CASE (6)
          STOP 'I1MACH: input arg = 6 is obsolete'
        CASE (7)
          I1MACH = RADIX(1)
        CASE (8)
          I1MACH = BIT_SIZE(I) - 1
        CASE (9)
          I1MACH = HUGE(1)
        CASE (10)
          I1MACH = RADIX(X_single)
        CASE (11)
          I1MACH = DIGITS(X_single)
        CASE (12)
          I1MACH = MINEXPONENT(X_single)
        CASE (13)
          I1MACH = MAXEXPONENT(X_single)
        CASE (14)
          I1MACH = DIGITS(X_double)
        CASE (15)
          I1MACH = MINEXPONENT(X_double)
        CASE (16)
          I1MACH = MAXEXPONENT(X_double) 
        CASE DEFAULT
          STOP 'I1MACH: input argument out of bounds'
      END SELECT

      RETURN
      END FUNCTION I1MACH

c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: /home/paul/rt/sbrt/RCS/disutil.f,v 1.7 2002/08/14 22:56:35 paul Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine  errmsg(msgnum,messag)

c        print out a warning or error messages to file disort_warning.??
c
c input:
c  msgnum    message number, used to identify whether a given message
c            has been issued more than once (0<msgnum<=nmsg)
c            or to idicate a fatal error (msgnum=0)
c
c  messag    a warning or fatal error message (character string)

      implicit none
      integer, parameter :: nmsg=20
      character*(*) messag
      character (len=2)  :: num
      character (len=132) :: line
      integer ::    msgnum
      integer,save :: msgset(nmsg)=0
      logical,save :: first=.true.

      if(msgnum.eq.0) then
        line='ERROR  >>>>>>'
      else
        if(msgset(msgnum).eq.1) return
        line='WARNING >>>>>'
      endif

      write(num,'(i2.2)') msgnum

      open(16,file='SBDART_WARNING.'//num,status='unknown',
     &     form='formatted')    
  
      write(16,'(a,1x,a)') trim(line),messag
      write(16,'(/70("#")/)')

      do 
        read(11,'(a)',end=100) line
        write(16,'(a)') trim(line)
      enddo
 100  continue

      rewind 11
      close(16)

      if(msgnum.eq.0) stop
      msgset(msgnum)=1
      return
      end

      LOGICAL FUNCTION  WrtBad ( VarNam )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

      implicit none
      CHARACTER*(*)  VarNam
      INTEGER        MaxMsg, NumMsg
      SAVE  NumMsg, MaxMsg
      DATA  NumMsg / 0 /,  MaxMsg / 50 /


      WrtBad = .TRUE.
      NumMsg = NumMsg + 1
      WRITE ( *, '(3A)' )  ' ****  Input variable  ', VarNam,
     &                     '  in error  ****'
      IF ( NumMsg.EQ.MaxMsg )
     &   call  errmsg (12,'Too many input errors.  Aborting...')

      RETURN
      END

      LOGICAL FUNCTION  WrtDim ( DimNam, MinVal )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

      implicit none
      CHARACTER*(*)  DimNam
      INTEGER        MinVal


      WRITE ( *, '(/,3A,I7)' )  ' ****  Symbolic dimension  ', DimNam,
     &                     '  should be increased to at least ', MinVal
      WrtDim = .TRUE.

      RETURN
      END

      LOGICAL FUNCTION  TstBad( VarNam, RelErr )

c       Write name (VarNam) of variable failing self-test and its
c       percent error from the correct value;  return  'FALSE'.

      use params, only: kr
      CHARACTER*(*)  VarNam
      REAL(KR)           RelErr


      TstBad = .FALSE.
      WRITE( *, '(/,3A,1P,E11.2,A)' )
     &       ' Output variable ', VarNam,' differed by ', 100.*RelErr,
     &       ' per cent from correct value.  Self-test failed.'

      RETURN
      END

c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c RCS version control information:
c $Header: /home/paul/rt/sbrt/RCS/disutil.f,v 1.7 2002/08/14 22:56:35 paul Exp $
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c Call tree:
c
c    SGBCO
c       SASUM
c       SDOT
c       SAXPY
c       SGBFA
c           ISAMAX
c           SAXPY
c           SSCAL
c       SSCAL
c   SGBSL
c       SDOT
c       SAXPY
c   SGECO
c       SASUM
c       SDOT
c       SAXPY
c       SGEFA
c           ISAMAX
c           SAXPY
c           SSCAL
c       SSCAL
c   SGESL
c       SDOT
c       SAXPY
c   SSWAP
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      SUBROUTINE SGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )

c         Factors a real band matrix by Gaussian elimination
c         and estimates the condition of the matrix.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     If  RCOND  is not needed, SGBFA is slightly faster.
c     To solve  A*X = B , follow SBGCO by SGBSL.

c     input:

C        ABD     REAL(LDA, N)
c                contains the matrix in band storage.  The columns
c                of the matrix are stored in the columns of  ABD  and
c                the diagonals of the matrix are stored in rows
c                ML+1 through 2*ML+MU+1 of  ABD .
c                See the comments below for details.

C        LDA     INTEGER
c                the leading dimension of the array  ABD .
c                LDA must be .GE. 2*ML + MU + 1 .

C        N       INTEGER
c                the order of the original matrix.

C        ML      INTEGER
c                number of diagonals below the main diagonal.
c                0 .LE. ML .LT. N .

C        MU      INTEGER
c                number of diagonals above the main diagonal.
c                0 .LE. MU .LT. N .
c                more efficient if  ML .LE. MU .

c     on return

c        ABD     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                The factorization can be written  A = L*U  where
c                L  is a product of permutation and unit lower
c                triangular matrices and  U  is upper triangular.

C        IPVT    INTEGER(N)
c                an integer vector of pivot indices.

C        RCOND   REAL
c                an estimate of the reciprocal condition of  A .
c                For the system  A*X = B , relative perturbations
c                in  A  and  B  of size  epsilon  may cause
c                relative perturbations in  X  of size  epsilon/RCOND .
c                If  RCOND  is so small that the logical expression
c                           1.0 + RCOND .EQ. 1.0
c                is true, then  A  may be singular to working
c                precision.  In particular,  RCOND  is zero  if
c                exact singularity is detected or the estimate
c                underflows.

C        Z       REAL(N)
c                a work vector whose contents are usually unimportant.
c                If  A  is close to a singular matrix, then  Z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .

c     Band storage

c           If  A  is a band matrix, the following program segment
c           will set up the input.

c                   ML = (band width below the diagonal)
c                   MU = (band width above the diagonal)
c                   M = ML + MU + 1
c                   DO 20 J = 1, N
c                      I1 = MAX(1, J-MU)
c                      I2 = MIN(N, J+ML)
c                      DO 10 I = I1, I2
c                         K = I - J + M
c                         ABD(K,J) = A(I,J)
c                10    CONTINUE
c                20 CONTINUE

c           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
c           In addition, the first  ML  rows in  ABD  are used for
c           elements generated during the triangularization.
c           The total number of rows needed in  ABD  is  2*ML+MU+1 .
c           The  ML+MU by ML+MU  upper left triangle and the
c           ML by ML  lower right triangle are not referenced.

c     Example:  if the original matrix is

c           11 12 13  0  0  0
c           21 22 23 24  0  0
c            0 32 33 34 35  0
c            0  0 43 44 45 46
c            0  0  0 54 55 56
c            0  0  0  0 65 66

c      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABD should contain

c            *  *  *  +  +  +  , * = not used
c            *  * 13 24 35 46  , + = used for pivoting
c            * 12 23 34 45 56
c           11 22 33 44 55 66
c           21 32 43 54 65  *

c --------------------------------------------------------------------


c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   LDA, ML, MU, N
      REAL(KR)      RCOND
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL(KR)      ABD( LDA, * ), Z( * )
c     ..
c     .. Local Scalars ..

      INTEGER   INFO, IS, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM
      REAL(KR)      ANORM, EK, S, SM, T, WK, WKM, YNORM
c     ..
c     .. External Functions ..

      REAL(KR)      SASUM, SDOT
      EXTERNAL  SASUM, SDOT
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY, SGBFA, SSCAL
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, MAX, MIN, SIGN
c     ..


c                       ** compute 1-norm of A
      ANORM  = 0.0E0
      L  = ML + 1
      IS = L + MU

      DO 10 J = 1, N

         ANORM  = MAX( ANORM, SASUM( L,ABD( IS,J ),1 ) )

         IF( IS.GT.ML + 1 ) IS = IS - 1

         IF( J.LE.MU ) L  = L + 1

         IF( J.GE.N - ML ) L  = L - 1

   10 CONTINUE
c                                               ** factor

      CALL SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

c     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
c     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E.
c     trans(A) is the transpose of A.  The components of E  are
c     chosen to cause maximum local growth in the elements of W  where
c     trans(U)*W = E.  The vectors are frequently rescaled to avoid
c     overflow.

c                     ** solve trans(U)*W = E
      EK = 1.0E0

      DO 20 J = 1, N
         Z( J ) = 0.0E0
   20 CONTINUE


      M  = ML + MU + 1
      JU = 0

      DO 50 K = 1, N

         IF( Z( K ).NE.0.0E0 ) EK = SIGN( EK, -Z( K ) )

         IF( ABS( EK - Z( K ) ).GT.ABS( ABD( M,K ) ) ) THEN

            S  = ABS( ABD( M,K ) ) / ABS( EK - Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            EK = S*EK

         END IF

         WK   = EK - Z( K )
         WKM  = -EK - Z( K )
         S    = ABS( WK )
         SM   = ABS( WKM )

         IF( ABD( M,K ).NE.0.0E0 ) THEN

            WK   = WK / ABD( M, K )
            WKM  = WKM / ABD( M, K )

         ELSE

            WK   = 1.0E0
            WKM  = 1.0E0

         END IF

         KP1  = K + 1
         JU   = MIN( MAX( JU,MU + IPVT( K ) ), N )
         MM   = M

         IF( KP1.LE.JU ) THEN

            DO 30 J = KP1, JU
               MM     = MM - 1
               SM     = SM + ABS( Z( J ) + WKM*ABD( MM,J ) )
               Z( J ) = Z( J ) + WK*ABD( MM, J )
               S      = S + ABS( Z( J ) )
   30       CONTINUE

            IF( S.LT.SM ) THEN

               T  = WKM - WK
               WK = WKM
               MM = M

               DO 40 J = KP1, JU
                  MM = MM - 1
                  Z( J ) = Z( J ) + T*ABD( MM, J )
   40          CONTINUE

            END IF

         END IF

         Z( K ) = WK

   50 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

c                         ** solve trans(L)*Y = W
      DO 60 KB = 1, N
         K  = N + 1 - KB
         LM = MIN( ML, N - K )

         IF( K.LT.N )
     &       Z( K ) = Z( K ) + SDOT( LM, ABD( M+1, K ), 1, Z( K+1 ), 1 )

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

         END IF

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T

   60 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = 1.0E0
c                         ** solve L*V = Y
      DO 70 K = 1, N

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
         LM     = MIN( ML, N - K )

         IF( K.LT.N )
     &       CALL SAXPY( LM, T, ABD( M+1, K ), 1, Z( K+1 ), 1 )

         IF( ABS( Z(K) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z(K) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

   70 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = S*YNORM

c                           ** solve  U*Z = W
      DO 80 KB = 1, N

         K  = N + 1 - KB

         IF( ABS( Z( K ) ).GT.ABS( ABD( M,K ) ) ) THEN

            S  = ABS( ABD( M,K ) ) / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

         IF( ABD( M,K ).NE.0.0E0 ) Z( K ) = Z( K ) / ABD( M, K )
         IF( ABD( M,K ).EQ.0.0E0 ) Z( K ) = 1.0E0

         LM = MIN( K, M ) - 1
         LA = M - LM
         LZ = K - LM
         T  = -Z( K )

         CALL SAXPY( LM, T, ABD( LA,K ), 1, Z( LZ ), 1 )

   80 CONTINUE
c                              ** make znorm = 1.0

      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = S*YNORM
      IF( ANORM.NE.0.0E0 ) RCOND  = YNORM / ANORM
      IF( ANORM.EQ.0.0E0 ) RCOND  = 0.0E0

      END

      SUBROUTINE SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

c         Factors a real band matrix by elimination.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     SGBFA is usually called by SBGCO, but it can be called
c     directly with a saving in time if  RCOND  is not needed.

c     Input:  same as SGBCO

c     On return:

c        ABD,IPVT    same as SGBCO

c        INFO    INTEGER
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  This is not an error
c                     condition for this subroutine, but it does
c                     indicate that SGBSL will divide by zero if
c                     called.  Use  RCOND  in SBGCO for a reliable
c                     indication of singularity.

c     (see SGBCO for description of band storage mode)

c ----------------------------------------------------------------


c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   INFO, LDA, ML, MU, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL(KR)      ABD( LDA, * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, I0, J, J0, J1, JU, JZ, K, KP1, L, LM, M, MM, NM1
      REAL(KR)      T
c     ..
c     .. External Functions ..

      INTEGER   ISAMAX
      EXTERNAL  ISAMAX
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY, SSCAL
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MAX, MIN
c     ..


      M    = ML + MU + 1
      INFO = 0
c                        ** zero initial fill-in columns
      J0 = MU + 2
      J1 = MIN( N, M ) - 1

      DO 20 JZ = J0, J1

         I0 = M + 1 - JZ

         DO 10 I = I0, ML
            ABD( I, JZ ) = 0.0E0
   10    CONTINUE

   20 CONTINUE

      JZ = J1
      JU = 0
c                       ** Gaussian elimination with partial pivoting
      NM1  = N - 1

      DO 50 K = 1, NM1

         KP1 = K + 1
c                                  ** zero next fill-in column
         JZ = JZ + 1

         IF( JZ.LE.N ) THEN

            DO 30 I = 1, ML
               ABD( I, JZ ) = 0.0E0
   30       CONTINUE

         END IF
c                                  ** find L = pivot index
         LM  = MIN( ML, N - K )
         L   = ISAMAX( LM + 1, ABD( M, K ), 1 ) + M - 1
         IPVT( K ) = L + K - M

         IF( ABD( L,K ).EQ.0.0E0 ) THEN
c                                      ** zero pivot implies this column
c                                      ** already triangularized
            INFO = K

         ELSE
c                                ** interchange if necessary
            IF( L.NE.M ) THEN

               T           = ABD( L, K )
               ABD( L, K ) = ABD( M, K )
               ABD( M, K ) = T
            END IF
c                                      ** compute multipliers
            T  = - 1.0E0 / ABD( M, K )

            CALL SSCAL( LM, T, ABD( M + 1,K ), 1 )

c                               ** row elimination with column indexing

            JU = MIN( MAX( JU,MU + IPVT( K ) ), N )
            MM = M

            DO 40 J = KP1, JU

               L  = L - 1
               MM = MM - 1
               T  = ABD( L, J )

               IF( L.NE.MM ) THEN

                  ABD( L, J ) = ABD( MM, J )
                  ABD( MM, J ) = T

               END IF

               CALL SAXPY( LM, T, ABD( M+1, K ), 1, ABD( MM+1, J ), 1)

   40       CONTINUE

         END IF

   50 CONTINUE


      IPVT( N ) = N
      IF( ABD( M,N ).EQ.0.0E0 ) INFO = N

      END

      SUBROUTINE SGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )

c         Solves the real band system
c            A * X = B  or  transpose(A) * X = B
c         using the factors computed by SBGCO or SGBFA.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     Input:

C        ABD     REAL(LDA, N)
c                the output from SBGCO or SGBFA.

C        LDA     INTEGER
c                the leading dimension of the array  ABD .

C        N       INTEGER
c                the order of the original matrix.

C        ML      INTEGER
c                number of diagonals below the main diagonal.

C        MU      INTEGER
c                number of diagonals above the main diagonal.

C        IPVT    INTEGER(N)
c                the pivot vector from SBGCO or SGBFA.

C        B       REAL(N)
c                the right hand side vector.

C        JOB     INTEGER
c                = 0         to solve  A*X = B ,
c                = nonzero   to solve  transpose(A)*X = B

c     On return

c        B       the solution vector  X

c     Error condition

c        A division by zero will occur if the input factor contains a
c        zero on the diagonal.  Technically, this indicates singularity,
c        but it is often caused by improper arguments or improper
c        setting of LDA .  It will not occur if the subroutines are
c        called correctly and if SBGCO has set RCOND .GT. 0.0
c        or SGBFA has set INFO .EQ. 0 .

c     To compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue

c --------------------------------------------------------

c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   JOB, LDA, ML, MU, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL(KR)      ABD( LDA, * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER   K, KB, L, LA, LB, LM, M, NM1
      REAL(KR)      T
c     ..
c     .. External Functions ..

      REAL(KR)      SDOT
      EXTERNAL  SDOT
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..


      M   = MU + ML + 1
      NM1 = N - 1

      IF( JOB.EQ.0 ) THEN
c                           ** solve  A * X = B

c                               ** first solve L*Y = B
         IF( ML.NE.0 ) THEN

            DO 10 K = 1, NM1

               LM = MIN( ML, N - K )
               L  = IPVT( K )
               T  = B( L )

               IF( L.NE.K ) THEN

                  B( L ) = B( K )
                  B( K ) = T

               END IF

               CALL SAXPY( LM, T, ABD( M + 1,K ), 1, B( K + 1 ), 1 )

   10       CONTINUE

         END IF

c                           ** now solve  U*X = Y
         DO 20 KB = 1, N

            K      = N + 1 - KB
            B( K ) = B( K ) / ABD( M, K )
            LM     = MIN( K, M ) - 1
            LA     = M - LM
            LB     = K - LM
            T      = -B( K )

            CALL SAXPY( LM, T, ABD( LA,K ), 1, B( LB ), 1 )

   20    CONTINUE


      ELSE
c                          ** solve  trans(A) * X = B

c                                  ** first solve  trans(U)*Y = B
         DO 30 K = 1, N

            LM     = MIN( K, M ) - 1
            LA     = M - LM
            LB     = K - LM
            T      = SDOT( LM, ABD( LA,K ), 1, B( LB ), 1 )
            B( K ) = ( B( K ) - T ) / ABD( M, K )

   30    CONTINUE

c                                  ** now solve trans(L)*X = Y
         IF( ML.NE.0 ) THEN

            DO 40 KB = 1, NM1

               K      = N - KB
               LM     = MIN( ML, N - K )
               B( K ) = B( K ) + SDOT( LM, ABD( M+1, K ), 1,
     &                                 B( K+1 ), 1 )
               L      = IPVT( K )

               IF( L.NE.K ) THEN

                  T    = B( L )
                  B( L ) = B( K )
                  B( K ) = T

               END IF

   40       CONTINUE

         END IF

      END IF

      END

      SUBROUTINE SGECO( A, LDA, N, IPVT, RCOND, Z )

c         Factors a real matrix by Gaussian elimination
c         and estimates the condition of the matrix.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c         If  RCOND  is not needed, SGEFA is slightly faster.
c         To solve  A*X = B , follow SGECO by SGESL.

c     On entry

c        A       REAL(LDA, N)
c                the matrix to be factored.

c        LDA     INTEGER
c                the leading dimension of the array  A .

c        N       INTEGER
c                the order of the matrix  A .

c     On return

c        A       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                The factorization can be written  A = L*U , where
c                L  is a product of permutation and unit lower
c                triangular matrices and  U  is upper triangular.

c        IPVT    INTEGER(N)
c                an integer vector of pivot indices.

c        RCOND   REAL
c                an estimate of the reciprocal condition of  A .
c                For the system  A*X = B , relative perturbations
c                in  A  and  B  of size  epsilon  may cause
c                relative perturbations in  X  of size  epsilon/RCOND .
c                If  RCOND  is so small that the logical expression
c                           1.0 + RCOND .EQ. 1.0
c                is true, then  A  may be singular to working
c                precision.  In particular,  RCOND  is zero  if
c                exact singularity is detected or the estimate
c                underflows.

C        Z       REAL(N)
c                a work vector whose contents are usually unimportant.
c                If  A  is close to a singular matrix, then  Z  is
c                an approximate null vector in the sense that
c                norm(A*Z) = RCOND*norm(A)*norm(Z) .

c ------------------------------------------------------------------

c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   LDA, N
      REAL(KR)      RCOND
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL(KR)      A( LDA, * ), Z( * )
c     ..
c     .. Local Scalars ..

      INTEGER   INFO, J, K, KB, KP1, L
      REAL(KR)      ANORM, EK, S, SM, T, WK, WKM, YNORM
c     ..
c     .. External Functions ..

      REAL(KR)      SASUM, SDOT
      EXTERNAL  SASUM, SDOT
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY, SGEFA, SSCAL
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, MAX, SIGN
c     ..


c                        ** compute 1-norm of A
      ANORM  = 0.0E0
      DO 10 J = 1, N
         ANORM  = MAX( ANORM, SASUM( N,A( 1,J ),1 ) )
   10 CONTINUE
c                                      ** factor

      CALL SGEFA( A, LDA, N, IPVT, INFO )

c     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
c     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E .
c     trans(A) is the transpose of A.  The components of E  are
c     chosen to cause maximum local growth in the elements of W  where
c     trans(U)*W = E.  The vectors are frequently rescaled to avoid
c     overflow.

c                        ** solve trans(U)*W = E
      EK = 1.0E0

      DO 20 J = 1, N
         Z( J ) = 0.0E0
   20 CONTINUE


      DO 50 K = 1, N

         IF( Z( K ).NE.0.0E0 ) EK = SIGN( EK, -Z( K ) )

         IF( ABS( EK - Z( K ) ).GT.ABS( A( K,K ) ) ) THEN

            S  = ABS( A( K,K ) ) / ABS( EK - Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            EK = S*EK

         END IF

         WK   = EK - Z( K )
         WKM  = -EK - Z( K )
         S    = ABS( WK )
         SM   = ABS( WKM )

         IF( A( K,K ).NE.0.0E0 ) THEN

            WK   = WK / A( K, K )
            WKM  = WKM / A( K, K )

         ELSE

            WK   = 1.0E0
            WKM  = 1.0E0

         END IF

         KP1  = K + 1

         IF( KP1.LE.N ) THEN

            DO 30 J = KP1, N
               SM     = SM + ABS( Z( J ) + WKM*A( K,J ) )
               Z( J ) = Z( J ) + WK*A( K, J )
               S      = S + ABS( Z( J ) )
   30       CONTINUE

            IF( S.LT.SM ) THEN

               T  = WKM - WK
               WK = WKM

               DO 40 J = KP1, N
                  Z( J ) = Z( J ) + T*A( K, J )
   40          CONTINUE

            END IF

         END IF

         Z( K ) = WK

   50 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )
c                                ** solve trans(L)*Y = W
      DO 60 KB = 1, N
         K  = N + 1 - KB

         IF( K.LT.N )
     &       Z( K ) = Z( K ) + SDOT( N - K, A( K+1, K ), 1, Z( K+1 ), 1)

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

         END IF

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
   60 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )
c                                 ** solve L*V = Y
      YNORM  = 1.0E0

      DO 70 K = 1, N
         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T

         IF( K.LT.N ) CALL SAXPY( N - K, T, A( K + 1,K ), 1, Z( K + 1 ),
     &                            1 )

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM
         END IF

   70 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )
c                                  ** solve  U*Z = V
      YNORM  = S*YNORM

      DO 80 KB = 1, N

         K  = N + 1 - KB

         IF( ABS( Z( K ) ).GT.ABS( A( K,K ) ) ) THEN

            S  = ABS( A( K,K ) ) / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

         IF( A( K,K ).NE.0.0E0 ) Z( K ) = Z( K ) / A( K, K )

         IF( A( K,K ).EQ.0.0E0 ) Z( K ) = 1.0E0

         T  = -Z( K )

         CALL SAXPY( K - 1, T, A( 1,K ), 1, Z( 1 ), 1 )

   80 CONTINUE
c                                   ** make znorm = 1.0
      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = S*YNORM

      IF( ANORM.NE.0.0E0 ) RCOND = YNORM / ANORM
      IF( ANORM.EQ.0.0E0 ) RCOND = 0.0E0

      END

      SUBROUTINE SGEFA( A, LDA, N, IPVT, INFO )

c         Factors a real matrix by Gaussian elimination.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     SGEFA is usually called by SGECO, but it can be called
c     directly with a saving in time if  RCOND  is not needed.
c     (time for SGECO) = (1 + 9/N) * (time for SGEFA) .

c     Input:  same as SGECO

c     On return:

c        A,IPVT  same as SGECO

c        INFO    INTEGER
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  This is not an error
c                     condition for this subroutine, but it does
c                     indicate that SGESL or SGEDI will divide by zero
c                     if called.  Use  RCOND  in SGECO for a reliable
c                     indication of singularity.

c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   INFO, LDA, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL(KR)      A( LDA, * )
c     ..
c     .. Local Scalars ..

      INTEGER   J, K, KP1, L, NM1
      REAL(KR)      T
c     ..
c     .. External Functions ..

      INTEGER   ISAMAX
      EXTERNAL  ISAMAX
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY, SSCAL
c     ..


c                      ** Gaussian elimination with partial pivoting
      INFO = 0
      NM1  = N - 1

      DO 20 K = 1, NM1

         KP1  = K + 1
c                                            ** find L = pivot index

         L  = ISAMAX( N - K + 1, A( K,K ), 1 ) + K - 1
         IPVT( K ) = L

         IF( A( L,K ).EQ.0.0E0 ) THEN
c                                     ** zero pivot implies this column
c                                     ** already triangularized
            INFO = K

         ELSE
c                                     ** interchange if necessary
            IF( L.NE.K ) THEN

               T         = A( L, K )
               A( L, K ) = A( K, K )
               A( K, K ) = T

            END IF
c                                     ** compute multipliers
            T  = -1.0E0 / A( K, K )

            CALL SSCAL( N - K, T, A( K + 1,K ), 1 )

c                              ** row elimination with column indexing
            DO 10 J = KP1, N

               T  = A( L, J )

               IF( L.NE.K ) THEN

                  A( L, J ) = A( K, J )
                  A( K, J ) = T

               END IF

               CALL SAXPY( N-K, T, A( K+1, K ), 1, A( K+1, J ), 1 )

   10       CONTINUE

         END IF

   20 CONTINUE


      IPVT( N ) = N
      IF( A( N,N ) .EQ. 0.0E0 ) INFO = N

      END

      SUBROUTINE SGESL( A, LDA, N, IPVT, B, JOB )

c         Solves the real system
c            A * X = B  or  transpose(A) * X = B
c         using the factors computed by SGECO or SGEFA.

c         Revision date:  8/1/82
c         Author:  Moler, C. B. (U. of New Mexico)

c     On entry

c        A       REAL(LDA, N)
c                the output from SGECO or SGEFA.

c        LDA     INTEGER
c                the leading dimension of the array  A

c        N       INTEGER
c                the order of the matrix  A

c        IPVT    INTEGER(N)
c                the pivot vector from SGECO or SGEFA.

c        B       REAL(N)
c                the right hand side vector.

c        JOB     INTEGER
c                = 0         to solve  A*X = B ,
c                = nonzero   to solve  transpose(A)*X = B

c     On return

c        B       the solution vector  X

c     Error condition

c        A division by zero will occur if the input factor contains a
c        zero on the diagonal.  Technically, this indicates singularity,
c        but it is often caused by improper arguments or improper
c        setting of LDA.  It will not occur if the subroutines are
c        called correctly and if SGECO has set RCOND .GT. 0.0
c        or SGEFA has set INFO .EQ. 0 .

c     To compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue

c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   JOB, LDA, N
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL(KR)      A( LDA, * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER   K, KB, L, NM1
      REAL(KR)      T
c     ..
c     .. External Functions ..

      REAL(KR)      SDOT
      EXTERNAL  SDOT
c     ..
c     .. External Subroutines ..

      EXTERNAL  SAXPY
c     ..


      NM1  = N - 1

      IF( JOB.EQ.0 ) THEN
c                                 ** solve  A * X = B

c                                     ** first solve  L*Y = B
         DO 10 K = 1, NM1

            L  = IPVT( K )
            T  = B( L )

            IF( L.NE.K ) THEN

               B( L ) = B( K )
               B( K ) = T

            END IF

            CALL SAXPY( N - K, T, A( K+1, K ), 1, B( K+1 ), 1 )

   10    CONTINUE
c                                    ** now solve  U*X = Y
         DO 20 KB = 1, N

            K      = N + 1 - KB
            B( K ) = B( K ) / A( K, K )
            T      = - B( K )

            CALL SAXPY( K-1, T, A( 1, K ), 1, B(1), 1 )

   20    CONTINUE


      ELSE
c                         ** solve  trans(A) * X = B

c                                    ** first solve  trans(U)*Y = B
         DO 30 K = 1, N

            T      = SDOT( K - 1, A( 1,K ), 1, B( 1 ), 1 )
            B( K ) = ( B( K ) - T ) / A( K, K )

   30    CONTINUE

c                                    ** now solve  trans(l)*x = y
         DO 40 KB = 1, NM1

            K      = N - KB
            B( K ) = B( K ) + SDOT( N - K, A( K+1, K ), 1, B( K+1 ), 1)
            L      = IPVT( K )

            IF( L.NE.K ) THEN

               T      = B( L )
               B( L ) = B( K )
               B( K ) = T

            END IF

   40    CONTINUE

      END IF

      END

      FUNCTION SASUM( N, SX, INCX )

c  INPUT--    N  Number of elements in vector to be summed
c            SX  Sing-prec array, length 1+(N-1)*INCX, containing vector
c          INCX  Spacing of vector elements in SX

c  OUTPUT-- SASUM   Sum from 0 to N-1 of  ABS(SX(1+I*INCX))
c ----------------------------------------------------------

c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   INCX, N
c     ..
c     .. Array Arguments ..

      REAL(KR)      SX( * ), sasum
c     ..
c     .. Local Scalars ..

      INTEGER   I, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, MOD
c     ..

      SASUM  = 0.0

      IF( N.LE.0 ) RETURN

      IF( INCX.NE.1 ) THEN
c                                          ** non-unit increments
         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SASUM  = SASUM + ABS( SX( I ) )
   10    CONTINUE

      ELSE
c                                          ** unit increments
         M  = MOD( N, 6 )

         IF( M.NE.0 ) THEN
c                             ** clean-up loop so remaining vector
c                             ** length is a multiple of 6.
            DO 20 I = 1, M
               SASUM  = SASUM + ABS( SX( I ) )
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 6
            SASUM  = SASUM + ABS( SX( I ) ) + ABS( SX( I + 1 ) ) +
     &               ABS( SX( I + 2 ) ) + ABS( SX( I + 3 ) ) +
     &               ABS( SX( I + 4 ) ) + ABS( SX( I + 5 ) )
   30    CONTINUE

      END IF

      END

      SUBROUTINE SAXPY( N, SA, SX, INCX, SY, INCY )

c          Y = A*X + Y  (X, Y = vectors, A = scalar)

c  INPUT--
c        N  Number of elements in input vectors X and Y
c       SA  Single precision scalar multiplier A
c       SX  Sing-prec array containing vector X
c     INCX  Spacing of elements of vector X in SX
c       SY  Sing-prec array containing vector Y
c     INCY  Spacing of elements of vector Y in SY

c OUTPUT--
c       SY   For I = 0 to N-1, overwrite  SY(LY+I*INCY) with
c                 SA*SX(LX+I*INCX) + SY(LY+I*INCY),
c            where LX = 1          if INCX .GE. 0,
c                     = (-INCX)*N  if INCX .LT. 0
c            and LY is defined analogously using INCY.
c ------------------------------------------------------------

c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   INCX, INCY, N
      REAL(KR)      SA
c     ..
c     .. Array Arguments ..

      REAL(KR)      SX( * ), SY( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX, IY, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      IF( N.LE.0 .OR. SA.EQ.0.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SY( I ) = SY( I ) + SA*SX( I )
   10    CONTINUE

      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

c                                        ** equal, unit increments
         M  = MOD( N, 4 )

         IF( M.NE.0 ) THEN
c                            ** clean-up loop so remaining vector length
c                            ** is a multiple of 4.
            DO 20 I = 1, M
               SY( I ) = SY( I ) + SA*SX( I )
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 4
            SY( I ) = SY( I ) + SA*SX( I )
            SY( I + 1 ) = SY( I + 1 ) + SA*SX( I + 1 )
            SY( I + 2 ) = SY( I + 2 ) + SA*SX( I + 2 )
            SY( I + 3 ) = SY( I + 3 ) + SA*SX( I + 3 )
   30    CONTINUE


      ELSE
c               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1
         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            SY( IY ) = SY( IY ) + SA*SX( IX )
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE

      END IF

      END

      FUNCTION SDOT( N, SX, INCX, SY, INCY )

c        Single-prec dot product of vectors  X  and  Y

c  INPUT--
c        N  Number of elements in input vectors X and Y
c       SX  Sing-prec array containing vector X
c     INCX  Spacing of elements of vector X in SX
c       SY  Sing-prec array containing vector Y
c     INCY  Spacing of elements of vector Y in SY

c OUTPUT--
c     SDOT   Sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
c            where  LX = 1          if INCX .GE. 0,
c                      = (-INCX)*N  if INCX .LT. 0,
c            and LY is defined analogously using INCY.
c ------------------------------------------------------------------

c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   INCX, INCY, N
c     ..
c     .. Array Arguments ..

      REAL(KR)      SX( * ), SY( * ), sdot
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX, IY, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      SDOT = 0.0

      IF( N.LE.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SDOT = SDOT + SX( I )*SY( I )
   10    CONTINUE


      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

c                                        ** equal, unit increments
         M  = MOD( N, 5 )

         IF( M.NE.0 ) THEN
c                            ** clean-up loop so remaining vector length
c                            ** is a multiple of 4.
            DO 20 I = 1, M
               SDOT = SDOT + SX( I )*SY( I )
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 5
            SDOT = SDOT + SX( I )*SY( I ) + SX( I + 1 )*SY( I + 1 ) +
     &               SX( I + 2 )*SY( I + 2 ) + SX( I + 3 )*SY( I + 3 ) +
     &               SX( I + 4 )*SY( I + 4 )
   30    CONTINUE

      ELSE
c               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1

         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            SDOT = SDOT + SX( IX )*SY( IY )
            IX   = IX + INCX
            IY   = IY + INCY
   40    CONTINUE

      END IF

      END

      SUBROUTINE SSCAL( N, SA, SX, INCX )

c         Multiply vector SX by scalar SA

c  INPUT--  N  Number of elements in vector
c          SA  Single precision scale factor
c          SX  Sing-prec array, length 1+(N-1)*INCX, containing vector
c        INCX  Spacing of vector elements in SX

c OUTPUT-- SX  Replace  SX(1+I*INCX)  with  SA * SX(1+I*INCX)
c                for I = 0 to N-1
c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   INCX, N
      REAL(KR)      SA
c     ..
c     .. Array Arguments ..

      REAL(KR)      SX( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, M
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      IF( N.LE.0 ) RETURN

      IF( INCX.NE.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SX( I ) = SA*SX( I )
   10    CONTINUE


      ELSE

         M  = MOD( N, 5 )

         IF( M.NE.0 ) THEN
c                           ** clean-up loop so remaining vector length
c                           ** is a multiple of 5.
            DO 20 I = 1, M
               SX( I ) = SA*SX( I )
   20       CONTINUE

         END IF
c                             ** unroll loop for speed
         DO 30 I = M + 1, N, 5
            SX( I ) = SA*SX( I )
            SX( I + 1 ) = SA*SX( I + 1 )
            SX( I + 2 ) = SA*SX( I + 2 )
            SX( I + 3 ) = SA*SX( I + 3 )
            SX( I + 4 ) = SA*SX( I + 4 )
   30    CONTINUE

      END IF

      END

      SUBROUTINE SSWAP( N, SX, INCX, SY, INCY )

c          Interchange s.p vectors  X  and  Y, as follows:

c     For I = 0 to N-1, interchange  SX(LX+I*INCX) and SY(LY+I*INCY),
c     where LX = 1          if INCX .GE. 0,
c              = (-INCX)*N  if INCX .LT. 0
c     and LY is defined analogously using INCY.


c  INPUT--
c        N  Number of elements in input vectors X and Y
c       SX  Sing-prec array containing vector X
c     INCX  Spacing of elements of vector X in SX
c       SY  Sing-prec array containing vector Y
c     INCY  Spacing of elements of vector Y in SY

c OUTPUT--
c       SX  Input vector SY (unchanged if N .LE. 0)
c       SY  Input vector SX (unchanged IF N .LE. 0)
c --------------------------------------------------------------

c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   INCX, INCY, N
c     ..
c     .. Array Arguments ..

      REAL(KR)      SX( * ), SY( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, IX, IY, M
      REAL(KR)      STEMP1, STEMP2, STEMP3
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..


      IF( N.LE.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N-1 )*INCX, INCX
            STEMP1 = SX( I )
            SX( I ) = SY( I )
            SY( I ) = STEMP1
   10    CONTINUE


      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

c                                        ** equal, unit increments
         M  = MOD( N, 3 )

         IF( M.NE.0 ) THEN
c                            ** clean-up loop so remaining vector length
c                            ** is a multiple of 3.
            DO 20 I = 1, M
               STEMP1 = SX( I )
               SX( I ) = SY( I )
               SY( I ) = STEMP1
   20       CONTINUE

         END IF
c                              ** unroll loop for speed
         DO 30 I = M + 1, N, 3
            STEMP1 = SX( I )
            STEMP2 = SX( I + 1 )
            STEMP3 = SX( I + 2 )
            SX( I ) = SY( I )
            SX( I + 1 ) = SY( I + 1 )
            SX( I + 2 ) = SY( I + 2 )
            SY( I ) = STEMP1
            SY( I + 1 ) = STEMP2
            SY( I + 2 ) = STEMP3
   30    CONTINUE


      ELSE
c               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1

         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            STEMP1 = SX( IX )
            SX( IX ) = SY( IY )
            SY( IY ) = STEMP1
            IX   = IX + INCX
            IY   = IY + INCY
   40    CONTINUE

      END IF

      END

      INTEGER FUNCTION ISAMAX( N, SX, INCX )

c INPUT--  N     Number of elements in vector of interest
c          SX    Sing-prec array, length 1+(N-1)*INCX, containing vector
c          INCX  Spacing of vector elements in SX

c OUTPUT-- ISAMAX   First I, I = 1 to N, to maximize
c                         ABS(SX(1+(I-1)*INCX))
c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      use params, only: kr
      implicit none
      INTEGER   INCX, N
c     ..
c     .. Array Arguments ..

      REAL(KR)      SX( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I, II
      REAL(KR)      SMAX, XMAG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS
c     ..
c     2017 new initialization
      ISAMAX = 0

      IF( N.LE.0 ) THEN

         ISAMAX = 0

      ELSE IF( N.EQ.1 ) THEN

         ISAMAX = 1

      ELSE

         SMAX = 0.0
         II   = 1

         DO 10 I = 1, 1 + ( N-1 )*INCX, INCX

            XMAG = ABS( SX( I ) )

            IF( SMAX.LT.XMAG ) THEN

               SMAX   = XMAG
               ISAMAX = II

            END IF

            II = II + 1

   10    CONTINUE

      END IF

      END

C================================================================================       
C        CALCULATE PHASE FUNCTION LEGENDRE EXPANSION COEFFICIENTS
C        IN VARIOUS SPECIAL CASES
C--------------------------------------------------------------------------------
C       INPUT: IPHAS   PHASE FUNCTION OPTIONS
C                      1 : ISOTROPIC
C                      2 : RAYLEIGH
C                      3 : HENYEY-GREENSTEIN WITH ASYMMETRY FACTOR -GG-
C                      4 : HAZE L AS SPECIFIED BY GARCIA/SIEWERT
C                      5 : CLOUD C.1 AS SPECIFIED BY GARCIA/SIEWERT
C              GG      ASYMMETRY FACTOR FOR HENYEY-GREENSTEIN CASE
C              NMOM    INDEX OF HIGHEST LEGENDRE COEFFICIENT NEEDED
C                        ( = NUMBER OF STREAMS 'NSTR'  CHOSEN
C                         FOR THE DISCRETE ORDINATE METHOD)
C--------------------------------------------------------------------------------
C      OUTPUT: PMOM(K)  LEGENDRE EXPANSION COEFFICIENTS (K=0 TO NMOM)
C                         (BE SURE TO DIMENSION '0:maxval' IN CALLING PROGRAM)
C--------------------------------------------------------------------------------
C      REFERENCE:  GARCIA, R. AND C. SIEWERT, 1985: BENCHMARK RESULTS
C                     IN RADIATIVE TRANSFER, TRANSP. THEORY AND STAT.
C                     PHYSICS 14, 437-484, TABLES 10 AND 17
C--------------------------------------------------------------------------------
      SUBROUTINE  GETMOM( IPHAS, GG, NMOM, PMOM )
      
      use params, only: kr
      implicit none

      INTEGER  IPHAS, NMOM, k

      REAL(KR)  GG, PMOM(0:*), HAZELM(82), CLDMOM(299)
      
      DATA HAZELM /  2.41260, 3.23047, 3.37296, 3.23150, 2.89350,
     A               2.49594, 2.11361, 1.74812, 1.44692, 1.17714,
     B               0.96643, 0.78237, 0.64114, 0.51966, 0.42563,
     C               0.34688, 0.28351, 0.23317, 0.18963, 0.15788,
     D               0.12739, 0.10762, 0.08597, 0.07381, 0.05828,
     E               0.05089, 0.03971, 0.03524, 0.02720, 0.02451,
     F               0.01874, 0.01711, 0.01298, 0.01198, 0.00904,
     G               0.00841, 0.00634, 0.00592, 0.00446, 0.00418,
     H               0.00316, 0.00296, 0.00225, 0.00210, 0.00160,
     I               0.00150, 0.00115, 0.00107, 0.00082, 0.00077,
     J               0.00059, 0.00055, 0.00043, 0.00040, 0.00031,
     K               0.00029, 0.00023, 0.00021, 0.00017, 0.00015,
     L               0.00012, 0.00011, 0.00009, 0.00008, 0.00006,
     M               0.00006, 0.00005, 0.00004, 0.00004, 0.00003,
     N               0.00003, 3*0.00002, 8*0.00001 /

      DATA  ( CLDMOM(K), K = 1, 159 ) /
     A  2.544,  3.883,  4.568,  5.235,  5.887,  6.457,  7.177,  7.859,
     B  8.494,  9.286,  9.856, 10.615, 11.229, 11.851, 12.503, 13.058,
     C 13.626, 14.209, 14.660, 15.231, 15.641, 16.126, 16.539, 16.934,
     D 17.325, 17.673, 17.999, 18.329, 18.588, 18.885, 19.103, 19.345,
     E 19.537, 19.721, 19.884, 20.024, 20.145, 20.251, 20.330, 20.401,
     F 20.444, 20.477, 20.489, 20.483, 20.467, 20.427, 20.382, 20.310,
     G 20.236, 20.136, 20.036, 19.909, 19.785, 19.632, 19.486, 19.311,
     H 19.145, 18.949, 18.764, 18.551, 18.348, 18.119, 17.901, 17.659,
     I 17.428, 17.174, 16.931, 16.668, 16.415, 16.144, 15.883, 15.606,
     J 15.338, 15.058, 14.784, 14.501, 14.225, 13.941, 13.662, 13.378,
     K 13.098, 12.816, 12.536, 12.257, 11.978, 11.703, 11.427, 11.156,
     L 10.884, 10.618, 10.350, 10.090,  9.827,  9.574,  9.318,  9.072,
     M  8.822, 8.584, 8.340, 8.110, 7.874, 7.652, 7.424, 7.211, 6.990,
     N  6.785, 6.573, 6.377, 6.173, 5.986, 5.790, 5.612, 5.424, 5.255,
     O  5.075, 4.915, 4.744, 4.592, 4.429, 4.285, 4.130, 3.994, 3.847,
     P  3.719, 3.580, 3.459, 3.327, 3.214, 3.090, 2.983, 2.866, 2.766,
     Q  2.656, 2.562, 2.459, 2.372, 2.274, 2.193, 2.102, 2.025, 1.940,
     R  1.869, 1.790, 1.723, 1.649, 1.588, 1.518, 1.461, 1.397, 1.344,
     S  1.284, 1.235, 1.179, 1.134, 1.082, 1.040, 0.992, 0.954, 0.909 /
      DATA  ( CLDMOM(K), K = 160, 299 ) /
     T  0.873, 0.832, 0.799, 0.762, 0.731, 0.696, 0.668, 0.636, 0.610,
     U  0.581, 0.557, 0.530, 0.508, 0.483, 0.463, 0.440, 0.422, 0.401,
     V  0.384, 0.364, 0.349, 0.331, 0.317, 0.301, 0.288, 0.273, 0.262,
     W  0.248, 0.238, 0.225, 0.215, 0.204, 0.195, 0.185, 0.177, 0.167,
     X  0.160, 0.151, 0.145, 0.137, 0.131, 0.124, 0.118, 0.112, 0.107,
     Y  0.101, 0.097, 0.091, 0.087, 0.082, 0.079, 0.074, 0.071, 0.067,
     Z  0.064, 0.060, 0.057, 0.054, 0.052, 0.049, 0.047, 0.044, 0.042,
     A  0.039, 0.038, 0.035, 0.034, 0.032, 0.030, 0.029, 0.027, 0.026,
     B  0.024, 0.023, 0.022, 0.021, 0.020, 0.018, 0.018, 0.017, 0.016,
     C  0.015, 0.014, 0.013, 0.013, 0.012, 0.011, 0.011, 0.010, 0.009,
     D  0.009, 3*0.008, 2*0.007, 3*0.006, 4*0.005, 4*0.004, 6*0.003,
     E  9*0.002, 18*0.001 /

      IF ( IPHAS.LT.1 .OR. IPHAS.GT.5 )
     $     call errmsg(13,'GETMOM--BAD INPUT VARIABLE IPHAS')
      IF ( IPHAS.EQ.3 .AND. (GG.LE.-1.0 .OR. GG.GE.1.0) )
     $     call errmsg(14,'GETMOM--BAD INPUT VARIABLE GG')
      IF ( NMOM.LT.2 )
     $     call errmsg(15,'GETMOM--BAD INPUT VARIABLE NMOM')
      
      PMOM(0) = 1.0
      DO  10  K = 1, NMOM
        PMOM(K) = 0.0
 10   CONTINUE
      
      IF ( IPHAS.EQ.2 )  THEN
C---------------------------
C    RAYLEIGH PHASE FUNCTION
C---------------------------

        PMOM(2) = 0.1
        
      ELSE IF ( IPHAS.EQ.3 ) THEN
C------------------------------
C   HENYEY-GREENSTEIN PHASE FCN
C------------------------------
        DO  20  K = 1, NMOM
          PMOM(K) = GG**K
 20     CONTINUE
        
      ELSE IF ( IPHAS.EQ.4 ) THEN
C---------------------------
C    HAZE-L PHASE FUNCTION
C---------------------------
        DO  30  K = 1, MIN0(82,NMOM)
          PMOM(K) = HAZELM(K) / ( 2*K+1 )
 30     CONTINUE
        
      ELSE IF ( IPHAS.EQ.5 ) THEN
C---------------------------
C   CLOUD C.1 PHASE FUNCTION
C---------------------------
        DO  40  K = 1, MIN0(298,NMOM)
          PMOM(K) = CLDMOM(K) / ( 2*K+1 )
 40     CONTINUE
        
      END IF
      
      RETURN
      END
