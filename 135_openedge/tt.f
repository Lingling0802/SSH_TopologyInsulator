	program test
	Implicit none

	Integer l, ld
	Integer maxl, maxld
        Integer maxcolors
        parameter( maxcolors = 64  )  ! I know this many colors
	parameter( maxl = 500   )     ! must match definition of
			  	      ! maxsize in WriteLattice
	parameter( maxld = maxl * maxl )

	Real site( maxld) 
	Integer isite(1:maxld) 
	Integer output(maxld), biggest, i
	Integer error
	Real z, maxZ, minZ, norm

	Integer done
      	Integer numc, end, found
	Integer maxstring
	parameter( maxstring = 40 )
	Character*40 filename, outfile

	done = 0
	DO while ( done .ne. 1 )

	  Read(5,'(a)') filename
	  IF( filename .ne. 'end' ) THEN 
	    error = 0
	    Open(66,file=filename,status='old',iostat=error,err=101 )
	    i = 0
	    DO while( ( error .eq. 0 ) .and. ( i .lt. maxld ))
	      i = i + 1
	      Read(66,*,end=101) site(i)
	      IF( i .eq. 1 ) THEN
		maxZ = site(i)
		minZ = site(i)
	      ELSE
	        maxZ = max( site(i), maxZ ) 
	        minZ = min( site(i), minZ ) 
	      ENDIF
	    ENDDO
 101        CONTINUE
            IF( i .eq. maxld) THEN
                Write(6,*)'System is too big'
                Write(6,*)'Increase maxl and recompile'
                stop
            ENDIF

	    ld = i -1
	    l = sqrt(float(ld))
	    Write(6,*)'Number of points is ', ld
	    IF ( maxZ .eq. minZ ) THEN
		Write(6,*)'Max and Min are the same. Making adjustment'
		IF( maxZ .gt. 0 ) THEN
		   maxZ = 1.5*maxZ
		   minZ = 0.5*maxZ
		ELSE
		   maxZ = 0.5*maxZ
		   minZ = 1.5*maxZ
		ENDIF
	    ENDIF
	    close(66)
	    found = 0
! Let's bin the numbers in to maxcolours bins.
	    norm =1.0e0/( maxZ - minZ )
	    DO i = 1, ld
		isite(i) = int((maxcolors-1)*( site(i) - minZ ) * norm)
	    ENDDO
	   

!  Here we create the output file name
	    end=1
	    DO i = 1, maxstring-4
	       IF ( found .ne. 1 )  THEN
	         IF ( filename(i:i) .ne. " " ) THEN
		   end = end + 1
	         ELSE
		   found = 1
	         ENDIF
	       ENDIF
	    ENDDO
	    outfile = filename(1:end-1)//".xpm"


	    Call WriteLattice( isite(1:ld),l,ld,0,2147483647,outfile )	

	    ELSE
	       done = 1
	    ENDIF
	ENDDO
	stop
	end

 	Subroutine WriteLattice(latin,l,ld,number,background,filename) 
	IMPLICIT none 

	Integer l, ld        ! denote the size of the lattice
	Integer latin( ld )  ! The lattice itself
	Integer background   ! label of back ground sites.
	Integer number       ! the largest number of a cluster.
	Character*40 filename



	Integer i,j
	Integer maxcolors, maxsize
	parameter( maxsize = 500 )    ! I can deal with l this big
	parameter( maxcolors = 64  )  ! I know this many colors
	Real temp
	character color(0:maxcolors-1)
	character*(maxsize) spin
	character*7 colormap(0:maxcolors-1)
	Integer numc, pick, site

        data color / "a","b","c","d","e","f","g","h","i","j",
     &               "k","l","m","n","o","p","q","r","s","t",
     &               "u","v","w","x","y","z","1","2","3","4",
     &               "5","6","7","8","9","0","!","@","#","$",
     &               "A","B","C","D","E","F","G","H","I","J",
     &               "K","L","M","N","O","P","Q","R","S","T",
     &               "U","V","W","X"/



	data colormap /
     &  "#00008F","#00009F","#0000AF","#0000BF","#0000CF","#0000DF",
     &  "#0000EF","#0000FF","#000FFF","#001FFF","#002FFF","#003FFF",
     &  "#004FFF","#005FFF","#006FFF","#007FFF","#008FFF","#009FFF",
     &  "#00AFFF","#00BFFF","#00CFFF","#00DFFF","#00EFFF","#00FFFF",
     &  "#0FFFFF","#1FFFEF","#2FFFDF","#3FFFCF","#4FFFBF","#5FFFAF",
     &  "#6FFF9F","#7FFF8F","#8FFF7F","#9FFF6F","#AFFF5F","#BFFF4F",
     &  "#CFFF3F","#DFFF2F","#EFFF1F","#FFFF0F","#FFFF00","#FFEF00",
     &  "#FFDF00","#FFCF00","#FFBF00","#FFAF00","#FF9F00","#FF8F00",
     &  "#FF7F00","#FF6F00","#FF5F00","#FF4F00","#FF3F00","#FF2F00",
     &  "#FF1F00","#FF0F00","#FF0000","#EF0000","#DF0000","#CF0000",
     &  "#BF0000","#AF0000","#9F0000","#8F0000"/
	IF( l .gt. maxsize ) THEN
		Write(6,*) "Lattice is too big. Recompile"
		return
	ENDIF

	Open( 66, file = filename, status = 'unknown' )

c	numc = MIN( number, maxcolors )  ! only use the colors needed
	numc =  maxcolors 

c The header stuff for the Xpm file

	Write(66,10)
	Write(66,20)
	Write(66,30)
	Write(66,40) l, l, numc, 1 
	Write(66,50)

c Save the needed colors

	DO i = 0, numc-1
		Write(66,60) color(i), colormap(i) 
	ENDDO
	Write(66,80)

c Write out the data.  
	Do i = 1, l
	  DO j = 1, l
		site = (i-1)*l+j
		pick= latin(site)-numc*((latin(site)-1)/numc)
		IF( latin(site) .eq. background ) pick = 0
		spin(j:j) = color( pick )
	  ENDDO
	      IF( i .ne. l ) THEN
	         Write(66,100)spin(1:l)
 100             format('"',(a),'",' )
	      ELSE
	         Write(66,110)spin(1:l)
 110             format('"',(a),'"' )
	      ENDIF
	ENDDO

	Write(66,90)

	Close( 66 )

	
 10    format( '/* XPM */' )
 20    format( 'static char *test[] = {' )
 30    format( '/* width height num_colors chars_per_pixel */ ' )
 40    format( '"',2x, i3, 2x, i3, 2x , i3, 2x , i3, '",' )
 50    format( '/* colors */ ' )
 60    format( '"',a1, '  c  ', a7,'",' )
 80    format( '/* pixels */'  )
 90    format( '};' )

	return
	end

