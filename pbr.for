!Travel time database
      Module ray
      parameter(maxnlat=100)
      parameter(maxnlon=100)
      parameter(maxndep=50)
! 20150130_Let the matrix be allocatalbe to safe the memory.
      real*8, ALLOCATABLE :: lat_c(:), lon_c(:), dep_c(:)
      real*8, ALLOCATABLE :: vel_p(:,:,:),vel_s(:,:,:)
      real*8 bld3,bld4,ro
      integer nlat_c,nlon_c,ndep_c,nxyz_c,nxy_c,nx_c,ips
      parameter(ilatdeg=1000000)
      parameter(ilondeg=1000000)
      parameter(idepkm=1000000)
      real*8 lat1_c,lon1_c,dep1_c
      integer ilonloc_c(ilondeg),ilatloc_c(ilatdeg),ideploc_c(idepkm)
      End Module ray
c      module observ
c      parameter (maxsta=3000)
c      parameter (maxobs=3000)
c      character*4 stn(maxsta)
c      integer ltds(maxsta),lnds(maxsta)
c      real sltm(maxsta),slnm(maxsta),stc(3,maxsta)
c      real pcor(maxsta),scor(maxsta),spcor(maxsta)
c      integer isto(maxobs),ifm(maxobs),inten(maxobs),isp(maxobs)
c      real secp(maxobs,2),epdis(maxobs),azimuth(maxobs),takeoff(maxobs)
c      real wa(maxobs),xpga(maxobs),xweio(maxobs,2),res(maxobs,2)
c      real xmls(maxobs)
c      real wa1(maxobs),xmls1(maxobs)
c      end module observ
      program main
      use ray
      implicit real*8 (a-h,o-z)
      parameter (msg=16384)
      real*8 w(3,msg+1)
      integer np
      real*8 tt
      REAL :: tempo(3), temp
      REAL(8),ALLOCATABLE :: evlo(:), evla(:), evdp(:)
      REAL(8),ALLOCATABLE :: stlo(:), stla(:), stel(:)
      REAL(8),ALLOCATABLE :: afs_stlo(:,:), afs_stla(:,:), afs_stel(:,:)
      REAL(8),ALLOCATABLE :: evt_la(:), sta_la(:)
      REAL(8),ALLOCATABLE :: ptt(:,:)
      INTEGER :: num_sta, num_sou
      character*20 evlo_tmp,evla_tmp,evdp_tmp,stlo_tmp,stla_tmp,stel_tmp
      character*80 ORDER
      real*8 evlo_tmp1
      call input_vel
! Read input data
      num_sta = 0
      num_sou = 0
      evlo_tmp=' '
      evla_tmp=' '
      evdp_tmp=' '
      stlo_tmp=' '
      stla_tmp=' '
      stel_tmp=' '
      nmarg=iargc()
      if (nmarg.eq.6) then
        num_sou=1
        num_sta=1

        ALLOCATE( evlo(num_sou), evla(num_sou), evdp(num_sou) )
        ALLOCATE( stlo(num_sta), stla(num_sta), stel(num_sta) )
        ALLOCATE( afs_stlo(num_sou,num_sta), afs_stla(num_sou,num_sta))
        ALLOCATE( afs_stel(num_sou,num_sta) )
        ALLOCATE( evt_la(num_sou), sta_la(num_sta) )
        ALLOCATE( ptt(num_sou,num_sta) )
        call getarg(1,evlo_tmp)
        read(evlo_tmp,*)evlo(1)
        call getarg(2,evla_tmp)
        read(evla_tmp,*)evla(1)
        call getarg(3,evdp_tmp)
        read(evdp_tmp,*)evdp(1)
        call getarg(4,stlo_tmp)
        read(stlo_tmp,*)stlo(1)
        call getarg(5,stla_tmp)
        read(stla_tmp,*)stla(1)
        call getarg(6,stel_tmp)
        read(stel_tmp,*)stel(1)
        evt_la(1) = geog_to_geoc(evla(1))
        stel(1) = -1.d0*stel(1)/1000.d0
        sta_la(1) = geog_to_geoc(stla(1))
      else
        write(*,*) 'input format: pbr [source_lon] [lat] [dep] [sta_lon]
     &[lat] [elv]'
        stop
      endif
cc      OPEN(11,status="old",file="sou_location.txt")
cc      OPEN(12,status="old",file="sta_location.txt")
c      DO WHILE ( .NOT. eof(11) )
cc      do i=1,1000
cc       READ(11,*,end=13) tempo(:)
cc        num_sou = num_sou + 1
cc      END DO
cc13    continue
cc      REWIND(11)
c      DO WHILE ( .NOT. eof(12) )
cc      do i=1,1000
cc        READ(12,*,end=14) tempo(:)
cc        num_sta = num_sta + 1
cc      END DO
cc14    continue
cc      REWIND(12)
cc      ALLOCATE( evlo(num_sou), evla(num_sou), evdp(num_sou) )
cc      ALLOCATE( stlo(num_sta), stla(num_sta), stel(num_sta) )
cc      ALLOCATE( afs_stlo(num_sou,num_sta), afs_stla(num_sou,num_sta))
cc      ALLOCATE( afs_stel(num_sou,num_sta) )
cc      ALLOCATE( evt_la(num_sou), sta_la(num_sta) )
cc      ALLOCATE( ptt(num_sou,num_sta) )
cc      DO i = 1,num_sou
cc        READ(11,*) evlo(i), evla(i), evdp(i)
cc        evt_la(i) = geog_to_geoc(evla(i))
cc      END DO
cc      DO i = 1,num_sta
cc        READ(12,*) stlo(i), stla(i), stel(i)
cc       stel(i) = -1.d0*stel(i)/1000.d0
cc        sta_la(i) = geog_to_geoc(stla(i))
cc      END DO
cc      write(*,*) 'evt_la, sta_la=',evt_la,sta_la
cc      stop
      ips=1
      DO i = 1,num_sou
        DO j = 1,num_sta
          call pbr(evt_la(i),evlo(i),evdp(i),sta_la(j),stlo(j),stel(j),w
     :,np,tt)
          write(*,'(a,f10.4)') ' P wave travel time: ',tt
          ptt(i,j) = tt
          ORDER=' '
          ORDER='mv out_route out_route_P'
          kerr=system(ORDER)
        END DO
      END DO
      ips=2
cc      write(*,*) 'evt_la,lo,dp,sta_la,lo,el=',evt_la,evlo,evdp,sta_la,st
cc     &lo,stel
cc      stop
      DO i = 1,num_sou
        DO j = 1,num_sta
          call pbr(evt_la(i),evlo(i),evdp(i),sta_la(j),stlo(j),stel(j),w
     &,np,tt)
          write(*,'(a,f10.4)') ' S wave travel time: ',tt
          ORDER=' '
          ORDER='mv out_route out_route_S'
          kerr=system(ORDER)
        END DO
      END DO
! Bubble sorting and output
      DO i = 1,num_sou
        afs_stlo(i,:) = stlo(:)
        afs_stla(i,:) = stla(:)
        afs_stel(i,:) = stel(:)
        DO j = num_sta-1,1,-1
          DO k = 1,j
            IF ( ptt(i,k) > ptt(i,k+1) ) THEN
              temp = ptt(i,k)
              tempo(1) = afs_stlo(i,k)
              tempo(2) = afs_stla(i,k)
              tempo(3) = afs_stel(i,k)
              ptt(i,k) = ptt(i,k+1)
              afs_stlo(i,k) = afs_stlo(i,k+1)
              afs_stla(i,k) = afs_stla(i,k+1)
              afs_stel(i,k) = afs_stel(i,k+1)
              ptt(i,k+1) = temp
              afs_stlo(i,k+1) = tempo(1)
              afs_stla(i,k+1) = tempo(2)
              afs_stel(i,k+1) = tempo(3)
            END IF
          END DO
        END DO
      END DO
      OPEN(21,file="RESULTS.txt")
      DO i = 1,num_sou
        WRITE(21,"(A4,2F10.4,F11.4,A13)") "0", evlo(i), evla(i),evdp(i),
     &"10.00"
        DO j = 1,num_sta
          WRITE(21,"(I4,2F10.4,F11.4,F13.6)") j, afs_stlo(i,j), afs_stla
     &(i,j),afs_stel(i,j)*(-1000.d0), ptt(i,j)
        END DO
      END DO
      end program main
      real*8 function geog_to_geoc(xla)
      implicit none
      real*8 xla,RAD_PER_DEG,B2A_SQ
      RAD_PER_DEG=0.0174532925199432955
      B2A_SQ=0.993305521
      geog_to_geoc = atan(B2A_SQ*tan(RAD_PER_DEG*xla)) / RAD_PER_DEG
      return
      end function geog_to_geoc
      real*8 function geoc_to_geog(xla)
      implicit none
      real*8 xla,RAD_PER_DEG,B2A_SQ
      RAD_PER_DEG=0.0174532925199432955
      B2A_SQ=0.993305521
      geoc_to_geog = atan(tan(RAD_PER_DEG*xla)/B2A_SQ) / RAD_PER_DEG
      return
      end function geoc_to_geog
      subroutine input_vel
      use ray
      implicit real*8 (a-h,o-z)
c      open(1,file="MOD_H13",status='old')
      open(1,file="MOD_CWB1D",status='old')
      read(1,*)bld3,bld4,nlon_c,nlat_c,ndep_c
! 20150130_Allocate the matirx size.
      ALLOCATE(lat_c(nlat_c),lon_c(nlon_c),dep_c(ndep_c),vel_p(nlon_c,nl
     &at_c,ndep_c),vel_s(nlon_c,nlat_c,ndep_c))
      read(1,*)(lon_c(i),i=1,nlon_c)
      read(1,*)(lat_c(i),i=1,nlat_c)
      read(1,*)(dep_c(i),i=1,ndep_c)
!-- read P velocity model
      do k=1,ndep_c
        do j=1,nlat_c
          read(1,*)(vel_p(i,j,k),i=1,nlon_c)
        enddo
      enddo
!-- read S velocity model
      do k=1,ndep_c
        do j=1,nlat_c
          read(1,*)(vel_s(i,j,k),i=1,nlon_c)
        enddo
      enddo
      close(1)
      call bldmap
      nxyz_c=nlon_c*nlat_c*ndep_c
      nxy_c=nlon_c*nlat_c
      nx_c=nlon_c
      nxyz2_c=(nlon_c-2)*(nlat_c-2)*(ndep_c-2)
      nxy2_c=(nlon_c-2)*(nlat_c-2)
      nx2_c=nlon_c-2
      ave=0.0
      do i=1,nlat_c
        ave=ave+lat_c(i)
      enddo
      ave=ave/real(nlat_c)
      ro=earthr(ave)
      end subroutine input_vel
      subroutine bldmap
      use ray
      implicit real*8 (a-h,o-z)
      real*8 lon_now,lat_now,dep_now
!-- for crustal velocity
      lon1_c=bld3-lon_c(1)
      ilonmax=(1e-10)+(lon_c(nlon_c)+lon1_c)/bld3
      lat1_c=bld3-lat_c(1)
c      write(*,*)'nlat_c=',nlat_c,'bld3=',bld3,'lat1_c=',lat1_c
      ilatmax=(1e-10)+(lat_c(nlat_c)+lat1_c)/bld3
      dep1_c=bld4-dep_c(1)
      idepmax=(1e-10)+(dep_c(ndep_c)+dep1_c)/bld4
      if ((ilonmax.gt.ilondeg).or.(ilatmax.gt.ilatdeg).or.(idepmax.gt.id
     &epkm))then
        print*,"Error, model dimension out of range!"
        stop
      endif
      ilon=1
      do i=1,ilonmax
        ilon1=ilon+1
        lon_now=float(i)*bld3-lon1_c
        if (lon_now.ge.lon_c(ilon1)) ilon=ilon1
        ilonloc_c(i)=ilon
      enddo
      do i=ilonmax+1,ilondeg
        ilonloc_c(i)=0
      enddo
      ilat=1
cc      write(*,*) 'ilatmax=',ilatmax
      do i=1,ilatmax
        ilat1=ilat+1
        lat_now=float(i)*bld3-lat1_c
        if (lat_now.ge.lat_c(ilat1)) ilat=ilat1
        ilatloc_c(i)=ilat
      enddo
      do i=ilatmax+1,ilatdeg
        ilatloc_c(i)=0
      enddo
      idep=1
      do i=1,idepmax
        idep1=idep+1
        dep_now=float(i)*bld4-dep1_c
        if (dep_now.ge.dep_c(idep1)) idep=idep1
        ideploc_c(i)=idep
      enddo
      do i=idepmax+1,idepkm
        ideploc_c(i)=0
      enddo
      end subroutine bldmap
      subroutine intmap_3d(lon,lat,dep,ip,jp,kp)
      use ray
      implicit real*8 (a-h,o-z)
      real*8 lon,lat,dep
      ip=int(1e-10+(lon+lon1_c)/bld3)
      jp=int(1e-10+(lat+lat1_c)/bld3)
      kp=int(1e-10+(dep+dep1_c)/bld4)
      if ((ip.le.0).or.(jp.le.0).or.(kp.le.0)) then
        print*,"Error,lon,lat,dep out of range!"
        print*,"lon=",lon,"lat=",lat,"dep=",dep
        print*,"ip,jp,kp",ip,jp,kp
! 20150130_Take "stop" away and skip the calculation.
        GOTO 130
      endif
      ip=ilonloc_c(ip)
      jp=ilatloc_c(jp)
      kp=ideploc_c(kp)
      if ((ip.eq.0).or.(jp.eq.0).or.(kp.eq.0)) then
        print*,"Error,crust lon,lat out of range!"
        print*,"lon=",lon,"lat=",lat,"dep=",dep
        print*,"ip,jp,kp",ip,jp,kp
! 20150130 Take "stop" away.
      endif
 130  return
      end subroutine intmap_3d
      function velocity(r,pa,ra)
      use ray
      implicit real*8 (a-h,o-z)
      real*8 lat,lon,dep,shiftlo
      real*8 r,pa,ra,r2d
      real*8 v,velocity
      common /coord/ shiftlo
      r2d = 90./asin(1.)
      lat=geoc_to_geog(90.0-pa*r2d)
      lon=ra*r2d+shiftlo
      dep=ro-r
cc      write(*,*) 'lon,lat,dep,v@vel3=',lon,lat,dep,v
      call vel3(lon,lat,dep,v)
cc      write(*,*) 'lon,lat,dep,v@@vel3=',lon,lat,dep,v
      velocity=v
      return
      end function velocity
      
      subroutine vel3(lon,lat,dep,v)
      use ray
      implicit real*8 (a-h,o-z)
      real*8 lon,lat,dep,v
      real*8 lonf,lonf1,latf,latf1,depf,depf1
      real*8 wv(2,2,2)
      common /weight/ wv,ip,jp,kp
      call intmap_3d(lon,lat,dep,ip,jp,kp)
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
      if((ip1.gt.nlon_c).or.(jp1.gt.nlat_c).or.(kp1.gt.ndep_c))then
        print*,"Error, ip1,jp1,kp1 out of range!"
        print*,"ip1=",ip1,"jp1=",jp1,"kp1=",kp1
! 20150130_Take "stop" away and skip the calculation.
        GOTO 140
      endif
! 20150130_Add jugdement to avoid "out of range".
      if ((ip.le.0).or.(jp.le.0).or.(kp.le.0)) then
        GOTO 140
      end if
      lonf=(lon-lon_c(ip))/(lon_c(ip1)-lon_c(ip))
      latf=(lat-lat_c(jp))/(lat_c(jp1)-lat_c(jp))
      depf=(dep-dep_c(kp))/(dep_c(kp1)-dep_c(kp))
      lonf1=1.0-lonf
      latf1=1.0-latf
      depf1=1.0-depf
      wv(1,1,1)=lonf1*latf1*depf1
      wv(2,1,1)=lonf*latf1*depf1
      wv(1,2,1)=lonf1*latf*depf1
      wv(2,2,1)=lonf*latf*depf1
      wv(1,1,2)=lonf1*latf1*depf
      wv(2,1,2)=lonf*latf1*depf
      wv(1,2,2)=lonf1*latf*depf
      wv(2,2,2)=lonf*latf*depf
      if(ips.eq.2)then
        v= wv(1,1,1)*vel_s(ip,jp,kp)+wv(2,1,1)*vel_s(ip1,jp,kp)+wv(1,2,1
     &)*vel_s(ip,jp1,kp) +wv(2,2,1)*vel_s(ip1,jp1,kp)+wv(1,1,2)*vel_s(ip
     &,jp,kp1) +wv(2,1,2)*vel_s(ip1,jp,kp1)+wv(1,2,2)*vel_s(ip,jp1,kp1)+
     &wv(2,2,2)*vel_s(ip1,jp1,kp1)
      else
        v= wv(1,1,1)*vel_p(ip,jp,kp)+wv(2,1,1)*vel_p(ip1,jp,kp)+wv(1,2,1
     &)*vel_p(ip,jp1,kp) +wv(2,2,1)*vel_p(ip1,jp1,kp)+wv(1,1,2)*vel_p(ip
     &,jp,kp1) +wv(2,1,2)*vel_p(ip1,jp,kp1)+wv(1,2,2)*vel_p(ip,jp1,kp1)+
     &wv(2,2,2)*vel_p(ip1,jp1,kp1)
      endif
 140  return
      end subroutine vel3
      subroutine pbr(evla,evlo,evdp,stla,stlo,stel,w,np,tk)
      use ray
      implicit real*8(a-h,o-z)
      parameter (msg=16384)
      real*8 w(3,msg+1)
      real*8 r(msg+1), a(msg+1), b(msg+1)
      integer ni,i
      real*8 shiftlo
      real*8 aas,bbs,hs,aar,bbr,hr
      real*8 xfac,flim,mins
      real*8 dpi,r2d
      real*8 velocity,rtim
      real*8 tk
      real*8 as,ar
      real*8 bre,bso,dlo
      real*8 ad,rs,rr
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,dx,dy,dz
      real*8 r1,a1,b1,r2,a2,b2,r3,a3,b3
      real*8 x,y,z,acosa,sina,cosa,to,tp
      real*8 dn,ddn,dr,da,db
      real*8 dseg,ddseg
      real*8 v1,v2,v3
      real*8 upz,dwz
      real*8 vr1,vr2,vr,vb1,vb2,vb,va1,va2,va
      real*8 pr,pa,pb
      real*8 vrd,rvr,rva,rvb,rvs
      real*8 cc,rcur,rdr,rda,rdb,rpr,ap,bp
      real*8 adV,bdV,rdV
      real*8 RNULL
      common /coord/ shiftlo
      data RNULL /0.0e10/
! right now force the receiver at elevation of 0
      aas=evla
      bbs=evlo
      hs=evdp
      aar=stla
      bbr=stlo
      hr=stel
      ni = msg+1
      xfac = 1.5
      n1 = 2
      n2 = msg
      nloop = 12800
      flim = 1.e-4/100.
      mins = 2.
      dpi = asin(1.)/ 90.
      r2d = 90./asin(1.)
!-- Check coordinates
      if(aas.LT.-90.OR.aas.GT.90.)then
        write(*,*)'Latitude of source is out of range'
        stop
      endif
      if(aar.LT.-90.OR.aar.GT.90.)then
        write(*,*)'Latitude of station is out of range'
        stop
      endif
      if(bbs.LT.-180.OR.bbs.GT.180.)then
        write(*,*)'Longitude of source is out of range'
        stop
      endif
      if(bbr.LT.-180.OR.bbr.GT.180.)then
        write(*,*)'Longitude of station is out of range'
        stop
      endif
!-- longitude and latitude range from 0 to 180.
!-- This program does not work with angles
!-- greater than 180.
!-- Pass from latitude to colatitude
      as = (90.00-aas) * dpi
      ar = (90.00-aar) * dpi
      if(bbr.LT.0.0)then
        bre=360.+bbr
      else
        bre=bbr
      endif
      if(bbs.LT.0.0)then
        bso=360.+bbs
      else
        bso=bbs
      endif
      dlo=abs(bso-bre)
      if(dlo.LT.180.)then
        shiftlo=0.0e10
        if(bso.LT.bre)then
          shiftlo=bso-(180.-dlo)/2.
          bbs=(180.-dlo)/2.
          bbr=bbs+dlo
        else
          shiftlo=bre-(180.-dlo)/2.
          bbr=(180.-dlo)/2.
          bbs=bbr+dlo
        endif
      else
        dlo=360.0000-dlo
        shiftlo=0.0e10
        if(bso.LT.bre)then
          shiftlo=bso-(dlo+(180.-dlo)/2.)
          bbs=(180.-dlo)/2.+dlo
          bbr=bbs-dlo
        else
          shiftlo=bre-(dlo+(180.-dlo)/2.)
          bbr=(180.-dlo)/2.+dlo
          bbs=bbr-dlo
        endif
      endif
      bs = bbs * dpi
      br = bbr * dpi
      ad = (as + ar) / 2.
      rs = ro - hs
      rr = ro - hr
! *** initial straight ray ***
! ni : number of ray segments
      ni = n1
      x1 = rs*sin(as)*cos(bs)
      y1 = rs*sin(as)*sin(bs)
      z1 = rs*cos(as)
      x2 = rr*sin(ar)*cos(br)
      y2 = rr*sin(ar)*sin(br)
      z2 = rr*cos(ar)
      dx = x2-x1
      dy = y2-y1
      dz = z2-z1
      dlen=sqrt(dx*dx+dy*dy+dz*dz)
      if (ni.lt.2) ni=2
      dx = (x2-x1) / ni
      dy = (y2-y1) / ni
      dz = (z2-z1) / ni
      do j=1,ni+1
        x = x1 + dx*(j-1)
        y = y1 + dy*(j-1)
        z = z1 + dz*(j-1)
        r(j) = sqrt(x**2 + y**2 + z**2)
        acosa=z/r(j)
        if(acosa.LT.-1.)acosa=-1.
        if(acosa.GT.1)acosa=1.
        a(j) = acos(acosa)
        acosa=x/r(j)/sin(a(j))
        if(acosa.LT.-1.)acosa=-1.
        if(acosa.GT.1)acosa=1.
        b(j) = acos(acosa)
        if(y.LT.0.00000)b(j)=360.00000*dpi-b(j)
      enddo
      to = rtim(ni+1,r,a,b)
      tp = to
      do i=1,ni+1
        w(1,i) = r(i)
        w(2,i) = a(i)
        w(3,i) = b(i)
      enddo
! *** number of points loop ***
      loops = 0
      do while(ni .le. n2)
! *** interation loop ***
        do l=1,nloop
          loops = loops + 1
          do kk=2,ni
!-- see um & thurber (1987) p.974.
            if(mod(kk,2) .eq. 0) then
              k = kk/2 + 1
            else
              k = ni+1 - (kk-1)/2
            endif
            r1 = r(k-1)
            a1 = a(k-1)
            b1 = b(k-1)
            x1 = r1*sin(a1)*cos(b1)
            y1 = r1*sin(a1)*sin(b1)
            z1 = r1*cos(a1)
            r3 = r(k+1)
            a3 = a(k+1)
            b3 = b(k+1)
            x3 = r3*sin(a3)*cos(b3)
            y3 = r3*sin(a3)*sin(b3)
            z3 = r3*cos(a3)
            dx = x3 - x1
            dy = y3 - y1
            dz = z3 - z1
            x2 = x1 + dx/2
            y2 = y1 + dy/2
            z2 = z1 + dz/2
            r2 = sqrt(x2**2 + y2**2 + z2**2)
            acosa=z2/r2
            if(acosa.LT.-1.)acosa=-1.
            if(acosa.GT.1)acosa=1.
            a2 = acos(acosa)
            sina = sin(a2)
            cosa = cos(a2)
            acosa=x2/r2/sina
            if(acosa.LT.-1.)acosa=-1.
            if(acosa.GT.1)acosa=1.
            b2 = acos(acosa)
            if(y.LT.0.00000)b2=360.00000*dpi-b2
            dn = dx**2 + dy**2 + dz**2
            ddn = sqrt(dn)
            dr = (r3-r1) / ddn
            da = (a3-a1) / ddn
            db = (b3-b1) / ddn
!-- Begin find the gradients and velocities
!-- first find the length of segment
            dseg=sqrt((dx/2)**2+(dy/2)**2+(dz/2)**2)
            ddseg=dseg/2.
! Now ddseg will be a distance to find dV
! along the coordinates
! Determine velocity at 3 points
            v1 = velocity(r1,a1,b1)
            v2 = velocity(r2,a2,b2)
            v3 = velocity(r3,a3,b3)
cc         write(*,*)'v1,v2,v3=',v1,v2,v3
cc         stop
!-- Begin to determine coordinates
!-- of pints surroundibg point a2,b2,r2
!-- at the distance ddseg
            upz = r2+ddseg
            dwz = r2-ddseg
            if(upz.gt.(ro+10.0))then !--- I guess it should be ro+10.0
              upz=ro+10.0
c            if(upz.gt.(ro))then ! make no nagetive r
c              upz=ro
              dwz=upz-dseg
            endif
            if(dwz.le.0.)then
              dwz=0.00000001
!-- set to ro, mistake?
              upz=ro
            endif
!-- The following if-endif is just for P & S, thus comment out for SKS &
!PKP !!!
!-- This gives the lowermost mantle Vp in the outer core
            vr1 = velocity(upz,a2,b2)
            vr2 = velocity(dwz,a2,b2)
            vr=(vr1-vr2)/dseg
            call km2deg(a2,b2,r2,ddseg,RNULL,adV,bdV,rdV)
            vb2 = velocity(rdV,adV,bdV)
            call km2deg(a2,b2,r2,-1.*ddseg,RNULL,adV,bdV,rdV)
            vb1 = velocity(rdV,adV,bdV)
            vb=-1.*(vb1-vb2)/dseg
            call km2deg(a2,b2,r2,RNULL,ddseg,adV,bdV,rdV)
            va2 = velocity(rdV,adV,bdV)
            call km2deg(a2,b2,r2,RNULL,-1.*ddseg,adV,bdV,rdV)
            va1 = velocity(rdV,adV,bdV)
            va=-1.*(va1-va2)/dseg
cc          stop
!-- spherical
!-- velocity gradient
!-- va = va / r2
!-- vb = vb / r2 / sina
!-- (tangential vector) = (slowness vector) / s
            pr = dr
            pa = r2 * da
            pb = r2 * sina * db
            vrd = pr*vr + pa*va + pb*vb
            rvr = vr - vrd*pr
            rva = va - vrd*pa
            rvb = vb - vrd*pb
            rvs = sqrt(rvr*rvr + rva*rva + rvb*rvb)
            if(rvs .eq. 0.) then
              r(k) = r2
              a(k) = a2
              b(k) = b2
            else
              rvr = rvr / rvs
              rva = rva / rvs
              rvb = rvb / rvs
cc           write(*,*)'r,a,b(1)=',r(1),a(1),b(1)
cc           write(*,*)'r,a,b(2)=',r(2),a(2),b(2)
cc           write(*,*)'r,a,b(3)=',r(3),a(3),b(3)
              cc = (1./v1+1./v3)/2.
              rcur = vr*rvr + va*rva + vb*rvb
! Tut esli rcur < 0.0 proishodit hernia
! poetomu postavlen abs. Ne yasno mozhno li eto delat
! ili net no rabotaet. Obichno oshibka poyavliaetsia
! ochen redko v nekotorih tochkah
! v etom sluchae abs prosto ne daet oshibki y posledniaya iteraciya
! uzhe ne imeet rcur negativnim y podgoniaet normalno reshenie
! ( mozhet bit)
              if(rcur.LE.0.0)then
                write(*,*)'Negative'
                rcur=abs(rcur)
              endif
              rcur = (cc*v2+1.) / (4.*cc*rcur)
              rcur = -1*rcur + sqrt(rcur**2+dn/(8.*cc*v2))
              rdr = rvr * rcur
              rda = rva * rcur
              rdb = rvb * rcur
              rpr = r2 + rdr
              ap = a2 + rda/r2
              bp = b2 + rdb/(r2*sina)
              r(k) = (rpr-r(k))*xfac + r(k)
! if r(k)>6371 then force it to the surface.
              if (r(k).gt.(ro+10.0)) r(k)=ro+10.0
              a(k) = (ap-a(k))*xfac + a(k)
              b(k) = (bp-b(k))*xfac + b(k)
            endif
          enddo
          idstn=ni
          do j=1,ni+1
            w(1,j) = r(j)
            w(2,j) = a(j)
            w(3,j) = b(j)
          enddo
          ni=idstn
cc          write(*,*) 'ni,r123,a123,b123=',ni,r(1:3),a(1:3),b(1:3)
          tk = rtim(ni+1,r,a,b)
          if(abs(to-tk) .le. to*flim) go to 310
          to = tk
        enddo
 310    continue
        to=tk
!-- skip increasing of segment number if minimum length
!-- of segment is exceed or maximum number of segments
!-- was reached
        if(dseg.lt.mins.or.ni.ge.n2) then
          igood=1
          go to 66666
        endif
!-- double the number of points.
        ni = ni * 2
        do i=1,ni/2+1
          r(i*2-1) = w(1,i)
          a(i*2-1) = w(2,i)
          b(i*2-1) = w(3,i)
        enddo
        do k=2,ni,2
          r1 = r(k-1)
          a1 = a(k-1)
          b1 = b(k-1)
          x1 = r1*sin(a1)*cos(b1)
          y1 = r1*sin(a1)*sin(b1)
          z1 = r1*cos(a1)
          r3 = r(k+1)
          a3 = a(k+1)
          b3 = b(k+1)
          x3 = r3*sin(a3)*cos(b3)
          y3 = r3*sin(a3)*sin(b3)
          z3 = r3*cos(a3)
          dx = x3 - x1
          dy = y3 - y1
          dz = z3 - z1
          x2 = x1 + dx/2
          y2 = y1 + dy/2
          z2 = z1 + dz/2
          r2 = sqrt(x2**2 + y2**2 + z2**2)
          acosa=z2/r2
          if(acosa.LT.-1.)acosa=-1.
          if(acosa.GT.1)acosa=1.
          a2 = acos(acosa)
          sina = sin(a2)
          acosa=x2/r2/sina
          if(acosa.LT.-1.)acosa=-1.
          if(acosa.GT.1)acosa=1.
          b2 = acos(acosa)
          if(y.LT.0.00000)b2=360.00000*dpi-b2
          r(k) = r2
          a(k) = a2
          b(k) = b2
        enddo
        tk = rtim(ni+1,r,a,b)
!-- here i change tp and put to
        if(abs(to-tk) .le. to*flim) then
          igood=1
          go to 99999
        endif
        to = tk
      enddo
99999 continue
cc      write(*,*) '-----------------'
cc      stop
      idstn=ni
      do i=1,ni+1
        w(1,i) = r(i)
        w(2,i) = a(i)
        w(3,i) = b(i)
      enddo
      ni=idstn
66666 continue
!-- Return coordinates to the origin
      idstn=ni
c        write(*,*)'ni=',ni
      do k=1,ni+1
        w(1,k) = ro-w(1,k)
        w(2,k) = w(2,k)*r2d
        w(2,k) = geoc_to_geog(90.0-w(2,k))
        w(3,k) = w(3,k)*r2d+shiftlo
        if(w(3,k).lt.0.)w(3,k)=360.+w(3,k)
      enddo
ccc output route of first arrival ray, JYHuang 20200423
      open(56,file='out_route')
      write(56,*)'line 1:n_route,s_treval_time, 2:source, 3:sta, 4:dep_f
     &rom_source, 5:lat, 6:lon'
      write(56,*) ni+1,tk
      write(56,*)evla,evlo,evdp
      write(56,*)stla,stlo,stel
      write(56,*)w(1,1:ni+1)
      write(56,*)w(2,1:ni+1)
      write(56,*)w(3,1:ni+1)
ccc
      ni=idstn
      np=ni+1
!-- convert ray point to cartesion coord.
      return
      end subroutine pbr
      subroutine km2deg(ala,alo,adp,dx,dy,bla,blo,bdp)
      implicit real*8(a-h,o-z)
      real*8 ala,alo,adp,dx,dy,bla,blo,bdp
      real*8 dpi,dps
!-- This subroutine calculate position of new point
!-- in polar coordinates basing on the coordinates
!-- of main point in radians ( la is colatitude) and dx and dy in
!kilometers
      dpi = asin(1.)/ 90.
      dps=adp*SIN(ala)
      blo=alo+atan2(dx,dps)
      bla=ala+atan2(dy,adp)
      if(bla.gt.(180.*dpi))then
        bla=360.*dpi-bla
        blo=blo+180.*dpi
      endif
      if(bla.lt.0.)then
        bla=abs(bla)
        blo=blo+180.*dpi
      endif
      if(blo.lt.0.)blo=360.*dpi+blo
      if(blo.gt.(360.*dpi))blo=blo-(360.*dpi)
      bdp=sqrt(adp**2+dx**2+dy**2)
      return
      end subroutine km2deg
      function rtim(m, r, a, b)
      implicit real*8(a-h,o-z)
      real*8 x1,y1,z1,x2,y2,z2,dl
      real*8 rv2,sm,rv1,rtim
      parameter (msg = 16384)
      real*8 r(msg+1), a(msg+1), b(msg+1)
      integer m
      if(m.GT.(msg+1))write(*,*)'*'
      rtim = 0.
      rv1 = 1./velocity(r(1),a(1),b(1))
      do j=1,m-1
        x1 = r(j)*sin(a(j))*cos(b(j))
        y1 = r(j)*sin(a(j))*sin(b(j))
        z1 = r(j)*cos(a(j))
        x2 = r(j+1)*sin(a(j+1))*cos(b(j+1))
        y2 = r(j+1)*sin(a(j+1))*sin(b(j+1))
        z2 = r(j+1)*cos(a(j+1))
        dl = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
        rv2 = 1./velocity(r(j+1),a(j+1),b(j+1))
        sm = (rv1 + rv2) / 2.
        rtim = rtim + sqrt(dl)*sm
        rv1 = rv2
      enddo
      end function rtim
      real*8 function earthr(xlat)
! this routine establishes the short distance conversion factors
! given the origin of coordinates
! the rotation angle is converted to radians also
! common block variables:
! local variables:
      double precision dlt1,dxlt,drad,drlt,xlat
      data re/6378.163/, ell/298.26/
      drad=1.7453292d-2
      drlt=9.9330647d-1
      dxlt=dble(xlat*60.0)
! conversion factor for latitude
      dlt1=datan(drlt*dtan(dxlt*drad/60.d0))
      earthr=re*(1.0-sngl(dsin(dlt1)**2)/ell)
      end function earthr
      



