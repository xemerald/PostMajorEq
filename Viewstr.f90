!############################################################################!
!filelist:
!         Revised_SMA_catalog-08-03-2014.filename
!         Revised_SMA_catalog-08-03-2014.evt
!         Revised_SMA_catalog-08-03-2014.stn
!         Revised_SMA_catalog-08-03-2014.evt_stn
!         list.txt
!
!############################################################################!
program ASCIIfile
character station_ID*4, filename*28(1:3), path*60,testmark*1, temp*10, date*8, eve_date*8, stname*4
real acc(32768,3), x(32768), y(32768), z(32768), sec, dt, parrive, sarrive, p_wei, s_wei
real vcc(32768,3)
logical file_exist
integer i, j, k, nps,stnum, m, n, ren, isps,kk
real pga, pga_time, pgv, pgv_time, pgd, pgd_time, sn, pa_aft, pa_bef,pd35,pga80,pd2,pga4
real pv(10), pa(10), pd(10),lon,lat
character st*3(20), filepath*120
integer event_num,sta_num,sta_sum, pd_test_type
integer sta_num_all
character yy*4, mm*2,dd*2 ,hh*2, min*2, ss*2, depth, st_name*3
character sta_date*8, sta_yy*4, sta_mm*2, sta_dd, sta_hh*2, sta_min*2, sta_ss*2, event_num_C*3, test_type_C*1
integer year, month, day, hour, minute
real second
integer st_year, st_month, st_day, st_hour, st_minute,test
real st_second
real rmark,tc(5),v_i2,d_i2,pvi2(5),pdi2(5),tcavg(5),tc_t(5)

real Magn, Evt_dept, Evt_lon, Evt_l, C1(1:5), C2(1:5), C3(1:5), tt, cwbm,r, ave_md
character info_temp*100,info*100(1:1000),info_t*100(1:1000)
real       MW_SUN, MW_s, MW_d,pi,t_temp,ttime(600),tt2,tt3,t_time2,t_time3,MW
character fn*128,stn*6,cdate*23,card*128,temp_data*200,filenameZ*90,filenameN*90,filenameE*90
character A_event_ID*2,A_event*200, data_test_point*1

integer istat,PGOPEN,np,event_ID,event_ID2,x_n,k_s

real Magn_s,SDV_s,Magn_all,Magn_all2
real leadtime

pi=4.0*atan(1.0)      

sta_sum=0
event_num=0
isps=100
dt=1.0/100


open(10,file="stalist.txt")
!open(11,file="all_use_events.txt")
open(12,file="dataoutput.txt")
open(13,file="lead.txt")

!call readsac(fn,sa,nps,dt,card)
read(10,*,END=10)Evt_l,Evt_lon
do while(.true.)
92  continue
    read(10,*,END=10)station_ID, lat, lon
    write(*,*)station_ID, lat, lon

87 continue


    filenameZ=TRIM(station_ID)//".HLZ.TW.--"
    filenameN=TRIM(station_ID)//".HLN.TW.--"
    filenameE=TRIM(station_ID)//".HLE.TW.--"

	inquire(file=filenameZ, exist=file_exist) 
    if (file_exist .EQ. .true. ) then

       call readsac(filenameZ,acc(:,3),nps,dt,card)
       call readsac(filenameN,acc(:,2),nps,dt,card)
       call readsac(filenameE,acc(:,1),nps,dt,card)
            
            do i=1,3
                do j=1,nps
	                if(acc(j,i) .eq. -12345.0)then
                        acc(j,i) = 0.0
	                endif
	            enddo
            enddo

       acc(:,3)=acc(:,3)/4303.0
       acc(:,2)=acc(:,2)/4303.0
       acc(:,1)=acc(:,1)/4303.0

                 if (nps .GT. 32768) nps=32768
                    if (nps .LE. 1000) then
                        write(*,*)"nps<1000"
                        !pause
                        goto 92
                    end if
                    is1=100
                    parrive=0
                    write(*,*)"autopick",nps,dt,is1
                    ren=0
                    call autopicking(acc,nps,dt,parrive,sarrive,p_wei,s_wei,is1,ren)

                    if(ren==1) write(4,*)eve_date,"pick err"

                    if (parrive>=100)then
                        write(*,*)"!!!!!!!!!!!",parrive
                        goto 90

                    end if

                    if (parrive<100)then
                        do i=1,nps
        
                            if (acc(i,1)**2>=1)then
                                parrive=i*dt
                                write(*,*)"nps=",nps," dt=",dt," is1",is1      
                                print*,'P S arrivals : ',parrive
                                write(4,*)eve_date,"pga**2>1.0"
                                if (parrive>=100)then
                                    goto 90
                                else
                                    goto 92
                                end if
                            end if

                        end do

                        write(4,*)eve_date,"no parr"

                        goto 92
                    end if
    		    else
                    write(4,*)eve_date,"have no file"
                    goto 92
                end if

40              continue

!#################################################################
90          close(2)

!-------------------------------------------------------------------------------------------
!-------------------------------------≠p∫‚PGA PGV PGD---------------------------------------
!-------------------------------------------------------------------------------------------

9528 continue
            pgv=0.0
            pgv_time=0.0
            pga=0.0
            pga_time=0.0
            pgd=0.0
            pgd_time=0.0
            iparr=parrive/dt

!-------------------------------------------------------------------------------------------








!-- Find PGA & PGA Time
            do i=1,3
                do j=iparr,nps
	                if((abs(acc(j,i)).gt.pga) .and. (abs(acc(j,i)).ne.12345.0))then
                        pga=abs(acc(j,i))
		                pga_time=j*dt
	                endif
	            enddo
            enddo

            pga4=600.0
            do i=1,3
                do j=iparr,nps
	                if(abs(acc(j,i)).gt. 4.0)then
                        if((j*dt).lt.pga4) pga4=j*dt
                        exit
	                endif
	            enddo
            enddo

            pga80=600.0
            do i=1,3
                do j=iparr,nps
	                if(abs(acc(j,i)).gt. 80.0)then
		                if ((j*dt).lt.pga80) pga80=j*dt
                        exit
	                endif
	            enddo
            enddo


!-- Find Pa1 to Pa10
556            pa=0.0
            do i=1,10
                do j=iparr,iparr+isps*i
	                if(abs(acc(j,1)).gt.pa(i)) pa(i)=abs(acc(j,1))
	            enddo
            enddo
            goto 334
            !--------!


 334        do i=1,3
	            y=acc(:,i)
	            call rmbsl(y,nps,200)
                !call remove_linear(np,y,dt)
                z(1)=0.
                do j=2,nps
                    z(j)=z(j-1)+(y(j)+y(j-1))*dt/2.
                enddo
	            call IIRFILT(z, nps,"BU      ",2,"HP      ",0.075,50.0,dt,1)
	            acc(:,i)=z
            enddo
           



!-- Find PGV & PGV Time
             do i=1,3
                do j=iparr,nps
	                if(abs(acc(j,i)).gt.pgv)then
	                    pgv=abs(acc(j,i))
		                pgv_time=j*dt
	                endif
	            enddo
            enddo
                    !--find pv**2©M
                v_i2=0
            do i=1,1000
                v_i2=v_i2+(((acc(i,1))**2+(acc(i,2))**2+(acc(i,3))**2))**0.5
                if(i==200)then
                    pvi2(1)=v_i2
                else if(i==400) then
                    pvi2(2)=v_i2
                else if(i==600) then
                    pvi2(3)=v_i2
                else if(i==800) then
                    pvi2(4)=v_i2
                else if(i==1000) then
                    pvi2(5)=v_i2
                end if

            end do



!-- Find Pv1 to Pv10
            pv=0.0
            do i=1,10
                do j=iparr,iparr+isps*i
	                if(abs(acc(j,1)).gt.pv(i)) pv(i)=abs(acc(j,1))
	            enddo
            enddo

            do i=1,3
	            y=acc(:,i)
                z(1)=0.
                do j=2,nps
                    z(j)=z(j-1)+(y(j)+y(j-1))*dt/2.
                enddo
	            call IIRFILT(z, nps,'BU      ',2,'HP      ',0.075,50.0,dt,1)
	            acc(:,i)=z
            enddo
    
!-- Find PGD & PGD Time
            do i=1,3
                do j=iparr,nps
	                if(abs(acc(j,i)).gt.pgd)then
	                    pgd=abs(acc(j,i))
		                pgd_time=j*dt
	                endif
	            enddo
            enddo

            pd2=600.0
            do i=1,3
                do j=iparr,nps
	                if(abs(acc(j,i)).gt.0.2)then
                        if((j*dt).lt.pd2)pd2=j*dt
                        exit
	                endif
	            enddo
            enddo

            pd35=600.0
            do i=1,3
                do j=iparr,nps
	                if(abs(acc(j,i)).gt.0.35)then
                        if((j*dt).lt.pd35)pd35=j*dt
                        exit
	                endif
	            enddo
            enddo


            !--find pd**2©M
557                d_i2=0
            do i=1,1000
                d_i2=d_i2+(((acc(i,1))**2+(acc(i,2))**2+(acc(i,3))**2))**0.5
                if(i==200)then
                    pdi2(1)=d_i2
                else if(i==400) then
                    pdi2(2)=d_i2
                else if(i==600) then
                    pdi2(3)=d_i2
                else if(i==800) then
                    pdi2(4)=d_i2
                else if(i==1000) then
                    pdi2(5)=d_i2
                end if

            end do

!-- Find Pd1 to Pd10
            pd=0.0
            do i=1,10
                do j=iparr,iparr+isps*i
	                if(abs(acc(j,1)).gt.pd(i)) pd(i)=abs(acc(j,1))
	            enddo
            enddo

91          continue
            close(2)

9529        continue

!--------------------- ≠p∫‚æ_∑Ω∂Z&æ_•°∂Z----------------------------------------



            call cal_delta(Evt_l,Evt_lon,lat,lon,sta_r) !æ_•°∂Z
            r=SQRT(sta_r**2+Evt_dept**2)

            if (pd2 < pga4)then
                if (pd35 > pga80)then
                        leadtime = pga80-pd2
                else
                        leadtime = pd35-pd2
                end if
            else
                if (pd35 > pga80)then
                        leadtime = pga80-pga4
                else
                        leadtime = pd35-pga4
                end if
            end if
            if (pd2<=0 .and. pga4<=0)leadtime = 0
            if (leadtime <=0)leadtime = 0
            !if ((pgv_time-pd35)<=(pgv_time-pga80))then
             !       leadtime = pgv_time-pga80
            !else
             !       leadtime = pgv_time-pd35
            !end if
            !if (pd35<=0 .and. pga80<=0)leadtime = 0
            !if (leadtime <=0)leadtime = 0

            write(12,"(a4,1x,f6.2,1x,f8.3,1x,f8.3,1x,f8.3,1x,e10.4,1x,f8.3,1x,e10.4,1x,f8.3,1x,e10.4,1x,f8.3,5(1x,e10.4),10(1x,e10.4))")station_ID,r,parrive,pd35,pga80,pga,pga_time,pgv,pgv_time,pgd,pgd_time,pa(1:5),pv(1:5),pd(1:5)
            !write(12,"(a4,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f8.3,1x,f8.3,1x,f8.3,1x,e10.4,1x,f8.3,1x,e10.4,1x,f8.3,1x,e10.4,1x,f8.3,5(1x,e10.4),10(1x,e10.4))")station_ID,lat,lon,leadtime,r,parrive,pd35,pga80,pga,pga_time,pgv,pgv_time,pgd,pgd_time,pa(1:5),pv(1:5),pd(1:5)
            if ( r .le. 50.0) then
                !write(13,"(a4,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2)")station_ID,lon,lat,leadtime,pd2,pga4,pd35,pga80
                !write(12,"(a4,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.2)")station_ID,lat,lon,pgv,pd(3)
            end if
            write(13,"(a4,1x,f6.2,1x,f6.2,1x,f8.3)")station_ID,lon,lat,pgv
95            continue

60  continue
end do

10	continue

558 write(*,*)"Pd done"
write(*,*)"dt=",dt,parrive
close(12)
close(13)
close(10)
pause
end program


!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################
!###################################################################################


SUBroutine remove_linear(ndata,y,dt)
      real*8 sx,sy,sxx,sxy,det,tc,a,b
      real x(32768),Y(32768)
!--------------------------------------------------------------------
      do 100 id=1,ndata
100   x(id)=dt*(float(id)-1.)
!-------------------------------------------------------------------------
      SXX=0.
      SXY=0.
      SX=0.
      SY=0.
      do 500 id=1,Ndata
      sxx=sxx+X(id)*x(id)
      SX=SX+X(ID)
      SY=SY+Y(ID)
      sxy=sxy+x(id)*y(id)
500   continue
      tc=float(ndata)
      DET=SXX*Tc-SX*SX
      A=(SXY*Tc-SX*SY)/DET
      B=(SXX*SY-SX*SXY)/DET

      do 510 id=1,ndata
      y(id)=y(id)-B-A*x(id)
510   continue
      return
end	SUBroutine remove_linear

subroutine rmbsl(x,n,nb)
!     remove base line defined by the first nb points
!     original trace is altered
  real x(n)
  s=0.
  do i=1,nb
    s=s+x(i)
  end do
  s=s/float(nb)
  do i=1,n
    x(i)=x(i)-s
  end do
end subroutine rmbsl 




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


