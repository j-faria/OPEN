      FUNCTION julday(mm,id,iyyy)
C julday returns the Julian Day Number that begins at noon of the calendar date specified by 
C month mm, day id, and year iyyy, all integer variables. Positive year signifies A.D.; 
C negative, B.C. Remember that the year after 1 B.C. was 1 A.D.
      INTEGER julday,id,iyyy,mm,IGREG
C Gregorian Calendar adopted Oct. 15, 1582.
      PARAMETER (IGREG=15+31*(10+12*1582))
      INTEGER ja,jm,jy
      jy=iyyy
      if (jy.eq.0) stop 'julday: there is no year zero'
      if (jy.lt.0) jy=jy+1
      if (mm.gt.2) then
        jm=mm+1
      else
        jy=jy-1
        jm=mm+13
      endif
      julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
C Test whether to change to Gregorian Calendar
      if (id+31*(mm+12*iyyy).ge.IGREG) then
        ja=int(0.01*jy)
        julday=julday+2-ja+int(0.25*ja)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .