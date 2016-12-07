
  SUBROUTINE sort_monotonic(lonin, latin, fieldin, lonout, latout, fieldout, npx, npy)

  IMPLICIT NONE

  REAL(8), DIMENSION(npy,npx), INTENT(in)      :: lonin, latin, fieldin
  REAL(8), DIMENSION(npy,npx), INTENT(out)     :: lonout, latout, fieldout
  REAL(8), DIMENSION(npx)                      :: ztmp, ztmplon, ztmplat, ztmpfield
  REAL(8), DIMENSION(npy,npx)                  :: zdist
  INTEGER, INTENT(in)                          :: npx, npy
  INTEGER                                      :: ji, jj, ibr, nmibr
  INTEGER, DIMENSION(2)                        :: jiind


  DO jj=1,npy
     ! continuity breaks at :
     ibr = MINLOC(lonin(jj,:),1)
     nmibr = npx - ibr + 1

     lonout(jj,1:nmibr)     = lonin(jj,ibr:npx)
     lonout(jj,nmibr+1:npx) = lonin(jj,1:ibr-1)

     latout(jj,1:nmibr)     = latin(jj,ibr:npx)
     latout(jj,nmibr+1:npx) = latin(jj,1:ibr-1)

     fieldout(jj,1:nmibr)     = fieldin(jj,ibr:npx)
     fieldout(jj,nmibr+1:npx) = fieldin(jj,1:ibr-1)

  ENDDO

  END SUBROUTINE

