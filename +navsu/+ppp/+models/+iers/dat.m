%
function [deltat,j]=dat(iy,im,id,fd)
%+
%  - - - -
%   D A T
%  - - - -
%
%  For a given UTC date, calculate delta(AT) = TAI-UTC.
%
%     :------------------------------------------:
%     :                                          :
%     :                 IMPORTANT                :
%     :                                          :
%     :  A new version of this routine must be   :
%     :  produced whenever a new leap second is  :
%     :  announced.  There are five items to     :
%     :  change on each such occasion:           :
%     :                                          :
%     :  1) The parameter NDAT must be           :
%     :     increased by 1.                      :
%     :                                          :
%     :  2) The set of DATA statements that      :
%     :     initialize the arrays IDAT and       :
%     :     DATS must be extended by one line.   :
%     :                                          :
%     :  3) The parameter IYV must be set to     :
%     :     the current year.                    :
%     :                                          :
%     :  4) The 'Latest leap second' comment     :
%     :     below must be set to the new leap    :
%     :     second date.                         :
%     :                                          :
%     :  5) The 'This revision' comment, later,  :
%     :     must be set to the current date.     :
%     :                                          :
%     :  Change (3) must also be carried out     :
%     :  whenever the routine is re-issued,      :
%     :  even if no leap seconds have been       :
%     :  added.                                  :
%     :                                          :
%     :  Latest leap second:  2016 December 31   :
%     :                                          :
%     :__________________________________________:
%
%  This routine is part of the International Astronomical Union's
%  SOFA (Standards of Fundamental Astronomy) software collection.
%
%  Status:  support routine.
%
%  Given:
%     IY       i     UTC:  year (Notes 1 and 2)
%     IM       i           month (Note 2)
%     ID       i           day (Notes 2 and 3)
%     FD       d           fraction of day (Note 4)
%
%  Returned:
%     DELTAT   d     TAI minus UTC, seconds
%     J        i     status (Note 5):
%                       1 = dubious year (Note 1)
%                       0 = OK
%                      -1 = bad year
%                      -2 = bad month
%                      -3 = bad day (Note 3)
%                      -4 = bad fraction (Note 4)
%                      -5 = internal error (Note 5)
%
%  Notes:
%
%  1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
%     to call the routine with an earlier date.  If this is attempted,
%     zero is returned together with a warning status.
%
%     Because leap seconds cannot, in principle, be predicted in
%     advance, a reliable check for dates beyond the valid range is
%     impossible.  To guard against gross errors, a year five or more
%     after the release year of the present routine (see parameter IYV)
%     is considered dubious.  In this case a warning status is returned
%     but the result is computed in the normal way.
%
%     For both too-early and too-late years, the warning status is J=+1.
%     This is distinct from the error status J=-1, which signifies a
%     year so early that JD could not be computed.
%
%  2) If the specified date is for a day which ends with a leap second,
%     the UTC-TAI value returned is for the period leading up to the
%     leap second.  If the date is for a day which begins as a leap
%     second ends, the UTC-TAI returned is for the period following the
%     leap second.
%
%  3) The day number must be in the normal calendar range, for example
%     1 through 30 for April.  The 'almanac' convention of allowing
%     such dates as January 0 and December 32 is not supported in this
%     routine, in order to avoid confusion near leap seconds.
%
%  4) The fraction of day is used only for dates before the introduction
%     of leap seconds, the first of which occurred at the end of 1971.
%     It is tested for validity (0 to 1 is the valid range) even if not
%     used;  if invalid, zero is used and status J=-4 is returned.  For
%     many applications, setting FD to zero is acceptable;  the
%     resulting error is always less than 3 ms (and occurs only
%     pre-1972).
%
%  5) The status value returned in the case where there are multiple
%     errors refers to the first error detected.  For example, if the
%     month and day are 13 and 32 respectively, J=-2 (bad month) will be
%     returned.  The 'internal error' status refers to a case that is
%     impossible but causes some compilers to issue a warning.
%
%  6) In cases where a valid result is not available, zero is returned.
%
%  References:
%
%  1) For dates from 1961 January 1 onwards, the expressions from the
%     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
%
%  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
%     the 1992 Explanatory Supplement.
%
%  Called:
%     CAL2JD   Gregorian calendar to JD
%
%  This revision:  2016 July 7
%
%  SOFA release 2016-05-03
%
%  Copyright (C) 2016 IAU SOFA Board.  See notes at end.
%
%-----------------------------------------------------------------------



%  Release year for this version of iau_DAT

% persistent da dats djm djm0 drift firstCall idat is iyv js m more n ndat nera1 ; if isempty(firstCall),firstCall=1;end;

iyv=2016 ;

%  Number of Delta(AT) changes (increase by 1 for each new leap second)
ndat=42 ;

%  Number of Delta(AT) expressions before leap seconds were introduced
nera1=14 ;

%  Dates (year, month) on which new Delta(AT) came into force

%  New Delta(AT) which came into force on the given dates
dats=zeros(1,ndat);

%  Reference dates (MJD) and drift rates (s/day), pre leap seconds
drift=zeros(2,nera1);

%  Miscellaneous local variables
more=false;
js=0;
m=0;
n=0;
is=0;
da=0;
djm0=0;
djm=0;

%  Dates, Delta(AT)s, reference dates, and drift rates
%  if firstCall, [ idat([1:2],[1:14]),dats([1:14]),drift([1:2],[1:14])]=reshape([1960,1,1.4178180d0,37300d0,0.001296d0,1961,1,1.4228180d0,37300d0,0.001296d0,1961,8,1.3728180d0,37300d0,0.001296d0,1962,1,1.8458580d0,37665d0,0.0011232d0,1963,11,1.9458580d0,37665d0,0.0011232d0,1964,1,3.2401300d0,38761d0,0.001296d0,1964,4,3.3401300d0,38761d0,0.001296d0,1964,9,3.4401300d0,38761d0,0.001296d0,1965,1,3.5401300d0,38761d0,0.001296d0,1965,3,3.6401300d0,38761d0,0.001296d0,1965,7,3.7401300d0,38761d0,0.001296d0,1965,9,3.8401300d0,38761d0,0.001296d0,1966,1,4.3131700d0,39126d0,0.002592d0,1968,2,4.2131700d0,39126d0,0.002592d0],2,ndat); end;
%
% %  Dates and Delta(AT)s
%  if firstCall, [ idat([1:2],[15:30]),dats([15:30])]=reshape([1972,1,10d0,1972,7,11d0,1973,1,12d0,1974,1,13d0,1975,1,14d0,1976,1,15d0,1977,1,16d0,1978,1,17d0,1979,1,18d0,1980,1,19d0,1981,7,20d0,1982,7,21d0,1983,7,22d0,1985,7,23d0,1988,1,24d0,1990,1,25d0],2,ndat); end;
%
%  if firstCall, [ idat([1:2],[31:ndat]),dats([31:ndat])]=reshape([1991,1,26d0,1992,7,27d0,1993,7,28d0,1994,7,29d0,1996,1,30d0,1997,7,31d0,1999,1,32d0,2006,1,33d0,2009,1,34d0,2012,7,35d0,2015,7,36d0,2017,1,37d0],2,ndat); end;
% firstCall=0;
datafull =  [1960,  1,  1.4178180D0, 37300D0, 0.001296D0;
    1961,  1,  1.4228180D0, 37300D0, 0.001296D0;
    1961,  8,  1.3728180D0, 37300D0, 0.001296D0;
    1962,  1,  1.8458580D0, 37665D0, 0.0011232D0;
    1963, 11,  1.9458580D0, 37665D0, 0.0011232D0;
    1964,  1,  3.2401300D0, 38761D0, 0.001296D0;
    1964,  4,  3.3401300D0, 38761D0, 0.001296D0;
    1964,  9,  3.4401300D0, 38761D0, 0.001296D0;
    1965,  1,  3.5401300D0, 38761D0, 0.001296D0;
    1965,  3,  3.6401300D0, 38761D0, 0.001296D0;
    1965,  7,  3.7401300D0, 38761D0, 0.001296D0;
    1965,  9,  3.8401300D0, 38761D0, 0.001296D0;
    1966,  1,  4.3131700D0, 39126D0, 0.002592D0;
    1968,  2,  4.2131700D0, 39126D0, 0.002592D0;
    1972,  1, 10D0 nan nan ;
    1972,  7, 11D0 nan nan ;
    1973,  1, 12D0 nan nan ;
    1974,  1, 13D0 nan nan ;
    1975,  1, 14D0 nan nan ;
    1976,  1, 15D0 nan nan ;
    1977,  1, 16D0 nan nan ;
    1978,  1, 17D0 nan nan ;
    1979,  1, 18D0 nan nan ;
    1980,  1, 19D0 nan nan ;
    1981,  7, 20D0 nan nan ;
    1982,  7, 21D0 nan nan ;
    1983,  7, 22D0 nan nan ;
    1985,  7, 23D0 nan nan ;
    1988,  1, 24D0 nan nan ;
    1990,  1, 25D0 nan nan ;
    1991,  1, 26D0 nan nan ;
    1992,  7, 27D0 nan nan ;
    1993,  7, 28D0 nan nan ;
    1994,  7, 29D0 nan nan ;
    1996,  1, 30D0 nan nan ;
    1997,  7, 31D0 nan nan ;
    1999,  1, 32D0 nan nan ;
    2006,  1, 33D0 nan nan ;
    2009,  1, 34D0 nan nan ;
    2012,  7, 35D0 nan nan ;
    2015,  7, 36D0 nan nan ;
    2017,  1, 37D0  nan nan ];
idat = datafull(:,1:2)';
dats = datafull(:,3)';
drift = datafull(1:14,4:5)';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%  Initialize the result to zero and the status to OK.
da = 0d0;
js = 0;

%  If invalid fraction of a day, set error status and give up.
if( fd<0d0 || fd>1d0 )
    js = -4;
    goto 100;
end

%  Convert the date into an MJD.
% [iy,im,id,djm0,djm,js]=iau_cal2jd(iy,im,id,djm0,djm,js);
[djm0,djm,js]=navsu.ppp.models.iers.cal2jd_iers(iy,im,id);

%  If invalid year, month, or day, give up.
if( js>=0 )
    
    %  If pre-UTC year, set warning status and give up.
    if( iy<idat(1,1) )
        js = 1;
        goto 100;
    end
    
    %  If suspiciously late year, set warning status but proceed.
    if( iy>iyv+5 )
        js = 1;
    end
    
    %  Combine year and month.
    m = fix(12.*iy + im);
    
    %  Find the most recent table entry.
    is = 0;
    more = true;
    for n = ndat : -1: 1
        if( more )
            is = fix(n);
            more = m<(12.*idat(1,n)+idat(2,n));
        end
    end; n = fix(1 + -1);
    
    %  Prevent underflow warnings.
    if( is<1 )
        js = -5;
        goto 100;
    end
    
    %  Get the Delta(AT).
    da = dats(is);
    
    %  If pre-1972, adjust for drift.
    if( is<=nera1 )
        da = da +(djm+fd-drift(1,is)).*drift(2,is);
    end
end

%  Return the Delta(AT) value and the status.
deltat = da;
j = fix(js);

%  Finished.

%+----------------------------------------------------------------------
%
%  Copyright (C) 2016
%  Standards Of Fundamental Astronomy Board
%  of the International Astronomical Union.
%
%  =====================
%  SOFA Software License
%  =====================
%
%  NOTICE TO USER:
%
%  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING SIX TERMS AND
%  CONDITIONS WHICH APPLY TO ITS USE.
%
%  1. The Software is owned by the IAU SOFA Board ('SOFA').
%
%  2. Permission is granted to anyone to use the SOFA software for any
%     purpose, including commercial applications, free of charge and
%     without payment of royalties, subject to the conditions and
%     restrictions listed below.
%
%  3. You (the user) may copy and distribute SOFA source code to others,
%     and use and adapt its code and algorithms in your own software,
%     on a world-wide, royalty-free basis.  That portion of your
%     distribution that does not consist of intact and unchanged copies
%     of SOFA source code files is a 'derived work' that must comply
%     with the following requirements:
%
%     a) Your work shall be marked or carry a statement that it
%        (i) uses routines and computations derived by you from
%        software provided by SOFA under license to you; and
%        (ii) does not itself constitute software provided by and/or
%        endorsed by SOFA.
%
%     b) The source code of your derived work must contain descriptions
%        of how the derived work is based upon, contains and/or differs
%        from the original SOFA software.
%
%     c) The names of all routines in your derived work shall not
%        include the prefix 'iau' or 'sofa' or trivial modifications
%        thereof such as changes of case.
%
%     d) The origin of the SOFA components of your derived work must
%        not be misrepresented;  you must not claim that you wrote the
%        original software, nor file a patent application for SOFA
%        software or algorithms embedded in the SOFA software.
%
%     e) These requirements must be reproduced intact in any source
%        distribution and shall apply to anyone to whom you have
%        granted a further right to modify the source code of your
%        derived work.
%
%     Note that, as originally distributed, the SOFA software is
%     intended to be a definitive implementation of the IAU standards,
%     and consequently third-party modifications are discouraged.  All
%     variations, no matter how minor, must be explicitly marked as
%     such, as explained above.
%
%  4. You shall not cause the SOFA software to be brought into
%     disrepute, either by misuse, or use for inappropriate tasks, or
%     by inappropriate modification.
%
%  5. The SOFA software is provided 'as is' and SOFA makes no warranty
%     as to its use or performance.   SOFA does not and cannot warrant
%     the performance or results which the user may obtain by using the
%     SOFA software.  SOFA makes no warranties, express or implied, as
%     to non-infringement of third party rights, merchantability, or
%     fitness for any particular purpose.  In no event will SOFA be
%     liable to the user for any consequential, incidental, or special
%     damages, including any lost profits or lost savings, even if a
%     SOFA representative has been advised of such damages, or for any
%     claim by any third party.
%
%  6. The provision of any version of the SOFA software under the terms
%     and conditions specified herein does not imply that future
%     versions will also be made available under the same terms and
%     conditions.
%
%  In any published work or commercial product which uses the SOFA
%  software directly, acknowledgement (see www.iausofa.org) is
%  appreciated.
%
%  Correspondence concerning SOFA software should be addressed as
%  follows:
%
%      By email:  sofa@ukho.gov.uk
%      By post:   IAU SOFA Center
%                 HM Nautical Almanac Office
%                 UK Hydrographic Office
%                 Admiralty Way, Taunton
%                 Somerset, TA1 2DN
%                 United Kingdom
%
%-----------------------------------------------------------------------

end

