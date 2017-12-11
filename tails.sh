#! /bin/csh

if ( `hostname` == brody-t2e ) then
    set hostnumber = "t2e"
else
    set hostnumber = `hostname | sed s/proanti//`
endif


set j=1
foreach f( `find ../Reports$hostnumber -name "report_*"` )
      grep "] [0-9][0-9]:" $f | tail -1
end



