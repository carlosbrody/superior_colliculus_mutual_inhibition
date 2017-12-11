#! /bin/csh

set j=1
while ( $j <= "30" )
      grep "] [0-9][0-9]:" Reports029/report_out_$j | tail -1
      @ j++
end



