if length(ARGS)<2
   println("Usage: remoteInspector  machinenum  belowlevel")
end

global mnum, bnum
try
   global mnum = parse(String, ARGS[1])
catch
   global mnum = "006"
end

try
   global bnum = parse(Int64, ARGS[2]);
catch
   global bnum = 4
end

remoteCmd = " \" cd A_REWRITE ; grep BELOW ju*\" "
ans = read(`gcloud compute ssh --command $remoteCmd carlosbrody@proanti006`, String)
