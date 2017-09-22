

function lockFile(filename)

    sleep 1

    lockname=filename*".lock";
    fnames = readdir()

    matched_filenames = Array{Bool}(length(fnames))
    matched_filenames2 = Array{Bool}(length(fnames))
    for i=1:length(fnames)
        matched_filenames[i] = ismatch(Regex(@sprintf("^%s", filename)), fnames[i])
        matched_filenames2[i] = ismatch(Regex(@sprintf("^%s", lockname)), fnames[i])
    end
    mynum = length(find(matched_filenames))
    if mynum>0
        locked=0
        return locked
    end
    mynum = length(find(matched_filenames2))
    if mynum>0
        locked=0
        return locked
    end

    cmd="touch"
    Base.run(`$cmd $lockname`)
    Base.run(`sync`)
    locked=1
    return locked

end



function releaseFile(filename)

    lockname=filename*".lock";
    fnames = readdir()

    matched_filenames2 = Array{Bool}(length(fnames))
    for i=1:length(fnames)
        matched_filenames2[i] = ismatch(Regex(@sprintf("^%s", lockname)), fnames[i])
    end
    mynum = length(find(matched_filenames2))
    if mynum>0
        cmd="rm"
        arg1="-f"
        Base.run(`$cmd $arg1 $lockname`)
    else
        println("No lock found")
    end
end


