
using JSON

"""
filenames_written = scrape_notebook(notebook_filename; verbose=false, includemagic="#@include_me")

Goes through a file notebook_filename, assuming it is an ipynb, and looks for code cells that start with
includemagic, followed by whitespace, followed by a string (which we shall call filename). 
When such a code cell is found, its contents are written into filename.

Returns an array with the written filenames

"""
function scrape_notebook(notebook_filename; verbose=false, includemagic="#@include_me")

    filenames = [];    
    A = JSON.parse(readstring(notebook_filename))
    for mycell in A["cells"]
        if mycell["cell_type"] == "code"
            lines = mycell["source"]
            if length(lines)>0
                m= match(Regex(@sprintf("(?<include>%s)\\s*(?<filename>\\S*)", includemagic)), lines[1])
                if typeof(m)!=Void && length(m["filename"])>0
                    if verbose; @printf("Writing out file %s\n", m["filename"]); end
                    filenames = [filenames ; m["filename"]]
                    f = open(m["filename"], "w")
                    for i=2:length(lines)
                        write(f, lines[i])
                    end
                    close(f)
                end
            end
        end
    end

    return filenames
end



"""
latest = latest_scrapedict(; scrapedir=".scrapedir", scrapefile="scrapelist")

Returns a dictionary that maps filenames to strings representing when they were
last scraped.  This information is stored in a human-readable text file, scrapedir/scrapefile

"""
function latest_scrapedict(; scrapedir=".scrapedir", scrapefile="scrapelist")
    if !isdir(scrapedir); mkdir(scrapedir); end;
    sfile = scrapedir * "/" * scrapefile
    if !isfile(sfile);
        return Dict()
    end
    
    answer = Dict()
    try
        A = readdlm(sfile, ',')
        for i=1:size(A,1)
            get!(answer, A[i,1], lstrip(A[i,2]))
        end
    catch
        answer = Dict()
    end
    return answer
end



"""
rescraped = scrape_all_notebooks(; scrapedir=".scrapedir", scrapefile="scrapelist", verbose=false)


"""
function scrape_all_notebooks(; scrapedir=".scrapedir", scrapefile="scrapelist", verbose=false)

    latest = latest_scrapedict(scrapedir=scrapedir, scrapefile=scrapefile)

    rescraped = []
    for f in filter(x -> endswith(x, ".ipynb"), readdir())
        if ~haskey(latest, f) || DateTime(latest[f]) < Dates.unix2datetime(stat(f).mtime) - Dates.Hour(4)
            if verbose; @printf("Will look into notebook %s\n", f); end
            rescraped = [rescraped; f]
            scrape_notebook(f)
            if haskey(latest, f)
                latest[f] = string(now())
            else
                get!(latest, f, string(now()))
            end
        end
    end

    sfile = scrapedir * "/" * scrapefile
    sf = open(sfile, "w")
    for k in keys(latest)
        write(sf, @sprintf("%s, %s\n", k, latest[k]))
    end
    close(sf)
    
    return rescraped
end

scrape_all_notebooks()
