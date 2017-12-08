# DON'T MODIFY THIS FILE -- the source is in file Scraper.ipynb. Look there for further documentation and examples of running the code.


using JSON


# DON'T MODIFY THIS FILE -- the source is in file Scraper.ipynb. Look there for further documentation and examples of running the code.


# Get database of scraped files

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
            get!(answer, lstrip(A[i,2]), A[i,1])
        end
    catch
        answer = Dict()
    end
    return answer
end


"""
write_scrapedict(latest; scrapedir=".scrapedir", scrapefile="scrapelist")

Writes a dictionary containing filename-latest_scrape_time_string pairs into a file

"""
function write_scrapedict(latest; scrapedir=".scrapedir", scrapefile="scrapelist")

    sfile = scrapedir * "/" * scrapefile
    sf = open(sfile, "w")
    for k in keys(latest)
        write(sf, @sprintf("%s, %s\n", latest[k], k))
    end
    close(sf)

end


# DON'T MODIFY THIS FILE -- the source is in file Scraper.ipynb. Look there for further documentation and examples of running the code.


# Scrape a notebook for julia code that should be written into an indicated file

"""
filenames_written = scrape_notebook(notebook_filename; verbose=false, includemagic="#@include_me",
    update_db=false)

Goes through a file notebook_filename, assuming it is an ipynb, and looks for code cells that start with
includemagic, followed by whitespace, followed by a string (which we shall call filename). 
When such a code cell is found, its contents are written into filename.

If more than one cell uses the same filename, then the first one starts the file, and subsequent cells
append to it.

Returns an array with the written filenames

"""
function scrape_notebook(notebook_filename; verbose=false, includemagic="#@include_me", update_db=false)

    filenames = [];    # List of output files found in this notebook
    A = []  # declare A outside the try/catch so it will be available as a variable outside the try/catch
    try 
        A = JSON.parse(readstring(notebook_filename))
    catch y
        @printf("\n=======\n\n   WARNING!!! Ran into trouble trying to JSON parse file %s\n\n", notebook_filename)
        @printf("Error was "); print(y); print("\n\n======\n")
        return filenames
    end
    
    
    for mycell in A["cells"]
        if mycell["cell_type"] == "code"   
            lines = mycell["source"]
            if length(lines)>0           # We only consider code cells that are not empty
                m= match(Regex(@sprintf("(?<include>%s)\\s*(?<filename>\\S*)", includemagic)), lines[1])
                if typeof(m)!=Void && length(m["filename"])>0  # proceed if we got a match and got a filename
                    if any(filenames .== m["filename"]); 
                        f = open(m["filename"], "a")           # we'll append if we already had that filename
                        if verbose; @printf("Appending to file %s\n", m["filename"]); end
                    else
                        f = open(m["filename"], "w")           # otherwise open fresh for writing
                        filenames = [filenames ; m["filename"]]
                        if verbose; @printf("Writing out file %s\n", m["filename"]); end
                    end
                    # Now write out the contents of the cell, with a warning at the top:
                    write(f, @sprintf("# DON'T MODIFY THIS FILE -- the source is in file %s. Look there for further documentation and examples of running the code.\n\n", notebook_filename))
                    for i=2:length(lines)
                        write(f, lines[i])
                    end
                    write(f, "\n\n\n")
                    close(f)
                end
            end
        end
    end
    
    if update_db && length(filenames)>0
        latest = latest_scrapedict();
        if haskey(latest, notebook_filename)
            latest[notebook_filename] = string(now())
        else
            get!(latest, notebook_filename, string(now()))
        end

        if verbose; @printf("Refreshing database with info about %s\n", notebook_filename) end;
        write_scrapedict(latest)
    end

    return filenames
end


# DON'T MODIFY THIS FILE -- the source is in file Scraper.ipynb. Look there for further documentation and examples of running the code.


# Go through all notebooks in directory and scrape them if they've been modified after their last
# scrape time.

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

    if length(rescraped)>0
        write_scrapedict(latest; scrapedir=scrapedir, scrapefile=scrapefile)
    end
    
    return rescraped
end




# DON'T MODIFY THIS FILE -- the source is in file Scraper.ipynb. Look there for further documentation and examples of running the code.


function scraperobot()
    while true
        scrape_all_notebooks()
        sleep(2)
    end
end


if length(ARGS)==0
    scraperobot()
elseif any(ARGS[1] .== ["-h", "--h", "-help", "--help"])
        @printf("\nUsage: julia scraper.jl &     to run the scraping robot in the background\n\n")
        @printf("   OR\n\n")
        @printf("julia scraper.jl notebook1.ipynb [notebook2.ipynb ...]   to scrape the indicated notebooks, then stop.\n\n")
else
    for nf in ARGS
        if endswith(nf, ".ipynb") && isfile(nf)
            scrape_notebook(nf; verbose=true, update_db=true)
        end
    end
end



