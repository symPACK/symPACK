from collections import defaultdict
import getopt 
import sys

def analyze():
    
    opts, args = getopt.getopt(sys.argv[1:], "n:m:")
    
    nodes = 0
    outfile_name = ""
    indir_name = ""
    for opt, arg in opts:
        if opt=="-n":
            nodes = int(arg)
        if opt=="-m":
            if arg=="cpu":
                outfile_name = "results_cpu.txt"
                indir_name = "buildt"
            elif arg=="gpu":
                outfile_name = "results_gpu.txt"
                indir_name = "build3"

    update_cpy_dat = defaultdict(lambda: 0)
    update_dat = defaultdict(lambda:0)
    trsm_dat = defaultdict(lambda:0)
    fact_dat = defaultdict(lambda:0)
    fetch_dat = defaultdict(lambda : 0)
    axpy_dat = defaultdict(lambda:0)
    alloc_dat = defaultdict(lambda:0)
    
    outfile = open(outfile_name, "w+")

    for p in range(nodes):
        
        file = open(f"{indir_name}/Statfile{p}", "r+")
        
        for line in file:
            if line.find("Update_cpy")!=-1:
                update_cpy_dat = parse_line(update_cpy_dat, line, p) 
            if line.find("Update: ")!=-1:
                update_dat = parse_line(update_dat, line, p)
            if line.find("TRSM")!=-1:
                trsm_dat = parse_line(trsm_dat, line, p)
            if line.find("Factorize")!=-1:
                fact_dat = parse_line(fact_dat, line, p)
            if line.find("Fetch")!=-1:
                fetch_dat = parse_line(fetch_dat, line, p)
            if line.find("Axpy")!=-1:
                axpy_dat = parse_line(axpy_dat, line, p)
            if line.find("Alloc")!=-1:
                alloc_dat = parse_line(alloc_dat, line, p)
        
        outfile.write(f"***********Stats for P{p}***********\n")
        outfile.write(f"Update_cpy total: {update_cpy_dat[p]}\n")  
        outfile.write(f"Update total: {update_dat[p]}\n")
        outfile.write(f"TRSM total: {trsm_dat[p]}\n")
        outfile.write(f"Factorize total: {fact_dat[p]}\n")
        outfile.write(f"Fetch total: {fetch_dat[p]}\n")
        outfile.write(f"Axpy total: {axpy_dat[p]}\n")
        outfile.write(f"Alloc total: {alloc_dat[p]}\n")
                        

def parse_line(data, line, p):
    time = line.split(": ")[1]
    time = int(time.strip("ns\n"))        
    data[p] += (time / 1e9)
    return data



if __name__=="__main__":
    analyze()
