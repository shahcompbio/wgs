import subprocess
import sys

def is_error(err):
    if len(err): return True

def got_output(out):
    if len(out): return True

def failure(err): 
    print ("couldn't complete testing, got error message: ", "\n\n\n", err)
    
def success(test):
    print ("test \"{}\" completed sucessfully".format(str(test)))

def handle_process_returns(out, err, process):
    if is_error(err): 
        failure(err)
        return True, 

    if not is_error(err) and got_output(out): 
        success(process)
        return False

def format_cmd(cmd): return " ".join(cmd)

def run_cmd(cmd):
    process = subprocess.Popen(cmd,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE)
    return process.communicate()

def extract_header_line(line):
    
    line = [ element.split(':') for element in line[1:] ]
    return {kvp[0]:kvp[1] for kvp in line}

def check_header_formation(bam, compare_bam = None, alignment=None, ref=None):
    
    cmd = ['samtools', 'view', '-H', bam]
    out, err = run_cmd(cmd)
    is_failure = handle_process_returns(out, err, format_cmd(cmd))
    if  is_failure: return
    bam_header = out

    if "@RG" in bam_header and "@PG" in bam_header: 
        readgroups = ""
        alignment_cmd = ""

        for line in bam_header.split("\n"):
            if line.startswith("@RG"): readgroups = line.split()
            if "bwa" in line: alignment_cmd = line

        if len(readgroups) <= 1 and len(alignment_cmd) <= 1: 
            failure("bam does not have read groups \
             and alignment command defined correctly in header \
             readgroups: {}".format(readgroups) + 
             "alignment_cmd: {}".format(alignment_cmd))
            return

        #check readgroups in header
        readgroups = extract_header_line(readgroups)
        if "SM" in readgroups and "ID" in readgroups and "CN" in readgroups:
            success("SM, ID and CN in header RG")
        else: failure("readgroups not defined correctly in bam header"); return
        
        #check alignment cmd and reference in header
        if not alignment: alignment = "bwa mem"
        if not ref: ref = ".fa"

        if alignment in alignment_cmd and ref in alignment_cmd: 
            success("header contains correct reference genome and alignment command")
        else:
            failure("header does not contain correct reference genome and alignment command")
            return
    else:
        failure("bam does not have read groups and alignment command defined in header")


    if compare_bam:
        cmd = ['samtools', 'view', '-H', compare_bam]
        out, err = run_cmd(cmd)
        is_failure = handle_process_returns(out, err, format_cmd(cmd) )
        if is_failure: return False
        compare_bam_header = out

        if bam_header == compare_bam_header: 
            success("bam header is the same as the comparison bam")            
        else: 
            failure("the two bams differ")
    else:
        success("bam header check")
        return True

def check_bam_body_formation(bam, compare_bam = None):
    cmd = ['samtools', 'view', '-c', '-F', '260', bam]
    out, err = run_cmd(cmd)
    is_failure = handle_process_returns(out, err, format_cmd(cmd) )
    if is_failure: return False
    
    bam_body = out
    
    if bam_body == 0: 
        failure("bam has no aligned reads")
        return False

    if compare_bam:
        cmd = ['samtools', 'view', '-c', '-F', '260', bam]
        out, err = run_cmd(cmd)
        is_failure = handle_process_returns(out, err, format_cmd(cmd))
        if is_failure: return False
        compare_bam_body = out

        if bam_body != compare_bam_body:
            failure("two bams have different # mapped reads")
            return  False
        else:
            success("two bams have sam # mapped reads")
    
    success("bam body check")
    return True

is_good_header = check_header_formation(sys.argv[1], sys.argv[2])
is_good_bod = check_bam_body_formation(sys.argv[1], sys.argv[2])

if is_good_bod and is_good_header:
    print ("this bam checks out!" )