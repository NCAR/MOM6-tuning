
def read_MOM6_time(logfilepath, clockname):
    logfile = open(logfilepath, 'r')
    for line in logfile:
        if line[0:len(clockname)]==clockname:
            return float(line.split()[-1])
    raise RuntimeError("Couldn't find clock entry in logfile")
    

if __name__ == "__main__":
    clock_continuity = read_MOM6_time("./logfile.000000.out", "(Ocean continuity PPM)")
    print(clock_continuity)
