from cheyenne.read_amplxe_results import read_amplxe_results

def read_MOM6_time(logfilepath, clockname):
    logfile = open(logfilepath, 'r')
    for line in logfile:
        if line[0:len(clockname)]==clockname:
            return float(line.split()[-1])
    raise RuntimeError("Couldn't find clock entry in logfile")


def get_MOM_continuity_PPM_time():

    totalTime = 0.0

    routines = read_amplxe_results("./result.txt")

    for r in [routines[r] for r in routines]:
        if r.name.startswith("mom_continuity_ppm"):
            totalTime += float(r.self)

    return totalTime


if __name__ == "__main__":
    # clock_continuity = read_MOM6_time("./logfile.000000.out", "(Ocean continuity PPM)")
    clock_continuity = get_MOM_continuity_PPM_time()
    print(clock_continuity)
