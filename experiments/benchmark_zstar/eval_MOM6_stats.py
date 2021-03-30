from cheyenne.read_amplxe_results import read_amplxe_results
import sys
import os

def get_MOM_continuity_PPM_time():

    totalTime = 0.0

    routines = read_amplxe_results(os.path.join(logPath,"vtune_report.txt"))

    for r in [routines[r] for r in routines]:
        if r.name.startswith("mom_continuity_ppm"):
            totalTime += float(r.self)

    return totalTime


if __name__ == "__main__":
    logPath = sys.argv[1]
    clock_continuity = get_MOM_continuity_PPM_time(logPath)
    print(clock_continuity)
