def main():
    import util

    import subprocess
    import sys


    result = subprocess.run(
        ["sh", "/Users/mathieuouellet/Desktop/prions/script_neb"], capture_output=True, text=True
        #['mpirun', '-np', '12', '--oversubscribe', 'lmp_mpi', '-partition', '12x1', '-in', '/Users/mathieuouellet/Desktop/prions/neb_task.lj'], capture_output=True, text=True
        #["mpirun", '-h'], capture_output=True, text=True
        #['mpirun', '-np', '4', 'lmp_mpi', '-partition', '4x1','-in', '/Users/mathieuouellet/lammps-29Oct20/examples/neb/in.neb.hop1'], capture_output=True, text=True, shell=True, check=True
        #['mpirun -np 4 pwd'], capture_output=True, text=True, shell=True
    )
    print(result)


    print("done")

if __name__ == "__main__":
    main()

    