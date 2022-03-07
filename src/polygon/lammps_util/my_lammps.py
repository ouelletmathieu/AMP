import random

class ID_counter:
    """simple class to count lammps call id
    """

    def __init__(self):
        self.active = 0   

    def next(self):
        self.active+=1
        return str(self.active)

class MyLammps:
    """class that handles call to Lammps  
    """
    
    def __init__(self, lammps, id_counter):
        self.id_counter = id_counter
        self.L = lammps
        
    #
    def create_molecule_2d(self,file_path):
        """read a Lammps file and create a molecule:
            file_path = path of the lammps file
        """
        self.L.units("lj")
        self.L.command("dimension 2")
        self.L.command("atom_style bond")
        self.L.command("boundary f f p")
        self.L.command("bond_style harmonic")
        self.L.command("read_data "+ file_path) 
        self.L.command("fix " + self.id_counter.next() + " all enforce2d")

    def create_molecule_3d(self,file_path):
        """read a Lammps file and create a molecule:
            file_path = path of the lammps file
        """
        self.L.units("lj")
        self.L.command("dimension 3")
        self.L.command("atom_style bond")
        self.L.command("boundary f f f")
        self.L.command("bond_style harmonic")
        self.L.command("read_data "+ file_path) 
        
    
    
    def create_molecule_2d_with_lj(self, cut_off, dist_emin, emin, file_path):
        """read a Lammps file and create a molecule with a lennard Jones potential of parameter:

        Args:
            cut_off ([type]): what distance to stop considering the potential
            dist_emin ([type]): distance to the minima  
            emin ([type]):energy at the minima
            file_path ([type])
        """       
        #TODO remove unused functions 
        print("cut_off, dist_emin, emin not implemented () lenard jones tranfered to another constructor")
        self.L.units("lj")
        self.L.command("dimension 2")
        self.L.command("atom_style full")
        self.L.command("boundary f f p")
        self.L.command("neighbor 1.0 nsq")
        self.L.command("bond_style harmonic")
        self.L.command("pair_style lj/cut 2.5")

        
        self.L.command("read_data "+ file_path) 
        self.L.command("fix " + self.id_counter.next() + " all enforce2d")
        self.L.command("fix " + self.id_counter.next() + " all nve")      
        #sigma = 0.89*dist_emin
        self.L.command("pair_coeff * * 1.0 1.0  2.5" )
        

    def run_minimize(self, max_iter, logger):
        """minimize the energy of the molecule. The initial molecule is stored in the logger at time 0 and the optimize one at time 1

        Args:
            max_iter ([int]): number of iteration to execute 
            logger ([Logger]):  logger class
        """
        
        logger.log(0)
        self.L.command("fix " + self.id_counter.next() + " all   nve")   ###NVE ensemble
        self.L.command("fix " + self.id_counter.next() + " all enforce2d")
        self.L.command("thermo " + str(100))
        self.L.command("minimize 0.000001 0.000001 "+ str(max_iter) + "  " + str(100*max_iter)  )
        
        logger.log(1)

    def run_brownian(self, temperature, damping, max_time, n_step, n_plot, logger, stop_change_conf = False, id_conf =[], square_angle_sum = 0, random_seed = True ):
        """run langevin dynamic 

        Args:
            temperature ([float]): temperature of the bath (see lammps langevin dynamic documentation)
            damping ([float]): damping of the bath (see lammps langevin dynamic documentation)    
            max_time ([type]): [description]
            n_step ([int]): number of step for integration for all the time 
            n_plot ([int]): number of time the position, angle etc are logged 
            logger ([type]): logger class that will contain the info (the output)
            stop_change_conf (bool, optional): if yes it will stop if there is a change of conformation detected at one of the plot time
            id_conf (list, optional): either  ["healthy", "prion"] if stop_change_conf is true. Defaults to [].
            square_angle_sum (int, optional): square_angle_sum min angular distance to be considered a change of conformation . Defaults to 0.
            random_seed (bool, optional): can put a given seed for testing and reproducibility Defaults to True.

        Returns:
            [float]: return the max time the simulation reached if stop_change_conf is true it return the time 
        """
             
        self.L.command("fix " + self.id_counter.next() + " all   nve")   ###NVE ensemble
        self.L.command("fix " + self.id_counter.next() + " all enforce2d")
        
        if random_seed:
            self.L.command("fix " + self.id_counter.next() + " all langevin "+ str(temperature) + " " + str(temperature)+ " " + str(damping)+ " " + str(random.randint(1,100000)))
        else:
            self.L.command("fix " + self.id_counter.next() + " all langevin "+ str(temperature) + " " + str(temperature)+ " " + str(damping)+ " " + str(1) )
               
        self.L.command("timestep " + str(max_time/n_step))

        step_per_plot = int(n_step/n_plot) 
        self.L.command("thermo " + str(int(step_per_plot)))
        
        for i in range(n_plot):

            try:
                self.L.command("run " + str(step_per_plot))
                time = i*step_per_plot*max_time/n_step
                logger.log(time)
                
                if stop_change_conf:
                    angle_distance = logger.get_last_distance( -1, id_to_compare = id_conf)
                    if angle_distance > square_angle_sum:
                        return logger.time_list[-1]

            except Exception as inst:
                print(type(inst))    # the exception instance
                print(inst.args)     # arguments stored in .args
                print(inst)
                return -1
        
        return max_time
        
    def command(self,command):
        """run a command in lammps 
        """
        self.L.command(command)

    def getEnergy(self):
        """get the potential energy of the system at the current time 
        """
        self.L.run(0)
        e = self.L.eval("pe")
        return e
    
    def getEnergyKin(self):
        """get the kinetic energy of the system at the current time   
        """
        self.L.run(0)
        e = self.L.eval("ke")
        return e    

    def lammps(self):
        """get the lammps instance
        """
        return self.L