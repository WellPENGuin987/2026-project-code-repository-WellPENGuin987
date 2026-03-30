import numpy as np


def ChooseFile():
    # placeholder for file input, not implemented yet
    pass


def SetUnits(Units):
    allowed_units_lwr = ["nm", "mim", "mm", "cm", "m"]
    if Units.lower() in allowed_units_lwr:
        Units = Units.lower()
    elif Units == "a":
        Units = "A"
    else:
        raise ValueError
    return Units


def SetTime(Time):
    allowed_units_lwr = ["ns", "mis", "ms", "s"]
    if Time.lower() in allowed_units_lwr:
        Time = Time.lower()
    else:
        raise ValueError
    return Time


def check_positive_float(parameter, values):
    try:
        values = [np.longdouble(v) for v in values]
        if any(v <= 0 for v in values):
            raise ValueError(f"{parameter} must be positive")
        return values
    except (ValueError, TypeError) as e:
        print("Invalid input: ", e)
        raise


def setMomentsOfInertia(degrees):
    degrees = int(degrees)
    if degrees == 3:
        d, MoI, vib = Moments_of_inertia_point()
    elif degrees == 4:
        d, MoI, vib = Moments_of_inertia_linear()
    elif degrees == 6:
        d, MoI, vib = Moments_of_inertia_non_linear()
    elif degrees > 6:
        d, MoI, vib = Moments_of_inertia_vibrational(degrees)
    else:
        raise ValueError("Invalid number of degrees of freedom, must be at least 3")

    return {
        "degrees": d,
        "moments_of_inertia": MoI,
        "vibrational_modes": vib,
    }


def Moments_of_inertia_point():
    return 3, [0, 0, 0], []  # for point-like particles, moments of inertia are all zero


def Moments_of_inertia_linear():
    while True:
        moments_of_inertia = input("Moment of inertia (I) in kg*m^2\t format: I_x,I_y (I_z is the minimum of I_x & I_y):\n").replace(" ", "").split(",")
        if len(moments_of_inertia) != 2:
            print("Invalid number of moments of inertia, expected 2 for symmetric top particles, got {len(moments_of_inertia)}")
        else:
            moments_of_inertia = [
                np.longdouble(str(moments_of_inertia[0])),
                np.longdouble(str(moments_of_inertia[1])),
                min(np.longdouble(str(moments_of_inertia[0])), np.longdouble(str(moments_of_inertia[1]))),
            ]  # for symmetric top particles, the 3rd moment is the min of the first two
            check_positive_float("Moments of inertia", moments_of_inertia)
            break
    return 5, moments_of_inertia, []


def Moments_of_inertia_non_linear():
    while True:
        moments_of_inertia = input("Moments of inertia (I) in kg*m^2\t format: I_x,I_y,I_z:\n").replace(" ", "").split(",")
        if len(moments_of_inertia) != 3:
            print("Invalid number of moments of inertia, expected 3 for rotational degrees of freedom, got {len(moments_of_inertia)}")
        else:
            check_positive_float("Moments of inertia", moments_of_inertia)
            break
    return 6, moments_of_inertia, []


def Moments_of_inertia_vibrational(degrees):
    while True:
        moments_of_inertia = input("Moments of inertia (I) in kg*m^2\t format: I_x,I_y,I_z:\n").replace(" ", "").split(",")
        if len(moments_of_inertia) != 3:
            print("Invalid number of moments of inertia, expected 3 for rotational degrees of freedom, got {len(moments_of_inertia)}")
        else:
            check_positive_float("Moments of inertia", moments_of_inertia)
            break
    while True:
        vibrational_modes = input().replace(" ", "").split(",")
        if len(vibrational_modes) != degrees - 6:
            print(f"Invalid number of vibrational modes, expected {degrees - 6}, got {len(vibrational_modes)}")
        else:
            check_positive_float("Vibrational modes", vibrational_modes)
            break
    return degrees, moments_of_inertia, vibrational_modes


def setPosDistType(pos_dist_type):
    allowed_pos_dist_types = [1, 2, 3]
    if pos_dist_type in allowed_pos_dist_types:
        return pos_dist_type
    else:
        raise ValueError


def setVelDistType(vel_dist_type):
    allowed_vel_dist_types = [1, 2, 3]
    if vel_dist_type in allowed_vel_dist_types:
        return vel_dist_type
    else:
        raise ValueError


def InitialiseBox():
    while True:
        Init_Box = (
            input(
                "Input parameters for the box:\n"
                "Len_X/Y/Z are the dimensions of the container\n"
                "UnitsX/Y/Z: supports A (Ångstroms), nm, mim (micrometres), mm, cm, m\n"
                "Format:\tLen_X,Units_X,Len_Y,Units_Y,Len_Z,UnitsZ\n"
            )
            .replace(" ", "")
            .split(",")
        )
        print(Init_Box)

        if Init_Box == ["F"] or Init_Box == ["f"]:
            # read paramaters from file, not implemented yet
            pass
        elif len(Init_Box) != 6:
            print("Invalid number of dimensions and units")
        else:
            try:
                Len_X, Units_X, Len_Y, Units_Y, Len_Z, Units_Z = Init_Box
                Len_X = np.longdouble(str(Len_X))
                Len_Y = np.longdouble(str(Len_Y))
                Len_Z = np.longdouble(str(Len_Z))
                Units_X = SetUnits(Units_X)
                Units_Y = SetUnits(Units_Y)
                Units_Z = SetUnits(Units_Z)
                Box_Params = [Len_X, Units_X, Len_Y, Units_Y, Len_Z, Units_Z]
                return Box_Params
            except (ValueError, TypeError):
                print("Invalid inputs")


def InitialiseTime():
    while True:
        Init_Time = (
            input(
                "Input paramenters for time:\n"
                "Time is the total run time\tdt is the time step between each test\n"
                "Units_Time/ dt/ T_plt: supports ns, mis (microseconds), ms, s\n"
                "T_plt is the time between data samples for plots\n"
                "Format:\tTime,Units_Time,dt,Units_dt,T_plt,Units_T_plt\n"
            )
            .replace(" ", "")
            .split(",")
        )

        if Init_Time == ["F"] or Init_Time == ["f"]:
            # read paramaters from file, not implemented yet
            pass
        elif len(Init_Time) != 6:
            print("Invalid number of dimensions and units")
        else:
            try:
                Time, Units_Time, dt, Units_dt, T_plt, Units_T_plt = Init_Time
                Time = np.longdouble(str(Time))
                dt = np.longdouble(str(dt))
                T_plt = np.longdouble(str(T_plt))
                Units_Time = SetTime(Units_Time)
                Units_dt = SetTime(Units_dt)
                Units_T_plt = SetTime(Units_T_plt)
                Time_Params = [Time, Units_Time, dt, Units_dt, T_plt, Units_T_plt]
                return Time_Params
            except (ValueError, TypeError):
                print("Invalid inputs")


def InitialiseParticles():
    Sep_Particles = []
    while True:
        Init_Particles = (
            input(
                "Input parametetrs of particles (X to escape):\n"
                "type is the type of particle, used to differentiate different types of particles in the graphs\n"
                "N is the number of particles of this type\tm is the mass of the particles (kg)\n"
                "D is how many degrees of freedom these particles have: 3 translational (minimum) + 1 or 3 rotational + 1 vibrational\n"
                "r is the radius of collision of these particles (Å)\tTemp is the initial temperature of these particles (K)\n"
                "Pos_dist is the type of distribution of the initial positions of these particles (Lattice = 1, Uniform pdf = 2, Custom gradient pdf = 3)\n"
                "Vel_dist is the type of distribution of the initial velocities of these particles (Identical RMS = 1, Uniform RMS pdf = 2, Maxwell-Boltzmann = 3)\n"
                "Format:\ttype,N,m,D,r,Temp,Pos_dist,Vel_dist\n"
            )
            .replace(" ", "")
            .split(",")
        )

        if Init_Particles == ["X"] or Init_Particles == ["x"]:
            print("Initialisation complete")
            break
        elif Init_Particles == ["F"] or Init_Particles == ["f"]:
            # read paramaters from file, not implemented yet
            pass
        else:
            print(f"{len(Init_Particles)} parameters inputted:\n{Init_Particles}")
        if len(Init_Particles) % 8 == 0:
            Current_Particles = []  # list to store parameters of current particle type, will be added to list of all particle parameters
            try:
                Current_Particles.append(str(Init_Particles[0]))  # particle type
                Current_Particles.append(int(Init_Particles[1]))  # number of particles of this type
                Current_Particles.append(np.longdouble(str(Init_Particles[2])))  # mass of particles of this type

                # keep degrees as numeric for simulation, and keep inertia structure as extra metadata
                inertia_data = setMomentsOfInertia(Init_Particles[3])
                Current_Particles.append(int(Init_Particles[3]))  # degrees of freedom
                Current_Particles.append(inertia_data)  # moments of inertia / vibrational modes metadata

                Current_Particles.append(np.longdouble(str(Init_Particles[4])))  # radius of collision
                Current_Particles.append(np.longdouble(str(Init_Particles[5])))  # initial temperature
                Current_Particles.append(setPosDistType(int(Init_Particles[6])))  # position distribution type
                Current_Particles.append(setVelDistType(int(Init_Particles[7])))  # velocity distribution type

                Sep_Particles.append(Current_Particles)  # add current particle parameters to list of all particle parameters

            except (ValueError, TypeError) as e:
                print("Invalid inputs: ", e)
        else:
            print("Invalid number of particle parameters")
    return Sep_Particles


def Initialise():
    # --------------initialising the parameters of the simulation--------------#
    Box_Params = InitialiseBox()
    Time_Params = InitialiseTime()
    Sep_Particles = InitialiseParticles()
    return (Box_Params, Time_Params, Sep_Particles)
